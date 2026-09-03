#pragma once

#include "../simulations/ISimulation.hpp"
#include "../types/ParticleCollection.hpp"
#include "../types/FluidCube.hpp"
#include "../types/MacGrid.hpp"
#include "../types/Scene.hpp"
#include "../Macros.hpp"
#include "../Maths.hpp"

namespace fsim {
    class FlipSimulation2d : public ISimulation {
    public:
        FlipSimulation2d(FluidCube *grid, ParticleCollection *particles, const Scene &scene) {
            this->grid = grid;
            this->particles = particles;
            this->scene = scene;
            this->numCellParticles = nullptr;
            this->firstCellParticles = nullptr;
            this->cellParticleIds = nullptr;
        }

        ~FlipSimulation2d() {
            if (grid) delete grid;
            if (particles) delete particles;
            if (numCellParticles) deleteIntArray(numCellParticles);
            if (firstCellParticles) deleteIntArray(firstCellParticles);
            if (cellParticleIds) deleteIntArray(cellParticleIds);
        }

        virtual void step() override {
            loopParticles(grid->time_step);
            transferVelocitiesToGrid();
            loopGrid(grid->time_step);
            transferVelocitiesToParticles();
        }

        FluidCube *getGrid() {
            return grid;
        }
        ParticleCollection *getParticles() {
            return particles;
        }

    protected:
        Scene scene;
        FluidCube *grid;
        ParticleCollection *particles;

        int *numCellParticles, *firstCellParticles, *cellParticleIds;

        virtual void loopParticles(float time) {
            moveParticles(time);
            pushParticles(time);
        }

        virtual void transferVelocitiesToGrid() {

        }

        virtual void transferVelocitiesToParticles() {

        }

        virtual void loopGrid(float time) {
            solvePressure(time);
        }

        virtual void solvePressure(float time) {
            MacGrid *mac = (MacGrid*) grid;
            int N = mac->size;
            float cp = mac->density * mac->h / time;

            for (int i = 0; i < N*N; i++) mac->p[i] = 0.0f;

            runForNSteps((int)scene.iterations) {
                for (int i = 1; i < N-1; i++) {
                    for (int j = 1; j < N-1; j++) {
                        if (!mac->isFluid(i, j)) continue;

                        float sx0 = mac->solidFlag(i-1, j);
                        float sx1 = mac->solidFlag(i+1, j);
                        float sy0 = mac->solidFlag(i, j-1);
                        float sy1 = mac->solidFlag(i, j+1);
                        float s   = sx0 + sx1 + sy0 + sy1;
                        if (s == 0.0f) continue;

                        float div = mac->Vx[mac->indexOf(i+1, j)] - mac->Vx[mac->indexOf(i, j)]
                                  + mac->Vy[mac->indexOf(i, j+1)] - mac->Vy[mac->indexOf(i, j)];

                        float p = -div / s * scene.overRelaxation;
                        mac->p[mac->indexOf(i, j)] += cp * p;

                        mac->Vx[mac->indexOf(i,   j)] -= sx0 * p;
                        mac->Vx[mac->indexOf(i+1, j)] += sx1 * p;
                        mac->Vy[mac->indexOf(i, j  )] -= sy0 * p;
                        mac->Vy[mac->indexOf(i, j+1)] += sy1 * p;
                    }
                }
            }
        }

        virtual void moveParticles(float time) {
            for (int p=0; p < particles->size; p++) {
                particles->vy[p] -= scene.gravity * time;
                particles->px[p] += particles->vx[p] * time;
                particles->py[p] += particles->vy[p] * time;
            }
        }

        virtual void pushParticles(float time) {
            float pInvSpacing = 1.0f / (2.2f * scene.particleRadius);
            int pNumX = floor(scene.width * pInvSpacing) + 1;
            int pNumY = floor(scene.height * pInvSpacing) + 1;
            int numCells = pNumX * pNumY;
            
            if (!numCellParticles) numCellParticles = makeIntArray(numCells); // number of particles per cell
            if (!firstCellParticles) firstCellParticles = makeIntArray(numCells+1); // partial sums
            if (!cellParticleIds) cellParticleIds = makeIntArray(particles->size);  // from grid to particle

            for (int i=0; i < particles->size; i++) {
                numCellParticles[i] = 0;
            }

            for (int i=0; i < particles->size; i++) {
                float x = particles->px[i];
                float y = particles->py[i];

                int cellX = fsim::clamp((int) floor(x * pInvSpacing), 0, pNumX - 1);
                int cellY = fsim::clamp((int) floor(y * pInvSpacing), 0, pNumY - 1);
                int cellNum = cellX * pNumY + cellY;
                numCellParticles[cellNum] ++;
            }

            int partialSum = 0;
            for (int i=0; i < numCells; i++) {
                partialSum += numCellParticles[i];
                firstCellParticles[i] = partialSum;
            }
            firstCellParticles[numCells] = partialSum;

            for (int i=0; i < particles->size; i++) {
                float x = particles->px[i];
                float y = particles->py[i];

                int cellX = fsim::clamp((int) floor(x * pInvSpacing), 0, pNumX - 1);
                int cellY = fsim::clamp((int) floor(y * pInvSpacing), 0, pNumY - 1);
                int cellNum = cellX * pNumY + cellY;

                firstCellParticles[cellNum] --;
                cellParticleIds[firstCellParticles[cellNum]] = i;
            }

            float minDist = 2.0f * scene.particleRadius;
            float minDistSq = minDist * minDist;

            runForNSteps(scene.pushParticlesIter) {
                for (int p=0; p < particles->size; p++) {
                    float px = particles->px[p];
                    float py = particles->py[p];
                    
                    int onGridX = (int) floor(px * pInvSpacing);
                    int onGridY = (int) floor(py * pInvSpacing);

                    float x0 = fsim::Max(onGridX-1, 0);
                    float y0 = fsim::Max(onGridY-1, 0);
                    float x1 = fsim::Min(onGridX+1, pNumX-1);
                    float y1 = fsim::Min(onGridY+1, pNumY-1);

                    for (int x=x0; x<=x1; x++) {
                        for (int y=y0; y<=y1; y++) {
                            int cellNum = x0 * pNumY + y;
                            int first = firstCellParticles[cellNum];
                            int last = firstCellParticles[cellNum+1];

                            for (int i=first; i<last; i++) {
                                int particleId = cellParticleIds[i];
                                if (particleId == p) continue;
                                float dx = particles->px[particleId] - px;
                                float dy = particles->py[particleId] - py;
                                float distSq = dx * dx + dy * dy;
                                if (distSq == 0 || distSq > minDistSq) continue;

                                float dist = fsim::Sqrt(distSq);
                                float scale = 0.5f * (minDist - dist) / dist;
                                dx *= scale;
                                dy *= scale;
                                particles->px[p] -= dx;
                                particles->py[p] -= dy;
                                particles->px[particleId] += dx;
                                particles->py[particleId] += dy;
                            }
                        }
                    }
                    
                }
            }
        }

    };
};
