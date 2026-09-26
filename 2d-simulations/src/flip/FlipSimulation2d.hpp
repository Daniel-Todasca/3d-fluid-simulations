#pragma once

#include "../shared/ISimulation.hpp"
#include "ParticleCollection.hpp"
#include "../stable-fluids/FluidCube.hpp"
#include "MacGrid.hpp"
#include "../shared/Scene.hpp"
#include "../defines.hpp"
#include "../shared/Maths.hpp"

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

            if (particleDensity) deleteFloatArray(particleDensity);
            if (du) deleteFloatArray(du);
            if (dv) deleteFloatArray(dv);
        }

        virtual void step() override {
            if (!particleDensity) particleDensity = makeFloatArray(((int)grid->volume()));
            if (!du) du = makeFloatArray(((int)grid->volume()));
            if (!dv) dv = makeFloatArray(((int)grid->volume()));

            float timeSpeedUp = 10;

            applyGravity(grid->time_step * timeSpeedUp);
            transferVelocitiesToGrid();
            loopGrid(grid->time_step * timeSpeedUp);
            transferVelocitiesToParticles();
            loopParticles(grid->time_step * timeSpeedUp);
        }

        FluidCube *getGrid() {
            return grid;
        }
        ParticleCollection *getParticles() {
            return particles;
        }
        const Scene &getScene() {
            return scene;
        }

    protected:
        Scene scene;
        FluidCube *grid;
        ParticleCollection *particles;

        int *numCellParticles, *firstCellParticles, *cellParticleIds, *cellType;
        float particleRestDensity = 0;
        float *particleDensity = nullptr;
        float *du = nullptr;
        float *dv = nullptr;

        virtual void applyGravity(float time) {
            float gravityBoost = 50;
            for (int p=0; p < particles->size; p++) {
                particles->vy[p] -= scene.gravity * gravityBoost * time;
            }
        }

        virtual void loopParticles(float time) {
            moveParticles(time);
            pushParticles(time);
            handleParticleCollisions(time);
        }

        virtual void loopGrid(float time) {
            updateParticleDensity();
            solvePressure(time);
        }

        virtual void moveParticles(float time) {
            for (int p=0; p < particles->size; p++) {
                particles->px[p] += particles->vx[p] * time;
                particles->py[p] += particles->vy[p] * time;
            }
        }

        virtual void pushParticles(float time) {
            // time is not used
            float pInvSpacing = 1.0f / (2.2f * scene.particleRadius);
            int pNumX = floor(scene.width * pInvSpacing) + 1;
            int pNumY = floor(scene.height * pInvSpacing) + 1;
            int numCells = pNumX * pNumY;
            
            if (!numCellParticles) numCellParticles = makeIntArray(numCells); // number of particles per cell
            if (!firstCellParticles) firstCellParticles = makeIntArray(numCells+1); // partial sums
            if (!cellParticleIds) cellParticleIds = makeIntArray(particles->size);  // from grid to particle

            for (int i=0; i < numCells; i++) {
                numCellParticles[i] = 0;
            }

            for (int i=0; i < particles->size; i++) {
                float x = particles->px[i];
                float y = particles->py[i];

                int cellX = fsim::Clamp((int) floor(x * pInvSpacing), 0, pNumX - 1);
                int cellY = fsim::Clamp((int) floor(y * pInvSpacing), 0, pNumY - 1);
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

                int cellX = fsim::Clamp((int) floor(x * pInvSpacing), 0, pNumX - 1);
                int cellY = fsim::Clamp((int) floor(y * pInvSpacing), 0, pNumY - 1);
                int cellNum = cellX * pNumY + cellY;

                firstCellParticles[cellNum] --;
                cellParticleIds[firstCellParticles[cellNum]] = i;
            }

            float minDist = 2.0f * scene.particleRadius;
            float minDistSq = minDist * minDist;

            // is this even important or can we skip pushParticlesIter
            runForNSteps(scene.pushParticlesIter) {
                for (int p = 0; p < particles->size; p++) {
                    float px = particles->px[p];
                    float py = particles->py[p];
                    
                    int onGridX = (int) floor(px * pInvSpacing);
                    int onGridY = (int) floor(py * pInvSpacing);

                    int x0 = fsim::Max(onGridX-1, 0);
                    int y0 = fsim::Max(onGridY-1, 0);
                    int x1 = fsim::Min(onGridX+1, pNumX-1);
                    int y1 = fsim::Min(onGridY+1, pNumY-1);

                    for (int x = x0; x<=x1; x++) {
                        for (int y = y0; y<=y1; y++) {
                            int cellNum = x * pNumY + y;
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

        virtual void handleParticleCollisions(float time) {
            // this doesn't use the time parameters
            MacGrid *mac = (MacGrid*) grid;
            const float h = mac->gridScale;

            // keep particle centres one solid cell + one radius inside the walls
            const float minX = h + scene.particleRadius;
            const float maxX = (mac->size - 1) * h - scene.particleRadius;
            const float minY = h + scene.particleRadius;
            const float maxY = (mac->size - 1) * h - scene.particleRadius;

            for (int i = 0; i < particles->size; i++) {
                if (particles->px[i] > maxX) {
                    particles->px[i] = maxX;
                    particles->vx[i] *= -scene.particleCollisionDamping;
                }
                else if (particles->px[i] < minX) {
                    particles->px[i] = minX;
                    particles->vx[i] *= -scene.particleCollisionDamping;
                }
                
                if (particles->py[i] > maxY) {
                    particles->py[i] = maxY;
                    particles->vy[i] *= -scene.particleCollisionDamping;
                }
                else if (particles->py[i] < minY) {
                    particles->py[i] = minY;
                    particles->vy[i] *= -scene.particleCollisionDamping;
                }
            }
        }

        virtual void transferVelocitiesToGrid() {
            MacGrid *mac = (MacGrid*) grid;
            const int numCells = grid->volume();

            for (int i = 0; i < numCells; i++) {
                mac->Vx[i] = 0.0f;
                mac->Vy[i] = 0.0f;
            }

            mac->resetFluidCells();

            for (int i = 0; i < particles->size; i++) {
                float x = particles->px[i];
                float y = particles->py[i];

                auto [cellX, cellY] = mac->worldCoordToGridCoord(x, y);

                mac->setFluid(cellX, cellY);
            }

            transferVelocityComponentToGrid(0);
            transferVelocityComponentToGrid(1);

            // FLIP delta must be measured against the freshly splatted field
            for (int i = 0; i < numCells; i++) {
                mac->prevVx[i] = mac->Vx[i];
                mac->prevVy[i] = mac->Vy[i];
            }
        }

        virtual void transferVelocitiesToParticles() {
            MacGrid *mac = (MacGrid*) grid;
            for (int component = 0; component < 2; component++) {
                AXIS axis = component == 0 ? X_AXIS : Y_AXIS;
                for (int i = 0; i < particles->size; i++) {
                    float x = particles->px[i];
                    float y = particles->py[i];

                    float picVelocity = mac->interpolateVelocity(x, y, axis);
                    float deltaVelocity = mac->interpolateVelocityDelta(x, y, axis, mac->prevVx, mac->prevVy);

                    float oldParticleVelocity = component == 0 ? particles->vx[i] : particles->vy[i];
                    float flipVelocity = oldParticleVelocity + deltaVelocity;

                    float newVelocity = (1.0f - scene.flipRatio) * picVelocity + scene.flipRatio * flipVelocity;

                    if (component == 0) {
                        particles->vx[i] = newVelocity;
                    }
                    else {
                        particles->vy[i] = newVelocity;
                    }
                }
            }
        }

        void transferVelocityComponentToGrid(int component) {
            MacGrid *mac = (MacGrid*) grid;
            float *velocity = component == 0 ? mac->Vx : mac->Vy;
            float *weights = makeFloatArray(((int)mac->volume()));

            for (int i = 0; i < particles->size; i++) {
                float x = particles->px[i];
                float y = particles->py[i];

                MacInterpolation s = mac->interpolation(x, y, component == 0 ? X_AXIS : Y_AXIS);

                int n0 = mac->indexOf(s.x0, s.y0);
                int n1 = mac->indexOf(s.x1, s.y0);
                int n2 = mac->indexOf(s.x1, s.y1);
                int n3 = mac->indexOf(s.x0, s.y1);

                float particleVelocity = component == 0 ? particles->vx[i] : particles->vy[i];

                velocity[n0] += particleVelocity * s.w0;
                weights[n0] += s.w0;
                velocity[n1] += particleVelocity * s.w1;
                weights[n1] += s.w1;
                velocity[n2] += particleVelocity * s.w2;
                weights[n2] += s.w2;
                velocity[n3] += particleVelocity * s.w3;
                weights[n3] += s.w3;
            }

            for (int i = 0; i < mac->volume(); i++) {
                if (weights[i] > 0.0f) {
                    velocity[i] /= weights[i];
                }
            }

            for (int x = 0; x < mac->size; x++) {
                for (int y = 0; y < mac->size; y++) {
                    int index = mac->indexOf(x, y);
                    if (component == 0) {
                        bool solid = mac->isSolid(x, y) || (
                                x > 0 &&
                                mac->isSolid(x - 1, y)
                            );

                        if (solid) velocity[index] = mac->prevVx[index];
                    }
                    else {
                        bool solid = mac->isSolid(x, y) ||
                            (
                                y > 0 &&
                                mac->isSolid(x, y - 1)
                            );

                        if (solid) velocity[index] = mac->prevVy[index];
                    }
                }
            }

            deleteFloatArray(weights);
        }

        virtual void updateParticleDensity() {
            MacGrid *mac = (MacGrid*) grid;
            const float h = mac->gridScale;
            const float invSpacing = 1.0f / h;

            for (int i = 0; i < grid->volume(); i++) {
                particleDensity[i] = 0.0f;
            }

            // deposit particle density onto the MAC grid using bilinear interpolation
            for (int i = 0; i < particles->size; i++) {
                float x = particles->px[i];
                float y = particles->py[i];

                x = fsim::Clamp(x, h, (grid->size - 1) * h);
                y = fsim::Clamp(y, h, (grid->size - 1) * h);

                int x0 = fsim::Floor((x - 0.5f * h) * invSpacing);
                float tx = ((x - 0.5f * h) - x0 * h) * invSpacing;
                int x1 = fsim::Min(x0 + 1, grid->size - 2);

                int y0 = fsim::Floor((y - 0.5f * h) * invSpacing);
                float ty = ((y - 0.5f * h) - y0 * h) * invSpacing;
                int y1 = fsim::Min(y0 + 1, grid->size - 2);

                float sx = 1.0f - tx;
                float sy = 1.0f - ty;

                if (x0 >= 0 && x0 < grid->size && y0 >= 0 && y0 < grid->size) 
                    particleDensity[mac->indexOf(x0, y0)] += sx * sy;
                if (x1 >= 0 && x1 < grid->size && y0 >= 0 && y0 < grid->size) 
                    particleDensity[mac->indexOf(x1, y0)] += tx * sy;
                if (x1 >= 0 && x1 < grid->size && y1 >= 0 && y1 < grid->size) 
                    particleDensity[mac->indexOf(x1, y1)] += tx * ty;
                if (x0 >= 0 && x0 < grid->size && y1 >= 0 && y1 < grid->size) 
                    particleDensity[mac->indexOf(x0, y1)] += sx * ty;
            }

            if (particleRestDensity == 0) {
                updateParticleRestDensity();
            }
        }

        void updateParticleRestDensity() {
            MacGrid *mac = (MacGrid*) grid;
            float sum = 0.0f;
            int numFluidCells = 0;

            for (int x = 0; x < mac->size; x++) {
                for (int y = 0; y < mac->size; y++) {
                    if (!mac->isFluid(x, y)) {
                        continue;
                    }
                    int index = mac->indexOf(x, y);

                    sum += particleDensity[index];
                    numFluidCells++;
                }
            }

            // particleRestDensity should be a constants
            // start of the simulation might point to some invalid state
            if (numFluidCells > 0) {
                particleRestDensity = sum / numFluidCells;
            }
        }


        virtual void solvePressure(float time) {
            MacGrid *mac = (MacGrid*) grid;
            int N = mac->size;
            float pressureCorrection = mac->fluidDensity * mac->gridScale / time;

            for (int i = 0; i < N*N; i++) mac->pressure[i] = 0.0f;

            runForNSteps(scene.iterations) {
                for (int i = 1; i < N-1; i++) {
                    for (int j = 1; j < N-1; j++) {
                        if (!mac->isFluid(i, j)) continue;

                        float sx0 = mac->isSolid(i-1, j) ? 0.0f : 1.0f;
                        float sx1 = mac->isSolid(i+1, j) ? 0.0f : 1.0f;
                        float sy0 = mac->isSolid(i, j-1) ? 0.0f : 1.0f;
                        float sy1 = mac->isSolid(i, j+1) ? 0.0f : 1.0f;
                        float s   = sx0 + sx1 + sy0 + sy1;
                        if (s == 0.0f) continue;

                        float div = mac->divergence(i, j);

                        if (particleRestDensity > 0.0) {
                            float compression =
                                particleDensity[mac->indexOf(i, j)] - particleRestDensity;

                            if (compression > 0.0f) {
                                div -= compression;
                            }
                        }

                        float pressureDelta = -div / s * scene.overRelaxation;
                        // likely unused and can be removed
                        // compare with other pressure solvers
                        mac->pressure[mac->indexOf(i, j)] += pressureCorrection * pressureDelta;

                        // likely pressureCorrection not needed either?
                        mac->Vx[mac->indexOf(i,   j)] -= sx0 * pressureDelta;
                        mac->Vx[mac->indexOf(i+1, j)] += sx1 * pressureDelta;
                        mac->Vy[mac->indexOf(i, j  )] -= sy0 * pressureDelta;
                        mac->Vy[mac->indexOf(i, j+1)] += sy1 * pressureDelta;
                    }
                }
            }
        }

    };
};
