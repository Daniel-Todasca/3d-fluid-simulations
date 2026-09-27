#pragma once

#include "../shared/ISimulation.hpp"
#include "ParticleCollection.hpp"
#include "../stable-fluids/FluidCube.hpp"
#include "MacGrid.hpp"
#include "PushParticlesSolver.hpp"
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
        }

        ~FlipSimulation2d() {
            if (grid) delete grid;
            if (particles) delete particles;

            if (pushSolver) delete pushSolver;
            if (particleDensity) deleteFloatArray(particleDensity);
            if (transferWeights) deleteFloatArray(transferWeights);
        }

        virtual void step() override {
            if (!particleDensity) particleDensity = makeFloatArray(((int)grid->volume()));
            if (!transferWeights) transferWeights = makeFloatArray(((int)grid->volume()));
            if (!pushSolver) pushSolver = new PushParticlesSolver(scene.particleRadius, scene.width, scene.height);

            applyGravity(grid->time_step * this->timeSpeedUp);
            transferVelocitiesToGrid();
            loopGrid(grid->time_step * this->timeSpeedUp);
            transferVelocitiesToParticles();
            loopParticles(grid->time_step * timeSpeedUp);
        }

        FluidCube *getGrid() const {
            return grid;
        }
        ParticleCollection *getParticles() const {
            return particles;
        }
        const Scene &getScene() const {
            return scene;
        }

    protected:
        Scene scene;
        FluidCube *grid;
        ParticleCollection *particles;
        PushParticlesSolver *pushSolver = nullptr;

        float *particleDensity = nullptr, *transferWeights = nullptr;

        float particleRestDensity = 0, timeSpeedUp = 10.0f, mass = 50.0f;

        virtual void applyGravity(float time) {
            for (int p = 0; p < particles->size; p++) {
                particles->vy[p] -= scene.gravity * this->mass * time;
            }
        }

        virtual void loopParticles(float time) {
            moveParticles(time);
            pushParticles();
            handleParticleCollisions();
        }

        virtual void loopGrid(float time) {
            updateParticleDensity();
            solvePressure(time);
        }

    private:
        float blendFlipAndPic(float flipRatio, float flipVelocity, float picVelocity) {
            return (1.0f - flipRatio) * picVelocity + flipRatio * flipVelocity;
        }

        float compensateDrift(float divergence, float cellDensity) {
            return divergence;
            if (particleRestDensity <= 0.0f) return divergence;

            float compression = cellDensity - particleRestDensity;
            if (compression > 0.0f) {
                divergence -= compression;
            }
            return divergence;
        }

        void moveParticles(float time) {
            for (int p=0; p < particles->size; p++) {
                particles->px[p] += particles->vx[p] * time;
                particles->py[p] += particles->vy[p] * time;
            }
        }

        void pushParticles() {
            // notice it is not dependent on time
            pushSolver->solve(particles, scene.pushParticlesIter);
        }

        void handleParticleCollisions() {
            // notice it is not dependent on time
            MacGrid *mac = (MacGrid*) grid;
            const float gridCellSpacing = mac->gridScale;

            const float minX = gridCellSpacing + scene.particleRadius;
            const float maxX = (mac->size - 1) * gridCellSpacing - scene.particleRadius;
            const float minY = gridCellSpacing + scene.particleRadius;
            const float maxY = (mac->size - 1) * gridCellSpacing - scene.particleRadius;

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

        void transferVelocitiesToGrid() {
            MacGrid *mac = (MacGrid*) grid;

            mac->resetVelocity();
            mac->resetFluidCells();

            for (int i = 0; i < particles->size; i++) {
                float x = particles->px[i];
                float y = particles->py[i];

                auto [cellX, cellY] = mac->worldCoordToGridCoord(x, y);
                mac->setFluid(cellX, cellY);
            }

            transferVelocityComponentToGrid(X_AXIS);
            transferVelocityComponentToGrid(Y_AXIS);

            mac->saveVelocitySnapshot();
        }

        virtual void transferVelocitiesToParticles() {
            MacGrid *mac = (MacGrid*) grid;
            for (int component = 0; component < 2; component++) {
                AXIS axis = component == 0 ? X_AXIS : Y_AXIS;
                for (int p = 0; p < particles->size; p++) {
                    float x = particles->px[p];
                    float y = particles->py[p];

                    float picVelocity, deltaVelocity;
                    std::tie(picVelocity, deltaVelocity) = mac->interpolateStep(x, y, axis);

                    float oldParticleVelocity = component == 0 ? particles->vx[p] : particles->vy[p];
                    float flipVelocity = oldParticleVelocity + deltaVelocity;

                    float newVelocity = blendFlipAndPic(scene.flipRatio, flipVelocity, picVelocity);

                    if (component == 0) {
                        particles->vx[p] = newVelocity;
                    }
                    else {
                        particles->vy[p] = newVelocity;
                    }
                }
            }
        }

        void transferVelocityComponentToGrid(AXIS axis) {
            MacGrid *mac = (MacGrid*) grid;
            float *velocity = axis == X_AXIS ? mac->Vx : mac->Vy;
            for (int i = 0; i < mac->volume(); i++) {
                transferWeights[i] = 0.0f;
            }

            for (int p = 0; p < particles->size; p++) {
                float x = particles->px[p];
                float y = particles->py[p];

                MacInterpolation vel4 = mac->interpolation(x, y, axis);

                int n0 = mac->indexOf(vel4.x0, vel4.y0);
                int n1 = mac->indexOf(vel4.x1, vel4.y0);
                int n2 = mac->indexOf(vel4.x1, vel4.y1);
                int n3 = mac->indexOf(vel4.x0, vel4.y1);

                float particleVelocity = axis == X_AXIS ? particles->vx[p] : particles->vy[p];

                velocity[n0] += particleVelocity * vel4.w0;
                transferWeights[n0] += vel4.w0;
                velocity[n1] += particleVelocity * vel4.w1;
                transferWeights[n1] += vel4.w1;
                velocity[n2] += particleVelocity * vel4.w2;
                transferWeights[n2] += vel4.w2;
                velocity[n3] += particleVelocity * vel4.w3;
                transferWeights[n3] += vel4.w3;
            }

            for (int i = 0; i < mac->volume(); i++) {
                if (transferWeights[i] > 0.0f) {
                    velocity[i] /= transferWeights[i];
                }
            }

            float *prevVelocity = axis == X_AXIS ? mac->prevVx : mac->prevVy;
            int dx = axis == X_AXIS ? 1 : 0;
            int dy = axis == X_AXIS ? 0 : 1;

            for (int x = 0; x < mac->size; x++) {
                for (int y = 0; y < mac->size; y++) {
                    int index = mac->indexOf(x, y);
                    bool solid = mac->isSolid(x, y) || (x - dx >= 0 && y - dy >= 0 && mac->isSolid(x - dx, y - dy));

                    if (solid) velocity[index] = prevVelocity[index];
                }
            }
        }

        void updateParticleDensity() {
            MacGrid *mac = (MacGrid*) grid;
            const float gridCellSpacing = mac->gridScale;
            const float pInvSpacing = 1.0f / gridCellSpacing;

            for (int i = 0; i < grid->volume(); i++) {
                particleDensity[i] = 0.0f;
            }

            for (int p = 0; p < particles->size; p++) {
                float x = particles->px[p];
                float y = particles->py[p];

                x = fsim::Clamp(x, gridCellSpacing, (grid->size - 1) * gridCellSpacing);
                y = fsim::Clamp(y, gridCellSpacing, (grid->size - 1) * gridCellSpacing);

                int x0 = fsim::Floor((x - 0.5f * gridCellSpacing) * pInvSpacing);
                float tx = ((x - 0.5f * gridCellSpacing) - x0 * gridCellSpacing) * pInvSpacing;
                int x1 = fsim::Min(x0 + 1, grid->size - 2);

                int y0 = fsim::Floor((y - 0.5f * gridCellSpacing) * pInvSpacing);
                float ty = ((y - 0.5f * gridCellSpacing) - y0 * gridCellSpacing) * pInvSpacing;
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
            LOG("Updating particle rest density");

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

            // could change particleRestDensity to a constant
            if (numFluidCells > 0) {
                particleRestDensity = sum / numFluidCells;
            }

            LOG(particleRestDensity);
        }

        void solvePressure(float time) {
            MacGrid *mac = (MacGrid*) grid;

            runForNSteps(scene.iterations) {
                for (int x = 1; x < mac->size-1; x++) {
                    for (int y = 1; y < mac->size-1; y++) {
                        if (!mac->isFluid(x, y)) continue;

                        float sx0 = mac->isSolid(x-1, y) ? 0.0f : 1.0f;
                        float sx1 = mac->isSolid(x+1, y) ? 0.0f : 1.0f;
                        float sy0 = mac->isSolid(x, y-1) ? 0.0f : 1.0f;
                        float sy1 = mac->isSolid(x, y+1) ? 0.0f : 1.0f;
                        float s = sx0 + sx1 + sy0 + sy1;
                        if (s == 0.0f) continue;

                        float divergence = mac->divergence(x, y);

                        int index = mac->indexOf(x, y);
                        divergence = compensateDrift(divergence, particleDensity[index]);

                        float pressureDelta = -divergence / s * scene.overRelaxation;

                        mac->Vx[index] -= sx0 * pressureDelta;
                        mac->Vx[mac->indexOf(x+1, y)] += sx1 * pressureDelta;
                        mac->Vy[index] -= sy0 * pressureDelta;
                        mac->Vy[mac->indexOf(x, y+1)] += sy1 * pressureDelta;
                    }
                }
            }
        }

    };
};
