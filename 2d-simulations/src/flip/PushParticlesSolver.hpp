#pragma once

#include "ParticleCollection.hpp"
#include "../defines.hpp"
#include "../shared/Maths.hpp"

namespace fsim {

    // counts based on suffix indexes
    class PushParticlesSolver {
    public:
        PushParticlesSolver(float particleRadius, float width, float height) {
            pInvSpacing = 1.0f / (2.2f * particleRadius);
            pNumX = (int) floorf(width * pInvSpacing) + 1;
            pNumY = (int) floorf(height * pInvSpacing) + 1;
            numCells = pNumX * pNumY;

            minDist = 2.0f * particleRadius;
            minDistSq = minDist * minDist;

            numCellParticles = makeIntArray(numCells);
            firstCellParticles = makeIntArray(numCells + 1);
        }

        ~PushParticlesSolver() {
            if (numCellParticles) deleteIntArray(numCellParticles);
            if (firstCellParticles) deleteIntArray(firstCellParticles);
            if (cellParticleIds) deleteIntArray(cellParticleIds);
        }

        void solve(ParticleCollection *particles, int iterations) {
            if (!cellParticleIds) cellParticleIds = makeIntArray(particles->size);

            buildSpatialHash(particles);

            runForNSteps(iterations) {
                pushParticles(particles);
            }
        }

    private:
        void buildSpatialHash(ParticleCollection *particles) {
            for (int i = 0; i < numCells; i++) {
                numCellParticles[i] = 0;
            }

            for (int i = 0; i < particles->size; i++) {
                numCellParticles[cellIndex(particles->px[i], particles->py[i])] += 1;
            }

            int partialSum = 0;
            for (int i = 0; i < numCells; i++) {
                partialSum += numCellParticles[i];
                firstCellParticles[i] = partialSum;
            }
            firstCellParticles[numCells] = partialSum;

            for (int i = 0; i < particles->size; i++) {
                int cellNum = cellIndex(particles->px[i], particles->py[i]);
                firstCellParticles[cellNum]--;
                cellParticleIds[firstCellParticles[cellNum]] = i;
            }
        }

        void pushParticles(ParticleCollection *particles) {
            for (int p = 0; p < particles->size; p++) {
                float px = particles->px[p];
                float py = particles->py[p];

                int onGridX = (int) floorf(px * pInvSpacing);
                int onGridY = (int) floorf(py * pInvSpacing);

                int x0 = fsim::Max(onGridX - 1, 0);
                int y0 = fsim::Max(onGridY - 1, 0);
                int x1 = fsim::Min(onGridX + 1, pNumX - 1);
                int y1 = fsim::Min(onGridY + 1, pNumY - 1);

                for (int x = x0; x <= x1; x++) {
                    for (int y = y0; y <= y1; y++) {
                        int cellNum = x * pNumY + y;
                        int first = firstCellParticles[cellNum];
                        int last = firstCellParticles[cellNum + 1];

                        for (int i = first; i < last; i++) {
                            int particleId = cellParticleIds[i];
                            if (particleId == p) continue;

                            float dx = particles->px[particleId] - px;
                            float dy = particles->py[particleId] - py;
                            float distSq = dx * dx + dy * dy;
                            if (distSq == 0.0f || distSq > minDistSq) continue;

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

        int cellIndex(float x, float y) const {
            int cellX = fsim::Clamp((int) floorf(x * pInvSpacing), 0, pNumX - 1);
            int cellY = fsim::Clamp((int) floorf(y * pInvSpacing), 0, pNumY - 1);
            return cellX * pNumY + cellY;
        }

        float pInvSpacing;
        int pNumX, pNumY, numCells;
        int *numCellParticles = nullptr;
        int *firstCellParticles = nullptr;
        int *cellParticleIds = nullptr;
        float minDist, minDistSq;
    };

}
