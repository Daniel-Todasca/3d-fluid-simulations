#pragma once

#include <cmath>
#include <cstring>

#include "../types/FluidCube.hpp"
#include "../types/VoxelizedCube.hpp"
#include "../simulations/SmokeSimulation2d.hpp"

#include "../Constants.hpp"
#include "../Enums.hpp"

namespace fsim {

    class WaterSimulation2d : public SmokeSimulation2d {
    public:
        WaterSimulation2d(FluidCube *grid) : WaterSimulation2d(grid, WATER_SIMULATION_2D_GRAVITY_DEFAULT) { }
        WaterSimulation2d(FluidCube *grid, float gravity) : SmokeSimulation2d(grid), gravity(gravity) { }

        virtual void step() override {
            VoxelizedCube* cube = (VoxelizedCube*) getGrid();

            cube->classifyCells();

            applyGravity();

            SmokeSimulation2d::step();

            cube->classifyCells();
            smoothenEvaporatingWater();
        }

    private:
        float gravity;

        void applyGravity() {
            VoxelizedCube* cube = (VoxelizedCube*) getGrid();
            for (int x = 1; x < cube->size-1; x++) {
                for (int y = 1; y < cube->size-1; y++) {
                    if (cube->isFluid(x, y)) {
                        cube->addFrameVelocity(x, y, 0.0f, -gravity);
                    }
                }
            }
        }

        // counter semi-Lagrangian dissipation fully submerged fluid cells so the body of water does not disappear
        void smoothenEvaporatingWater() {
            VoxelizedCube* cube = (VoxelizedCube*) getGrid();
            int N = grid->size;
            for (int y = 1; y < N-1; y++) {
                for (int x = 1; x < N-1; x++) {
                    int idx = grid->indexOf(x, y);
                    if (grid->density[idx] < 0.0f)      grid->density[idx] = 0.0f;
                    else if (grid->density[idx] > 1.0f) grid->density[idx] = 1.0f;

                    // only top up cells that are still nearly full
                    if (grid->density[idx] > WATER_SHARPEN_THRESHOLD &&
                        cube->isFluid(x-1, y) && cube->isFluid(x+1, y) &&
                        cube->isFluid(x, y-1) && cube->isFluid(x, y+1)) {
                        grid->density[idx] = 1.0f;
                    }
                }
            }
        }

        virtual void diffuse(
            AXIS axis, 
            float *x, 
            float *field_prev, 
            float diffusion
        ) override {
            if (axis == SCALAR) {
                memcpy(x, field_prev, sizeof(float) * grid->volume());
            }
            else {
               SmokeSimulation2d::diffuse(axis, x, field_prev, diffusion);
            }
        }

        virtual void project(
            float *Vx,
            float *Vy,
            float *pressure,
            float *divergence
        ) override {
            int N = grid->size;
            VoxelizedCube* cube = (VoxelizedCube*) getGrid();

            for (int y = 1; y < N-1; y++) {
                for (int x = 1; x < N-1; x++) {
                    pressure[grid->indexOf(x, y)] = 0;
                    if (!cube->isFluid(x, y)) {
                        divergence[grid->indexOf(x, y)] = 0;
                        continue;
                    }
                    divergence[grid->indexOf(x, y)] = -0.5f * (
                            Vx[grid->indexOf(x+1, y)] -
                            Vx[grid->indexOf(x-1, y)] +
                            Vy[grid->indexOf(x, y+1)] -
                            Vy[grid->indexOf(x, y-1)]
                        ) / N;
                }
            }

            gaussSeidel(SCALAR, pressure, divergence, 1, 4, WATER_PRESSURE_ITER);

            for (int y = 1; y < N-1; y++) {
                for (int x = 1; x < N-1; x++) {
                    if (!cube->isFluid(x, y)) continue;
                    // neumann at solid faces (zero gradient into the wall), so we never
                    // read the uninitialised violent border pressure cells
                    float p  = pressure[grid->indexOf(x, y)];
                    float p1 = cube->isSolid(x-1, y) ? p : pressure[grid->indexOf(x-1, y)];
                    float p2 = cube->isSolid(x+1, y) ? p : pressure[grid->indexOf(x+1, y)];
                    float p3 = cube->isSolid(x, y-1) ? p : pressure[grid->indexOf(x, y-1)];
                    float p4 = cube->isSolid(x, y+1) ? p : pressure[grid->indexOf(x, y+1)];
                    Vx[grid->indexOf(x, y)] -= 0.5f * (p2 - p1) * N;
                    Vy[grid->indexOf(x, y)] -= 0.5f * (p4 - p3) * N;
                }
            }

            setBounds(X_AXIS, Vx, N);
            setBounds(Y_AXIS, Vy, N);
        }
        
        virtual void gaussSeidel (
            AXIS axis, 
            float *field, 
            float *field_prev, 
            float coeff, 
            float normalization, 
            int iters
        ) override {
            int N = grid->size;
            if (axis != SCALAR) {
                SmokeSimulation2d::gaussSeidel(axis, field, field_prev, coeff, normalization, iters);
                return;
            }

            // free-surface pressure solver
            // air pressure stays 0, so it acts as the Dirichlet term
            VoxelizedCube* cube = (VoxelizedCube*) grid;
            fsim::GaussSeidel(
                field, field_prev, coeff, iters, N,
                [cube](int x, int y) { return cube->isFluid(x, y); },
                [cube](int x, int y) { return cube->isSolid(x, y) ? 0.0f : 1.0f; }
            );
        }

        virtual void setBounds(
            AXIS axis, 
            float *field, 
            int N
        ) override {
            SmokeSimulation2d::setBounds(axis, field, N);
        }

    };

}
