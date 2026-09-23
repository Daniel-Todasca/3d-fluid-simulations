#pragma once

#include <cmath>
#include <cstring>

#include "../types/FluidCube.hpp"
#include "../types/VoxelizedCube.hpp"
#include "../types/Scene.hpp"

#include "../Constants.hpp"
#include "../Enums.hpp"

#include "../simulations/ISimulation.hpp"

namespace fsim {

    class SmokeSimulation2d : public ISimulation {

    public:
        SmokeSimulation2d(FluidCube *grid, const Scene &scene) {
            this->grid = grid;
            this->scene = scene;
        }

        ~SmokeSimulation2d() {
            if (grid) delete grid;
        }

        virtual void step() override {
            int N           = grid->size;
            float viscosity = grid->viscosity;
            float diffusion = grid->diffusion;
            float *density  = grid->density;
            float *Vx       = grid->Vx;
            float *Vy       = grid->Vy;
            float *Vx0      = grid->Vx0;
            float *Vy0      = grid->Vy0;
            float *dye      = grid->dye;

            // first part of the velocity equation
            diffuse(X_AXIS, Vx0, Vx, viscosity);
            diffuse(Y_AXIS, Vy0, Vy, viscosity);

            project(Vx0, Vy0, Vx, Vy);

            // second part of the velocity equation
            advect(X_AXIS, Vx, Vx0, Vx0, Vy0);
            advect(Y_AXIS, Vy, Vy0, Vx0, Vy0);

            project(Vx, Vy, Vx0, Vy0);

            // density equation
            diffuse(SCALAR, dye, density, diffusion);
            advect(SCALAR, density, dye, Vx, Vy);
        }

        FluidCube *getGrid() {
            return grid;
        }
        Scene getScene() {
            return scene;
        }

    protected:
        Scene scene;
        FluidCube *grid;

        virtual void diffuse (
            AXIS axis, 
            float *x, 
            float *field_prev, 
            float diffusion
        ) {
            float time_step = grid->time_step;
            int N = grid->size;
            int iters = scene.iterations;

            float coeff = time_step * diffusion * (N - 2) * (N - 2);
            gaussSeidel(axis, x, field_prev, coeff, 1 + 4 * coeff, iters);
        }

        virtual void project(
            float *Vx,
            float *Vy,
            float *pressure,
            float *divergence
        ) {
            int N = grid->size;
            int iters = scene.iterations;

            for (int x = 1; x < N-1; x++) {
                for (int y = 1; y < N-1; y++) {
                    divergence[grid->indexOf(x, y)] = -0.5f * (
                            Vx[grid->indexOf(x+1, y)] -
                            Vx[grid->indexOf(x-1, y)] +
                            Vy[grid->indexOf(x, y+1)] -
                            Vy[grid->indexOf(x, y-1)]
                        ) / N;
                    pressure[grid->indexOf(x, y)] = 0;
                }
            }

            setBounds(SCALAR, divergence, N);
            setBounds(SCALAR, pressure, N);
            gaussSeidel(SCALAR, pressure, divergence, 1, 4, iters);

            for (int x = 1; x < N-1; x++) {
                for (int y = 1; y < N-1; y++) {
                    Vx[grid->indexOf(x, y)] -= 0.5f * (pressure[grid->indexOf(x+1, y)] -
                                                         pressure[grid->indexOf(x-1, y)]) * N;
                    Vy[grid->indexOf(x, y)] -= 0.5f * (pressure[grid->indexOf(x, y+1)] -
                                                         pressure[grid->indexOf(x, y-1)]) * N;
                }
            }
            
            setBounds(X_AXIS, Vx, N);
            setBounds(Y_AXIS, Vy, N);
        }

        virtual void advect(
            AXIS axis,
            float *field,
            float *field_prev,
            float *Vx,
            float *Vy
        ) {
            int N       = grid->size;
            float dt    = grid->time_step;

            float Nf    = N-2;
            float dtx   = dt * Nf;
            float dty   = dt * Nf;

            for (int j = 1; j < N-1; j++) {
                for (int i = 1; i < N-1; i++) {
                    float x = ((float) i) - dtx * Vx[grid->indexOf(i, j)];
                    float y = ((float) j) - dty * Vy[grid->indexOf(i, j)];
                    
                    if (x < 0.5f)        x = 0.5f;
                    if (x > Nf + 0.5f)   x = Nf + 0.5f;
                    if (y < 0.5f)        y = 0.5f;
                    if (y > Nf + 0.5f)   y = Nf + 0.5f;

                    int i0 = floorf(x), i1 = i0 + 1;
                    int j0 = floorf(y), j1 = j0 + 1;

                    float s1 = x - i0, s0 = 1.0f - s1;
                    float t1 = y - j0, t0 = 1.0f - t1;

                    field[grid->indexOf(i, j)] = s0 * (
                        t0 * field_prev[grid->indexOf(i0, j0)] + t1 * field_prev[grid->indexOf(i0, j1)]
                    ) + s1 * (
                        t0 * field_prev[grid->indexOf(i1, j0)] + t1 * field_prev[grid->indexOf(i1, j1)]
                    );
                }
            }
            
            setBounds(axis, field, N);
        }

    protected:
        virtual void gaussSeidel (
            AXIS axis, 
            float *field, 
            float *field_prev, 
            float coeff, 
            float normalization, 
            int iters
        ) {
            fsim::GaussSeidel(field, field_prev, scene.overRelaxation, coeff, normalization, iters, grid->size);
            setBounds(axis, field, grid->size);
        }

        virtual void setBounds(
            AXIS axis, 
            float *field, 
            int N
        ) {
            for (int i = 1; i < N-1; i++) {
                field[grid->indexOf(i, 0)] = 
                    axis == 2 ? -field[grid->indexOf(i, 1)] 
                    : field[grid->indexOf(i, 1)];
                field[grid->indexOf(i, N-1)] = 
                    axis == 2 ? -field[grid->indexOf(i, N-2)] 
                    : field[grid->indexOf(i, N-2)];
            }

            for (int j = 1; j < N-1; j++) {
                field[grid->indexOf(0, j)] = 
                    axis == 1 ? -field[grid->indexOf(1, j)] 
                    : field[grid->indexOf(1, j)];
                field[grid->indexOf(N-1, j)] = 
                    axis == 1 ? -field[grid->indexOf(N-2, j)] 
                    : field[grid->indexOf(N-2, j)];
            }
            
            field[grid->indexOf(0, 0)]      = 0.5f * (field[grid->indexOf(1, 0)]     + field[grid->indexOf(0, 1)]);
            field[grid->indexOf(0, N-1)]    = 0.5f * (field[grid->indexOf(1, N-1)]   + field[grid->indexOf(0, N-2)]);
            field[grid->indexOf(N-1, 0)]    = 0.5f * (field[grid->indexOf(N-2, 0)]   + field[grid->indexOf(N-1, 1)]);
            field[grid->indexOf(N-1, N-1)]  = 0.5f * (field[grid->indexOf(N-2, N-1)] + field[grid->indexOf(N-1, N-2)]);
        }

    };

}
