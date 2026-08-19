#pragma once

#include <cmath>

#include "types.cpp"
#include "constants.cpp"
#include "enums.cpp"

namespace fsim {

    class GridSimulation2d {

    public:
        GridSimulation2d(FluidCube *grid) {
            this->grid = grid;
        }

        void step() {
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

    private:
        FluidCube *grid;

        void diffuse (
            AXIS axis, 
            float *x, 
            float *field_prev, 
            float diffusion
        ) {
            float time_step = grid->time_step;
            int N = grid->size;
            int iters = GRID_BASED_ITER;

            float coeff = time_step * diffusion * (N - 2) * (N - 2);
            gaussSeidel(axis, x, field_prev, coeff, 1 + 4 * coeff, iters, N);
        }

        void project(
            float *Vx,
            float *Vy,
            float *pressure,
            float *divergence
        ) {
            int N = grid->size;
            int iters = GRID_BASED_ITER;

            for (int j = 1; j < N-1; j++) {
                for (int i = 1; i < N-1; i++) {
                    divergence[grid->indexOf(i, j)] = -0.5f * (
                            Vx[grid->indexOf(i+1, j)] -
                            Vx[grid->indexOf(i-1, j)] +
                            Vy[grid->indexOf(i, j+1)] -
                            Vy[grid->indexOf(i, j-1)]
                        ) / N;
                    pressure[grid->indexOf(i, j)] = 0;
                }
            }

            setBounds(SCALAR, divergence, N);
            setBounds(SCALAR, pressure, N);
            gaussSeidel(SCALAR, pressure, divergence, 1, 4, iters, N);

            for (int j = 1; j < N-1; j++) {
                for (int i = 1; i < N-1; i++) {
                    Vx[grid->indexOf(i, j)] -= 0.5f * (pressure[grid->indexOf(i+1, j)] -
                                                         pressure[grid->indexOf(i-1, j)]) * N;
                    Vy[grid->indexOf(i, j)] -= 0.5f * (pressure[grid->indexOf(i, j+1)] -
                                                         pressure[grid->indexOf(i, j-1)]) * N;
                }
            }
            
            setBounds(X_AXIS, Vx, N);
            setBounds(Y_AXIS, Vy, N);
        }

        void advect(
            AXIS axis,
            float *field,
            float *field_prev,
            float *Vx,
            float *Vy
        ) {
            int N       = grid->size;
            float dt    = grid->time_step;

            float Nf    = N-2;
            float dtx   = dt * (N - 2);
            float dty   = dt * (N - 2);

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

    private:
        void gaussSeidel (
            AXIS axis, 
            float *field, 
            float *field_prev, 
            float coeff, 
            float normalization, 
            int iters,
            int N
        ) {
            runForNSteps(iters) {
                for (int i = 1; i < N-1; i++) {
                    for (int j = 1; j < N-1; j++) {
                        field[grid->indexOf(i, j)] = (
                            field_prev[grid->indexOf(i, j)] + coeff * (
                                field[grid->indexOf(i-1, j)] +
                                field[grid->indexOf(i+1, j)] +
                                field[grid->indexOf(i, j-1)] +
                                field[grid->indexOf(i, j+1)] 
                            ) 
                        ) * 1.0f / normalization;
                    }
                }

                setBounds(axis, field, N);
            }
        }

        void setBounds(
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
