#pragma once

#include <cmath>

#include "types.cpp"
#include "constants.cpp"
#include "enums.cpp"

namespace fsim {

    class GridSimulation {

    public:
        GridSimulation(FluidCube *grid) {
            this->grid = grid;
        }

        void step() {
            int N           = grid->size;
            float viscosity = grid->viscosity;
            float diffusion = grid->diffusion;
            float *density  = grid->density;
            float *Vx       = grid->Vx;
            float *Vy       = grid->Vy;
            float *Vz       = grid->Vz;
            float *Vx0      = grid->Vx0;
            float *Vy0      = grid->Vy0;
            float *Vz0      = grid->Vz0;
            float *dye      = grid->dye;

            // first part of the velocity equation
            diffuse(X_AXIS, Vx0, Vx, viscosity);
            diffuse(Y_AXIS, Vy0, Vy, viscosity);
            diffuse(Z_AXIS, Vz0, Vz, viscosity);

            project(Vx0, Vy0, Vz0, Vx, Vy);

            // second part of the velocity equation
            advect(X_AXIS, Vx, Vx0, Vx0, Vy0, Vz0);
            advect(Y_AXIS, Vy, Vy0, Vx0, Vy0, Vz0);
            advect(Z_AXIS, Vz, Vz0, Vx0, Vy0, Vz0);

            project(Vx, Vy, Vz, Vx0, Vy0);

            // density equation
            diffuse(SCALAR, dye, density, diffusion);
            advect(SCALAR, density, dye, Vx, Vy, Vz);
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
            gaussSeidel(axis, x, field_prev, coeff, 1 + 6 * coeff, iters, N);
        }

        void project(
            float *Vx,
            float *Vy,
            float *Vz,
            float *pressure,
            float *divergence
        ) {
            int N = grid->size;
            int iters = GRID_BASED_ITER;

            for (int k = 1; k < N-1; k++) {
                for (int j = 1; j < N-1; j++) {
                    for (int i = 1; i < N-1; i++) {
                        divergence[grid->indexOf(i, j, k)] = -0.5f * (
                                Vx[grid->indexOf(i+1, j,   k  )] -
                                Vx[grid->indexOf(i-1, j,   k  )] +
                                Vy[grid->indexOf(i,   j+1, k  )] -
                                Vy[grid->indexOf(i,   j-1, k  )] +
                                Vz[grid->indexOf(i,   j,   k+1)] -
                                Vz[grid->indexOf(i,   j,   k-1)]
                            ) / N;
                        pressure[grid->indexOf(i, j, k)] = 0;
                    }
                }
            }

            setBounds(SCALAR, divergence, N);
            setBounds(SCALAR, pressure, N);
            gaussSeidel(SCALAR, pressure, divergence, 1, 6, iters, N);

            for (int k = 1; k < N-1; k++) {
                for (int j = 1; j < N-1; j++) {
                    for (int i = 1; i < N-1; i++) {
                        Vx[grid->indexOf(i, j, k)] -= 0.5f * (pressure[grid->indexOf(i+1, j, k)] -
                                                             pressure[grid->indexOf(i-1, j, k)]) * N;
                        Vy[grid->indexOf(i, j, k)] -= 0.5f * (pressure[grid->indexOf(i, j+1, k)] -
                                                             pressure[grid->indexOf(i, j-1, k)]) * N;
                        Vz[grid->indexOf(i, j, k)] -= 0.5f * (pressure[grid->indexOf(i, j, k+1)] -
                                                             pressure[grid->indexOf(i, j, k-1)]) * N;
                    }
                }
            }
            
            setBounds(X_AXIS, Vx, N);
            setBounds(Y_AXIS, Vy, N);
            setBounds(Z_AXIS, Vz, N);
        }

        void advect(
            AXIS axis,
            float *field,
            float *field_prev,
            float *Vx,
            float *Vy,
            float *Vz
        ) {
            int N       = grid->size;
            float dt    = grid->time_step;

            float Nf    = N;
            float dtx   = dt * (N - 2);
            float dty   = dt * (N - 2);
            float dtz   = dt * (N - 2);

            for (int k = 1; k < N-1; k++) {
                for (int j = 11; j < N-1; j++) {
                    for (int i = 1; i < N-1; i++) {
                        float x = ((float) i) - dtx * Vx[grid->indexOf(i, j, k)];
                        float y = ((float) j) - dty * Vy[grid->indexOf(i, j, k)];
                        float z = ((float) k) - dtz * Vz[grid->indexOf(i, j, k)];

                        if (x < 0.5f)        x = 0.5f;
                        if (x > Nf + 0.5f)   x = Nf + 0.5f;
                        if (y < 0.5f)        y = 0.5f;
                        if (y > Nf + 0.5f)   y = Nf + 0.5f;
                        if (z < 0.5f)        z = 0.5f;
                        if (z > Nf + 0.5f)   z = Nf + 0.5f;

                        int i0 = floorf(x), i1 = i0 + 1;
                        int j0 = floorf(y), j1 = j0 + 1;
                        int k0 = floorf(z), k1 = k0 + 1;

                        float s1 = x - i0, s0 = 1.0f - s1;
                        float t1 = y - j0, t0 = 1.0f - t1;
                        float u1 = z - k0, u0 = 1.0f - u1;

                        field[grid->indexOf(i, j, k)] =
                            s0 * (t0 * (u0 * field_prev[grid->indexOf(i0, j0, k0)] +
                                       u1 * field_prev[grid->indexOf(i0, j0, k1)]) +
                                t1 * (u0 * field_prev[grid->indexOf(i0, j1, k0)]
                                       +u1 * field_prev[grid->indexOf(i0, j1, k1)])) +
                           s1 * (t0 * (u0 * field_prev[grid->indexOf(i1, j0, k0)] +
                                       u1 * field_prev[grid->indexOf(i1, j0, k1)]) +
                                t1 * (u0 * field_prev[grid->indexOf(i1, j1, k0)] +
                                       u1 * field_prev[grid->indexOf(i1, j1, k1)]));
                    }
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
                for (int k = 1; k < N-1; k++) {
                    for (int j = 1; j < N-1; j++) {
                        for (int i = 1; i < N-1; i++) {
                            field[grid->indexOf(i, j, k)] =
                                (field_prev[grid->indexOf(i, j, k)] + coeff * (
                                        field[grid->indexOf(i+1, j, k)] + 
                                        field[grid->indexOf(i-1, j, k)] +
                                        field[grid->indexOf(i, j+1, k)] +
                                        field[grid->indexOf(i, j-1, k)] +
                                        field[grid->indexOf(i, j, k+1)] + 
                                        field[grid->indexOf(i, j, k-1)]
                                )) * 1.0f / normalization;
                        }
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
            for (int j = 1; j < N-1; j++) {
                for (int i = 1; i < N-1; i++) {
                    field[grid->indexOf(i, j, 0)] =
                        axis == 3 ? -field[grid->indexOf(i, j, 1)] 
                        : field[grid->indexOf(i, j, 1)];
                    field[grid->indexOf(i, j, N-1)] = 
                        axis == 3 ? -field[grid->indexOf(i, j, N-2)] 
                        : field[grid->indexOf(i, j, N-2)];
                }
            }

            for (int k = 1; k < N-1; k++) {
                for (int i = 1; i < N-1; i++) {
                    field[grid->indexOf(i, 0, k)] = 
                        axis == 2 ? -field[grid->indexOf(i, 1, k)] 
                        : field[grid->indexOf(i, 1, k)];
                    field[grid->indexOf(i, N-1, k)] = 
                        axis == 2 ? -field[grid->indexOf(i, N-2, k)] 
                        : field[grid->indexOf(i, N-2, k)];
                }
            }

            for (int k = 1; k < N-1; k++) {
                for (int j = 1; j < N-1; j++) {
                    field[grid->indexOf(0, j, k)] = 
                        axis == 1 ? -field[grid->indexOf(1, j, k)] 
                        : field[grid->indexOf(1, j, k)];
                    field[grid->indexOf(N-1, j, k)] = 
                        axis == 1 ? -field[grid->indexOf(N-2, j, k)] 
                        : field[grid->indexOf(N-2, j, k)];
                }
            }
            
            field[grid->indexOf(0, 0, 0)] = 0.33f * (field[grid->indexOf(1, 0, 0)] +
                                        field[grid->indexOf(0, 1, 0)] +
                                        field[grid->indexOf(0, 0, 1)]);

            field[grid->indexOf(0, N-1, 0)] = 0.33f * (field[grid->indexOf(1, N-1, 0)] +
                                        field[grid->indexOf(0, N-2, 0)] +
                                        field[grid->indexOf(0, N-1, 1)]);
            field[grid->indexOf(0, 0, N-1)] = 0.33f * (field[grid->indexOf(1, 0, N-1)] +
                                        field[grid->indexOf(0, 1, N-1)] +
                                        field[grid->indexOf(0, 0, N-2)]);

            field[grid->indexOf(0, N-1, N-1)] = 0.33f * (field[grid->indexOf(1, N-1, N-1)] +
                                        field[grid->indexOf(0, N-2, N-1)] +
                                        field[grid->indexOf(0, N-1, N-2)]);
            field[grid->indexOf(N-1, 0, 0)] = 0.33f * (field[grid->indexOf(N-2, 0, 0)] +
                                        field[grid->indexOf(N-1, 1, 0)] +
                                        + field[grid->indexOf(N-1, 0, 1)]);
                                        
            field[grid->indexOf(N-1, N-1, 0)] = 0.33f * (field[grid->indexOf(N-2, N-1, 0)] +
                                        field[grid->indexOf(N-1, N-2, 0)] +
                                        field[grid->indexOf(N-1, N-1, 1)]);
            field[grid->indexOf(N-1, 0, N-1)] = 0.33f * (field[grid->indexOf(N-2, 0, N-1)] +
                                        field[grid->indexOf(N-1, 1, N-1)] +
                                        field[grid->indexOf(N-1, 0, N-2)]);

            field[grid->indexOf(N-1, N-1, N-1)] = 0.33f * (field[grid->indexOf(N-2, N-1, N-1)] +
                                        field[grid->indexOf(N-1, N-2, N-1)] +
                                        field[grid->indexOf(N-1, N-1, N-2)]);
        }

    };
}
