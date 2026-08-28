#pragma once

#include "../Macros.hpp"
#include "../Constants.hpp"
#include "../Maths.hpp"

namespace fsim {
    
    class FluidCube {
    public:
        inline float volume() {
            return size * size;
        }

        void addDensity(int x, int y, float amount) {
            density[indexOf(x, y)] += amount;
        }

        void addDensityToCircle(int x, int y, float radius, float amount) {
            int N = size;
            for (int x2 = 1; x2 < N-1; x2++) {
                for (int y2 = 1; y2 < N-1; y2++) {
                    float Distance = fsim::Distance((float) x, (float) y, (float) x2, (float) y2);
                    if (Distance > radius) continue;
                    density[indexOf(x2, y2)] += amount / (1.0f + Distance);
                }
            }
        }

        void addVelocity(int x, int y, float vx, float vy) {
            Vx[indexOf(x, y)] += vx;
            Vy[indexOf(x, y)] += vy;
        }

        void addFrameVelocity(int x, int y, float vx, float vy) {
            Vx[indexOf(x, y)] += vx * time_step;
            Vy[indexOf(x, y)] += vy * time_step;
        }

        inline int indexOf(int x, int y) {
            // if (x < 0 || y < 0 || x >= size || y >= size) return 0;
            return x + y*size;
        }

    public:        
        int size;                   // all cubes will have the same size
        float diffusion;            // all cubes will have the same diffusion
        float viscosity;            // all cubes will have the same viscosity
        
        float time_step;            // length of timestep

        float *dye;          
        float *density;             // density array

        float *Vx, *Vy;             // velocity
        float *Vx0, *Vy0;           // previous velocity

        FluidCube() : FluidCube(
            CUBE_SIZE_DEFAULT,
            CUBE_DIFFUSION_DEFAULT,
            CUBE_VISCOSITY_DEFAULT,
            CUBE_TIMESTEP_DEFAULT
        ) { }

        FluidCube(
            int size,
            float diffusion,
            float viscosity,
            float time_step
        ) {
            this->size      = size;
            this->diffusion = diffusion;
            this->viscosity = viscosity;
            this->time_step = time_step;

            int N   = size;
            int N2  = N*N;

            dye     = makeFloatArray(N2);
            density = makeFloatArray(N2);
            Vx      = makeFloatArray(N2);
            Vy      = makeFloatArray(N2);
            Vx0     = makeFloatArray(N2);
            Vy0     = makeFloatArray(N2);
        }

        ~FluidCube() {
            deleteFloatArray(dye);
            deleteFloatArray(density);
            deleteFloatArray(Vx);
            deleteFloatArray(Vy);
            deleteFloatArray(Vx0);
            deleteFloatArray(Vy0);
        }

    };

}
