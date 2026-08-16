#pragma once

#include "macros.cpp"
#include "constants.cpp"

namespace fsim {
    class FluidCube {

    public:
        void addDensity(int x, int y, int z, float amount) {
            density[indexOf(x, y, z)] += amount;
        }

        void addVelocity(int x, int y, int z, float vx, float vy, float vz) {
            Vx[indexOf(x, y, z)] += vx;
            Vy[indexOf(x, y, z)] += vy;
            Vz[indexOf(x, y, z)] += vz;
        }

        inline int indexOf(int x, int y, int z) {
            return x + y*size + z*size*size;
        }

    public:        
        int size;                   // all cubes will have the same size
        float diffusion;            // all cubes will have the same diffusion
        float viscosity;            // all cubes will have the same viscosity
        
        float time_step;            // length of timestep

        float *dye;          
        float *density;             // density array

        float *Vx, *Vy, *Vz;        // velocity array
        float *Vx0, *Vy0, *Vz0;     // previous velocity array

        FluidCube() : FluidCube(
            CUBE_SIZE_DEFAULT,
            CUBE_DIFFUSION_DEFAULT,
            CUBE_VISCOSITY_DEFAULT,
            CUBE_TIMESTEP_DEFAULT
        ) {}

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
            int N3  = N*N*N;

            dye     = makeFloatArray(N3);
            density = makeFloatArray(N3);
            Vx      = makeFloatArray(N3);
            Vy      = makeFloatArray(N3);
            Vz      = makeFloatArray(N3);
            Vx0     = makeFloatArray(N3);
            Vy0     = makeFloatArray(N3);
            Vz0     = makeFloatArray(N3);
        }

        ~FluidCube() {
            deleteFloatArray(dye);
            deleteFloatArray(density);
            deleteFloatArray(Vx);
            deleteFloatArray(Vy);
            deleteFloatArray(Vz);
            deleteFloatArray(Vx0);
            deleteFloatArray(Vy0);
            deleteFloatArray(Vz0);
        }

    };

}
