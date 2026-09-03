#pragma once

#include "VoxelizedCube.hpp"

namespace fsim {

    class MacGrid : public VoxelizedCube {
    public:
        float *p;       // pressure
        float density;  // fluid density (kg/m^3)
        float h;        // cell spacing (domain width / grid size)

        float solidFlag(int x, int y) {
            return isSolid(x, y) ? 0.0f : 1.0f;
        }

        MacGrid(int size, float diffusion, float viscosity, float time_step,
                float density = 1000.0f, float h = -1.0f)
            : VoxelizedCube(size, diffusion, viscosity, time_step),
              density(density), h(h > 0.0f ? h : 1.0f / size) {
            p = makeFloatArray(size * size);
        }

        ~MacGrid() {
            deleteFloatArray(p);
        }
    };

}
