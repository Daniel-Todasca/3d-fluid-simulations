#pragma once

#include "FluidCube.hpp"

namespace fsim {

    enum VoxelCellType {
        SOLID = 0,
        FLUID = 1,
        AIR = 2
    };

    class VoxelizedCube : public FluidCube {
    public:

        VoxelizedCube() : VoxelizedCube(
            CUBE_SIZE_DEFAULT,
            CUBE_DIFFUSION_DEFAULT,
            CUBE_VISCOSITY_DEFAULT,
            CUBE_TIMESTEP_DEFAULT
        ) { }

        VoxelizedCube(
            int size,
            float diffusion,
            float viscosity,
            float time_step
        ) : FluidCube(size, diffusion, viscosity, time_step) {
            cellType = new VoxelCellType[size*size];
            for (int index = 0; index < size * size; index++) {
                cellType[index] = AIR;
            }
            classifyCells();
        }

        ~VoxelizedCube() {
            delete[] cellType;
        }

        bool isSolid(int x, int y) {
            if (x < 0 || y < 0 || x >= size || y >= size) return false;
            return cellType[indexOf(x, y)] == SOLID;
        }
        bool isFluid(int x, int y) {
            if (x < 0 || y < 0 || x >= size || y >= size) return false;
            return cellType[indexOf(x, y)] == FLUID;
        }
        bool isAir(int x, int y) {
            if (x < 0 || y < 0 || x >= size || y >= size) return false;
            return cellType[indexOf(x, y)] == AIR;
        }

        bool isSurface(int x, int y) {
            // warn: there can be indexes out of bounds, caller must be careful
            return isFluid(x, y) && (
                isAir(x-1, y) || isAir(x+1, y) ||
                isAir(x, y-1) || isAir(x, y+1)
            );
        }

        void classifyCells() {
            for (int x = 1; x < size-1; x++) {
                for (int y = 1; y < size-1; y++) {
                    if (density[indexOf(x, y)] > MIN_DENSITY) {
                        cellType[indexOf(x, y)] = FLUID;
                    } else {
                        cellType[indexOf(x, y)] = AIR;
                    }
                }
            }

            for (int n = 0; n < size; n++) {
                cellType[indexOf(n, 0)] = SOLID;
                cellType[indexOf(n, size-1)] = SOLID;
                cellType[indexOf(0, n)] = SOLID;
                cellType[indexOf(size-1, n)] = SOLID;
            }
        }

        private:
        VoxelCellType *cellType;
    };

}
