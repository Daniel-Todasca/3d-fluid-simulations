#pragma once

#include <tuple>

#include "../stable-fluids/VoxelizedCube.hpp"
#include "../shared/Enums.hpp"

namespace fsim {

    struct MacInterpolation {
        int x0, x1, y0, y1;
        float w0, w1, w2, w3;
    };

    class MacGrid : public VoxelizedCube {
    public:
        float *prevVx;   // grid velocity snapshot before pressure equation
        float *prevVy;
        float fluidDensity;   // fluid density (kg/m^3)
        float gridScale; // cell spacing (domain width / grid size)

        MacGrid(int size, float diffusion, float viscosity, float time_step, float fluidDensity = 10.0f, float h = -1.0f)
            : VoxelizedCube(size, diffusion, viscosity, time_step),
              fluidDensity(fluidDensity), gridScale(h > 0.0f ? h : 1.0f / size) {
            this->prevVx = makeFloatArray(size * size);
            this->prevVy = makeFloatArray(size * size);

            this->viscosity = 0; // currently unused
            this->fluidDensity = 0;   // currently unused
        }

        ~MacGrid() {
            deleteFloatArray(this->prevVx);
            deleteFloatArray(this->prevVy);
        }

        std::tuple<int, int> worldCoordToGridCoord(float worldX, float worldY) const {
            int gridX = fsim::Clamp((int)fsim::Floor(worldX / gridScale), 0, size - 1);
            int gridY = fsim::Clamp((int)fsim::Floor(worldY / gridScale), 0, size - 1);
            return { gridX, gridY };
        }

        std::tuple<float, float> gridCoordToWorldCoord(int gridX, int gridY) const {
            return { (gridX + 0.5f) * gridScale, (gridY + 0.5f) * gridScale };
        }

        void saveVelocitySnapshot() {
            for (int i = 0; i < volume(); i++) {
                prevVx[i] = Vx[i];
                prevVy[i] = Vy[i];
            }
        }

        void resetVelocity() {
            for (int i = 0; i < volume(); i++) {
                Vx[i] = 0.0f;
                Vy[i] = 0.0f;
            }
        }

        float divergence(int x, int y) const {
            // check off by one boundary indexes
            int xRight = (x + 1 < size) ? x + 1 : x;
            int yTop = (y + 1 < size) ? y + 1 : y;
            return Vx[indexOf(xRight, y)] - Vx[indexOf(x, y)] + Vy[indexOf(x, yTop)] - Vy[indexOf(x, y)];
        }

        MacInterpolation interpolation(float x, float y, AXIS axis) const {
            return axis == X_AXIS ? interpolation(x, y, 0.0f, 0.5f * gridScale)
                                  : interpolation(x, y, 0.5f * gridScale, 0.0f);
        }

        // returns { picVelocity, deltaVelocity }
        std::tuple<float, float> interpolateStep(float x, float y, AXIS axis) const {
            MacInterpolation vel4 = interpolation(x, y, axis);
            const float *current = axis == X_AXIS ? Vx : Vy;
            const float *previous = axis == X_AXIS ? prevVx : prevVy;

            auto validCell = [this, axis](int cx, int cy) {
                return axis == X_AXIS ? validWaterCellHorizontal(cx, cy)
                                      : validWaterCellVertical(cx, cy);
            };

            int n0 = indexOf(vel4.x0, vel4.y0);
            int n1 = indexOf(vel4.x1, vel4.y0);
            int n2 = indexOf(vel4.x1, vel4.y1);
            int n3 = indexOf(vel4.x0, vel4.y1);

            float weight = 0.0f;
            float value = 0.0f;
            float deltaValue = 0.0f;

            if (validCell(vel4.x0, vel4.y0)) { value += current[n0] * vel4.w0; deltaValue += vel4.w0 * (current[n0] - previous[n0]); weight += vel4.w0; }
            if (validCell(vel4.x1, vel4.y0)) { value += current[n1] * vel4.w1; deltaValue += vel4.w1 * (current[n1] - previous[n1]); weight += vel4.w1; }
            if (validCell(vel4.x1, vel4.y1)) { value += current[n2] * vel4.w2; deltaValue += vel4.w2 * (current[n2] - previous[n2]); weight += vel4.w2; }
            if (validCell(vel4.x0, vel4.y1)) { value += current[n3] * vel4.w3; deltaValue += vel4.w3 * (current[n3] - previous[n3]); weight += vel4.w3; }

            if (weight == 0.0f) return { 0.0f, 0.0f };
            return { value / weight, deltaValue / weight };
        }

    private:
        bool validWaterCellHorizontal(int x, int y) const {
            return !isAir(x, y) || !isAir(x - 1, y);
        }

        bool validWaterCellVertical(int x, int y) const {
            return !isAir(x, y) || !isAir(x, y - 1);
        }
        
        MacInterpolation interpolation(float x, float y, float dx, float dy) const {
            x = fsim::Clamp(x, gridScale, (size - 1) * gridScale);
            y = fsim::Clamp(y, gridScale, (size - 1) * gridScale);

            int x0 = fsim::Min((int)fsim::Floor((x - dx) / gridScale), size - 2);
            int y0 = fsim::Min((int)fsim::Floor((y - dy) / gridScale), size - 2);

            int x1 = fsim::Min(x0 + 1, size - 2);
            int y1 = fsim::Min(y0 + 1, size - 2);

            float lerpX = ((x - dx) - x0 * gridScale) / gridScale;
            float lerpY = ((y - dy) - y0 * gridScale) / gridScale;

            float invLerpX = 1.0f - lerpX;
            float invLerpY = 1.0f - lerpY;

            return MacInterpolation{
                x0, x1, y0, y1,
                invLerpX * invLerpY, lerpX * invLerpY, lerpX * lerpY, invLerpX * lerpY
            };
        }

    };

}
