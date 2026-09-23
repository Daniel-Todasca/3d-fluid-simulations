#pragma once

#include "VoxelizedCube.hpp"

namespace fsim {

    struct MacInterpolation {
        int x0;
        int x1;
        int y0;
        int y1;

        float w0;
        float w1;
        float w2;
        float w3;
    };

    class MacGrid : public VoxelizedCube {
    public:
        float *pressure; // pressure
        float density;   // fluid density (kg/m^3)
        float gridScale; // cell spacing (domain width / grid size)

        float solidFlag(int x, int y) {
            return isSolid(x, y) ? 0.0f : 1.0f;
        }

        MacGrid(int size, float diffusion, float viscosity, float time_step,
                float density = 10.0f, float h = -1.0f)
            : VoxelizedCube(size, diffusion, viscosity, time_step),
              density(density), gridScale(h > 0.0f ? h : 1.0f / size) {
            pressure = makeFloatArray(size * size);
            this->viscosity = 0;
            this->density = 0;
        }

        ~MacGrid() {
            deleteFloatArray(pressure);
        }

        float divergence(int x, int y) {
            return Vx[indexOf(x + 1, y)]
                - Vx[indexOf(x, y)]
                + Vy[indexOf(x, y + 1)]
                - Vy[indexOf(x, y)];
        }

        bool validU(int x, int y) {
            return !isAir(x, y) || !isAir(x - 1, y);
        }

        bool validV(int x, int y) {
            return !isAir(x, y) || !isAir(x, y - 1);
        }

        MacInterpolation interpolation(float x, float y, int component) {
            const float halfH = 0.5f * gridScale;
            const float dx = component == 0 ? 0.0f : halfH;
            const float dy = component == 0 ? halfH : 0.0f;

            x = fsim::Clamp(x, gridScale, (size - 1) * gridScale);
            y = fsim::Clamp(y, gridScale, (size - 1) * gridScale);

            int x0 = fsim::Min((int)fsim::Floor((x - dx) / gridScale), size - 2);
            int y0 = fsim::Min((int)fsim::Floor((y - dy) / gridScale), size - 2);

            int x1 = fsim::Min(x0 + 1, size - 2);
            int y1 = fsim::Min(y0 + 1, size - 2);

            float tx = ((x - dx) - x0 * gridScale) / gridScale;
            float ty = ((y - dy) - y0 * gridScale) / gridScale;

            float sx = 1.0f - tx;
            float sy = 1.0f - ty;

            MacInterpolation result;

            result.x0 = x0;
            result.x1 = x1;
            result.y0 = y0;
            result.y1 = y1;

            result.w0 = sx * sy;
            result.w1 = tx * sy;
            result.w2 = tx * ty;
            result.w3 = sx * ty;

            return result;
        }

        float interpolateVelocity(float x, float y, int component) {
            MacInterpolation s = interpolation(x, y, component);
            const float *velocity = component == 0 ? Vx : Vy;

            int n0 = indexOf(s.x0, s.y0);
            int n1 = indexOf(s.x1, s.y0);
            int n2 = indexOf(s.x1, s.y1);
            int n3 = indexOf(s.x0, s.y1);

            float validWeight = 0.0f;
            float value = 0.0f;

            if (component == 0) {
                if (validU(s.x0, s.y0)) {
                    value += velocity[n0] * s.w0;
                    validWeight += s.w0;
                }
                if (validU(s.x1, s.y0)) {
                    value += velocity[n1] * s.w1;
                    validWeight += s.w1;
                }
                if (validU(s.x1, s.y1)) {
                    value += velocity[n2] * s.w2;
                    validWeight += s.w2;
                }
                if (validU(s.x0, s.y1)) {
                    value += velocity[n3] * s.w3;
                    validWeight += s.w3;
                }
            }
            else {
                if (validV(s.x0, s.y0)) {
                    value += velocity[n0] * s.w0;
                    validWeight += s.w0;
                }
                if (validV(s.x1, s.y0)) {
                    value += velocity[n1] * s.w1;
                    validWeight += s.w1;
                }
                if (validV(s.x1, s.y1)) {
                    value += velocity[n2] * s.w2;
                    validWeight += s.w2;
                }
                if (validV(s.x0, s.y1)) {
                    value += velocity[n3] * s.w3;
                    validWeight += s.w3;
                }
            }

            if (validWeight == 0.0f)
                return 0.0f;

            return value / validWeight;
        }

        float interpolateVelocityDelta(float x, float y, int component, const float *previousVx, const float *previousVy) {
            MacInterpolation s = interpolation(x, y, component);
            const float *current = component == 0 ? Vx : Vy;
            const float *previous = component == 0 ? previousVx : previousVy;

            int n0 = indexOf(s.x0, s.y0);
            int n1 = indexOf(s.x1, s.y0);
            int n2 = indexOf(s.x1, s.y1);
            int n3 = indexOf(s.x0, s.y1);

            float validWeight = 0.0f;
            float value = 0.0f;

            if (component == 0) {
                if (validU(s.x0, s.y0)) {
                    value += s.w0 * (current[n0] - previous[n0]);
                    validWeight += s.w0;
                }
                if (validU(s.x1, s.y0)) {
                    value += s.w1 * (current[n1] - previous[n1]);
                    validWeight += s.w1;
                }
                if (validU(s.x1, s.y1)) {
                    value += s.w2 * (current[n2] - previous[n2]);
                    validWeight += s.w2;
                }
                if (validU(s.x0, s.y1)) {
                    value += s.w3 * (current[n3] - previous[n3]);
                    validWeight += s.w3;
                }
            }
            else {
                if (validV(s.x0, s.y0)) {
                    value += s.w0 * (current[n0] - previous[n0]);
                    validWeight += s.w0;
                }
                if (validV(s.x1, s.y0)) {
                    value += s.w1 * (current[n1] - previous[n1]);
                    validWeight += s.w1;
                }
                if (validV(s.x1, s.y1)) {
                    value += s.w2 * (current[n2] - previous[n2]);
                    validWeight += s.w2;
                }
                if (validV(s.x0, s.y1)) {
                    value += s.w3 * (current[n3] - previous[n3]);
                    validWeight += s.w3;
                }
            }

            if (validWeight == 0.0f)
                return 0.0f;

            return value / validWeight;
        }

    };

}
