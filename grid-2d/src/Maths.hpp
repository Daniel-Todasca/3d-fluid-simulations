#pragma once

#include <cmath>

#include "Macros.hpp"

namespace fsim {
    template <typename T>
    T Clamp(T value, T lo, T hi) {
        return value < lo ? lo : (value > hi ? hi : value);
    }

    template <typename T>
    T Max(T a, T b) {
        return a > b ? a : b;
    }

    template <typename T>
    T Min(T a, T b) {
        return a < b ? a : b;
    }

    template <typename T>
    T Sqrt(T value) {
        return std::sqrt(value);
    }

    template <typename T>
    T Floor(T value) {
        return std::floor(value);
    }

    float Distance(float x, float y, float x2, float y2) {
        return sqrt((x-x2) * (x-x2) + (y-y2) * (y-y2));
    }

    void GaussSeidel (
        float *field, 
        float *field_prev, 
        float overRelaxation,
        float coeff, 
        float normalization, 
        int iters,
        int N
    ) {
        auto ToIndex = [N](int x, int y) {
            return x + y * N;
        };

        runForNSteps(iters) {
            for (int i = 1; i < N-1; i++) {
                for (int j = 1; j < N-1; j++) {
                    field[ToIndex(i, j)] = (
                        field_prev[ToIndex(i, j)] + coeff * (
                            field[ToIndex(i-1, j)] +
                            field[ToIndex(i+1, j)] +
                            field[ToIndex(i, j-1)] +
                            field[ToIndex(i, j+1)] 
                        ) 
                    ) * overRelaxation / normalization;
                }
            }
        }
    }

    // isActive(x, y) picks which cells to solve
    // weight(x, y) is a neighbour's coefficient (return 0 to drop it)
    template <typename IsActive, typename Weight>
    void GaussSeidel (
        float *field, 
        float *field_prev, 
        float overRelaxation,
        float coeff, 
        int iters,
        int N,
        IsActive isActive,
        Weight weight
    ) {
        auto ToIndex = [N](int x, int y) {
            return x + y * N;
        };

        runForNSteps(iters) {
            for (int y = 1; y < N-1; y++) {
                for (int x = 1; x < N-1; x++) {
                    if (!isActive(x, y)) continue;

                    float a1 = weight(x-1, y);
                    float a2 = weight(x+1, y);
                    float a3 = weight(x, y-1);
                    float a4 = weight(x, y+1);

                    float total = a1 + a2 + a3 + a4;
                    if (total == 0.0f) continue;

                    float sum = a1 * field[ToIndex(x-1, y)] + a2 * field[ToIndex(x+1, y)]
                              + a3 * field[ToIndex(x, y-1)] + a4 * field[ToIndex(x, y+1)];

                    field[ToIndex(x, y)] = overRelaxation * (field_prev[ToIndex(x, y)] + coeff * sum) / total;
                }
            }
        }
    }
}
