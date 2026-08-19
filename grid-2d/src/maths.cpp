#pragma once

#include <cmath>

namespace fsim {
    float distance(float x, float y, float x2, float y2) {
        return sqrt((x-x2) * (x-x2) + (y-y2) * (y-y2));
    }
}
