#pragma once

#include "../Macros.hpp"

namespace fsim {
    class ParticleCollection {
    public:
        ParticleCollection(int numParticles) {
            this->size = numParticles;
            this->px = makeFloatArray(numParticles);
            this->py = makeFloatArray(numParticles);
            this->vx = makeFloatArray(numParticles);
            this->vy = makeFloatArray(numParticles);
        }

        ~ParticleCollection() {
            if (px) deleteFloatArray(px);
            if (py) deleteFloatArray(py);
            if (vx) deleteFloatArray(vx);
            if (vy) deleteFloatArray(vy);
        }

        int size;
        float *px, *py; // particle positions
        float *vx, *vy; // particle velocities
    };
}
