#pragma once

#include "../Constants.hpp"

namespace fsim {
    enum SceneType {
        SMOKE_SIMULATION_2D = 0,
        WATER_SIMULATION_2D = 1,
        FLIP_SIMULATION_2D  = 2
    };

    class Scene {
    public:
        Scene() : Scene(
            WATER_SIMULATION_2D_GRAVITY_DEFAULT, CUBE_TIMESTEP_DEFAULT, 1.0f,
            CUBE_DIFFUSION_DEFAULT, CUBE_VISCOSITY_DEFAULT, GRID_BASED_ITER,
            MIN_DENSITY, CUBE_SIZE_DEFAULT, WATER_SHARPEN_THRESHOLD, NUM_PARTICLES_DEFAULT, 
            FLIP_PARTICLE_RADIUS, SCENE_WIDTH_2D, SCENE_HEIGHT_2D, GRID_BASED_ITER, SMOKE_SIMULATION_2D
        ) { }

        Scene(
            float gravity, float timestep, float overRelaxation,
            float diffusion, float viscosity, int iterations,
            float minDensityClass, float cubeSize, float waterSharpenThreshold, float numParticles,
            float particleRadius, float width, float height, float pushParticlesIter, SceneType type
        ) {
            this->gravity = gravity;
            this->timestep = timestep;
            this->overRelaxation = overRelaxation;
            this->diffusion = diffusion;
            this->viscosity = viscosity;
            this->iterations = iterations;
            this->minDensityClass = minDensityClass;
            this->cubeSize = cubeSize;
            this->waterSharpenThreshold = waterSharpenThreshold;
            this->numParticles = numParticles;
            this->particleRadius = particleRadius;
            this->width = width;
            this->height = height;
            this->pushParticlesIter = pushParticlesIter;
            this->type = type;
            this->flipRatio = 0.9f;
        }

        SceneType type;
        float gravity;
        float timestep;
        float overRelaxation;
        float diffusion;
        float viscosity;
        int iterations;
        float minDensityClass;
        float cubeSize;
        float waterSharpenThreshold;
        float numParticles;
        float particleRadius;
        float width;
        float height;
        float pushParticlesIter;
        float flipRatio;
        float particleCollisionDamping = 0.5f;
    };
}
