#pragma once

#include "../Constants.hpp"

namespace fsim {
    enum SceneType {
        SMOKE_SIMULATION_2D = 0,
        WATER_SIMULATION_2D = 1
    };

    class Scene {
    public:
        Scene() : Scene(
            WATER_SIMULATION_2D_GRAVITY_DEFAULT, CUBE_TIMESTEP_DEFAULT, 1.0f,
            CUBE_DIFFUSION_DEFAULT, CUBE_VISCOSITY_DEFAULT, GRID_BASED_ITER,
            MIN_DENSITY, CUBE_SIZE_DEFAULT, WATER_SHARPEN_THRESHOLD, 0, 
            SMOKE_SIMULATION_2D
        ) { }

        Scene(
            float gravity, float timestep, float overRelaxation,
            float diffusion, float viscosity, float iterations,
            float minDensityClass, float cubeSize, float waterSharpenThreshold, float numParticles,
            SceneType type
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
            this->type = type;
        }

        SceneType type;
        float gravity;
        float timestep;
        float overRelaxation;
        float diffusion;
        float viscosity;
        float iterations;
        float minDensityClass;
        float cubeSize;
        float waterSharpenThreshold;
        float numParticles;
    };
}
