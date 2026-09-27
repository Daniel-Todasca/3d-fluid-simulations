#pragma once

#include "Scene.hpp"

namespace fsim {
    class SceneFactory {
    public:
        Scene getFlipScene() {
            Scene scene = Scene();
            scene.type = FLIP_SIMULATION_2D;
            scene.gravity = 98.1f;
            scene.timestep = 0.0005f;
            scene.overRelaxation = 1.0f;
            scene.diffusion = 0.0f;
            scene.viscosity = 0.01f;
            scene.iterations = 40;
            scene.minDensityClass = 0.1f;
            scene.cubeSize = 64;
            scene.waterSharpenThreshold = 0.8f;
            scene.numParticles = 1000;
            scene.particleRadius = 1.5f;
            scene.width = 300;
            scene.height = 300;
            scene.pushParticlesIter = 40;
            scene.flipRatio = 0.8f;
            scene.particleCollisionDamping = 0.0f;
            scene.seed = 0;
            return scene;
        }

        Scene getSmokeScene() {
            Scene scene = Scene();
            scene.type = SMOKE_SIMULATION_2D;
            scene.gravity = 98.1f;
            scene.timestep = 0.01f;
            scene.overRelaxation = 1.0f;
            scene.diffusion = 0.0001f;
            scene.viscosity = 0.001f;
            scene.iterations = 4;
            scene.minDensityClass = 0.1f;
            scene.cubeSize = 64;
            scene.waterSharpenThreshold = 0.8f;
            scene.numParticles = 5000;
            scene.particleRadius = 1.5f;
            scene.width = 300;
            scene.height = 300;
            scene.pushParticlesIter = 4;
            scene.flipRatio = 0.8f;
            scene.particleCollisionDamping = 0.0f;
            scene.seed = 0;
            return scene;
        }

        Scene getWaterScene() {
            Scene scene = Scene();
            scene.type = WATER_SIMULATION_2D;
            scene.gravity = 98.1f;
            scene.timestep = 0.0005f;
            scene.overRelaxation = 1.0f;
            scene.diffusion = 0.0f;
            scene.viscosity = 0.01f;
            scene.iterations = 40;
            scene.minDensityClass = 0.1f;
            scene.cubeSize = 64;
            scene.waterSharpenThreshold = 0.8f;
            scene.numParticles = 5000;
            scene.particleRadius = 1.5f;
            scene.width = 300;
            scene.height = 300;
            scene.pushParticlesIter = 4;
            scene.flipRatio = 0.8f;
            scene.particleCollisionDamping = 0.0f;
            scene.seed = 0;
            return scene;
        }
    };
}
