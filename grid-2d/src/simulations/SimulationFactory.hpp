#pragma once

#include "../simulations/SmokeSimulation2d.hpp"
#include "../simulations/WaterSimulation2d.hpp"

#include "../types/Scene.hpp"
#include "../types/FluidCube.hpp"
#include "../types/VoxelizedCube.hpp"

namespace fsim {
    class SimulationFactory {
    public:
        ISimulation* createSimulation(Scene scene) {
            if (scene.type == SMOKE_SIMULATION_2D) {
                return createSmokeSimulation2d(scene);
            }
            else if (scene.type == WATER_SIMULATION_2D) {
                return createWaterSimulation2d(scene);
            }
            else {
                return nullptr;
            }
        }

        SmokeSimulation2d* createSmokeSimulation2d(Scene scene) {
            FluidCube *cube = new FluidCube(
                scene.cubeSize,
                scene.diffusion,
                scene.viscosity,
                scene.timestep
            );

            return new SmokeSimulation2d(cube, scene);
        }

        WaterSimulation2d* createWaterSimulation2d(Scene scene) {
            FluidCube *cube = new VoxelizedCube(
                scene.cubeSize,
                scene.diffusion,
                scene.viscosity,
                scene.timestep
            );

            return new WaterSimulation2d(cube, scene);
        }
    };
}
