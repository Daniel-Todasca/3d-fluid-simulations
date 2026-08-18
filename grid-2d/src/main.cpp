#include <iostream>
#include <iomanip>

#include "simulation.cpp"

int main() {
    fsim::FluidCube *cube = new fsim::FluidCube(16, 0.0f, 0.0f, 0.1f);
    fsim::GridSimulation2d *sim = new fsim::GridSimulation2d(cube);

    int cx = cube->size / 2;
    int cy = cube->size / 2;

    cube->addDensity(cx, cy, 100.0f);
    cube->addVelocity(cx, cy, 1.0f, 0.5f);

    runForNSteps(10) {
        sim->step();
    }

    for (int y = 0; y < cube->size; y++) {
        for (int x = 0; x < cube->size; x++) {
            std::cout << std::fixed << std::setprecision(2) << std::setw(7)
                      << cube->density[cube->indexOf(x, y)];
        }
        std::cout << '\n';
    }

    std::cout << "\nCenter density: " << cube->density[cube->indexOf(cx, cy)] << '\n';
    std::cout << "Center velocity: ("
              << cube->Vx[cube->indexOf(cx, cy)] << ", "
              << cube->Vy[cube->indexOf(cx, cy)] << ")\n";

    delete sim;
    delete cube;

    return 0;
}