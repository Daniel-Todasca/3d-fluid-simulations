#include <iostream>
#include <iomanip>
#include <algorithm>

#include <GLFW/glfw3.h>

#include "simulation.cpp"

int main() {
    int N = 16;

    fsim::FluidCube *cube = new fsim::FluidCube(N, 0.0f, 0.0f, 0.03f);
    fsim::GridSimulation2d *sim = new fsim::GridSimulation2d(cube);

    for (int i = 0; i < N*N; i++) {
        cube->dye[i] = 0;
        cube->density[i] = 0;
        cube->Vx[i] = 0;
        cube->Vy[i] = 0;
        cube->Vx0[i] = 0;
        cube->Vy0[i] = 0;
    }

    int cx = cube->size / 2;
    int cy = cube->size / 2;

    cube->addDensity(cx, cy, 10000.0f);
    cube->addVelocity(cx, cy, 0.0f, 20.0f);

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

    if (!glfwInit()) {
        std::cerr << "Failed to initialize GLFW\n";
        return 1;
    }

    GLFWwindow* window =
        glfwCreateWindow(300, 300, "Fluid Simulation", nullptr, nullptr);

    if (!window) {
        std::cerr << "Failed to create GLFW window\n";
        glfwTerminate();
        return 1;
    }

    glClear(GL_COLOR_BUFFER_BIT);
    glfwMakeContextCurrent(window);

    while (!glfwWindowShouldClose(window)) {
        glfwPollEvents();

        sim->step();
    
        glClear(GL_COLOR_BUFFER_BIT);

        float widthPx = 2.0f / cube->size;
        float heightPx = 2.0f / cube->size;

        glBegin(GL_QUADS);

        for (int y = 0; y < cube->size; ++y) {
            for (int x = 0; x < cube->size; ++x) {
                float density = cube->density[cube->indexOf(x, y)];
                float brightness = std::clamp(density / 100.0f, 0.0f, 1.0f);

                glColor3f(brightness, brightness, brightness);

                float x0 = -1.0f + x * widthPx;
                float y0 = -1.0f + y * heightPx;
                float x1 = x0 + widthPx;
                float y1 = y0 + heightPx;

                glVertex2f(x0, y0);
                glVertex2f(x1, y0);
                glVertex2f(x1, y1);
                glVertex2f(x0, y1);
            }
        }

        glEnd();

        glfwSwapBuffers(window);
    }

    glfwDestroyWindow(window);
    glfwTerminate();

    delete sim;
    delete cube;

    return 0;
}