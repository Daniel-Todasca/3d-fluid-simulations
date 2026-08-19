#include <iostream>
#include <iomanip>
#include <algorithm>

#include <GLFW/glfw3.h>

#include "simulation.cpp"

int main() {
    int N = CUBE_SIZE_DEFAULT;

    fsim::FluidCube *cube = new fsim::FluidCube();
    fsim::GridSimulation2d *sim = new fsim::GridSimulation2d(cube);

    int cx = cube->size / 2;
    int cy = cube->size / 2;

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

    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            float distance = sqrt((cx-i) * (cx-i) + (cy-j) * (cy-j));
            if (distance > 2) continue;
            cube->addDensity(i, j, 2.0f / (1.0f + distance));
        }
    }
    
    for (int i = 0; i< N/2; i++) {
        cube->addVelocity(N/2, 3*N/4-i, 0, 1.0f * CUBE_TIMESTEP_DEFAULT);
    }

    while (!glfwWindowShouldClose(window)) {
        glfwPollEvents();
        
        for (int i = 0; i < N; i++) {
            for (int j = 0; j < N; j++) {
                float distance = sqrt((cx-i) * (cx-i) + (cy-j) * (cy-j));
                if (distance > 2) continue;
                cube->addDensity(i, j, 0.01f * CUBE_TIMESTEP_DEFAULT / (1.0f + distance));
            }
        }
    
        for (int i = 0; i < N/2; i++) {
            cube->addVelocity(3*N/4, 3*N/4-i, 0, 0.001f * CUBE_TIMESTEP_DEFAULT);
        }
    
        for (int i = 0; i < N/2; i++) {
            cube->addVelocity(N/4, N/4-i, 0, -0.001f * CUBE_TIMESTEP_DEFAULT);
        }

        sim->step();
    
        glClear(GL_COLOR_BUFFER_BIT);

        float widthPx = 2.0f / cube->size;
        float heightPx = 2.0f / cube->size;

        glBegin(GL_QUADS);
        
        for (int i = 0; i < N/2; i++) {
            cube->addVelocity(N/2, 3*N/4-i, 0, 0.001f * CUBE_TIMESTEP_DEFAULT);
        }

        for (int y = 0; y < cube->size; ++y) {
            for (int x = 0; x < cube->size; ++x) {
                float density = cube->density[cube->indexOf(x, y)];
                float brightness = std::clamp(density / 1.0f, 0.0f, 1.0f);

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