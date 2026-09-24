#include <iostream>
#include <iomanip>
#include <algorithm>
#include <vector>
#include <cmath>

#include <GLFW/glfw3.h>

#include "types/Scene.hpp"
#include "simulations/WaterSimulation2d.hpp"
#include "simulations/SmokeSImulation2d.hpp"
#include "simulations/SimulationFactory.hpp"

#include <thread>

#define INITIAL_VELOCITY 2.0f
#define FRAME_VELOCITY 0.1f
#define DENSITY_RADIUS 3.0f

int POINT_COUNT = 1;

const int N = CUBE_SIZE_DEFAULT;

float pointsX[]         = { N/4, 3 * N/4 };
float pointsY[]         = { 3*N/4, N/4 };
float initialDensity[]  = { 2.0f, 2.0f };
float frameDensity[]    = { 2.0f, 2.0f };

std::ostream& operator<<(std::ostream& os, fsim::FluidCube& cube);

void drawGridAsTexture(fsim::FluidCube *grid) {
    int N = grid->size;

    std::vector<float> pixels(N * N);
    for (int y = 0; y < N; ++y)
        for (int x = 0; x < N; ++x)
            pixels[y * N + x] = fsim::Clamp(grid->density[grid->indexOf(x, y)], 0.0f, 1.0f);

    glTexImage2D(GL_TEXTURE_2D, 0, GL_LUMINANCE, N, N, 0, GL_LUMINANCE, GL_FLOAT, pixels.data());

    glEnable(GL_TEXTURE_2D);
    glBegin(GL_QUADS);
    glTexCoord2f(0, 0); glVertex2f(-1, -1);
    glTexCoord2f(1, 0); glVertex2f( 1, -1);
    glTexCoord2f(1, 1); glVertex2f( 1,  1);
    glTexCoord2f(0, 1); glVertex2f(-1,  1);
    glEnd();
    glDisable(GL_TEXTURE_2D);
}

void drawWater(fsim::VoxelizedCube *cube) {
    int N = cube->size;

    std::vector<float> pixels(N * N);
    for (int y = 0; y < N; ++y)
        for (int x = 0; x < N; ++x)
            pixels[y * N + x] = (cube->isFluid(x, y)) ? 1.0f : 0.0f;

    glTexImage2D(GL_TEXTURE_2D, 0, GL_LUMINANCE, N, N, 0, GL_LUMINANCE, GL_FLOAT, pixels.data());

    glEnable(GL_TEXTURE_2D);
    glBegin(GL_QUADS);
    glTexCoord2f(0, 0); glVertex2f(-1, -1);
    glTexCoord2f(1, 0); glVertex2f( 1, -1);
    glTexCoord2f(1, 1); glVertex2f( 1,  1);
    glTexCoord2f(0, 1); glVertex2f(-1,  1);
    glEnd();
    glDisable(GL_TEXTURE_2D);
}

void drawParticles(fsim::ParticleCollection *particles, fsim::Scene scene) {
    int N = particles->size;

    glPointSize(3.0f);
    glBegin(GL_POINTS);
    for (int i = 0; i < N; ++i) {
        float x = particles->px[i] / scene.width;
        float y = particles->py[i] / scene.height;
        glVertex2f(2.0f * x - 1.0f, 2.0f * y - 1.0f);
    }
    glEnd();
}

void addSwirlVelocity(fsim::FluidCube *grid, float velocity) {
    int distance = grid->size / 8;
    for (int i = 0; i < grid->size/2; i++) {
        grid->addVelocity(pointsX[0] - distance, pointsY[0] + i, 10, -velocity);
        grid->addVelocity(pointsX[0] + distance, pointsY[0] + i, -10, velocity);
        grid->addVelocity(pointsX[1] - distance, pointsY[1] + i, 10, velocity);
        grid->addVelocity(pointsX[1] + distance, pointsY[1] + i, -10, -velocity);
    }
}

void runSmokeSimulation(GLFWwindow* window, fsim::SmokeSimulation2d *simulation) {
    fsim::FluidCube *grid = simulation->getGrid();

    glfwMakeContextCurrent(window);
    glClear(GL_COLOR_BUFFER_BIT);

    GLuint texture;
    glGenTextures(1, &texture);
    glBindTexture(GL_TEXTURE_2D, texture);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP);

    addSwirlVelocity(grid, INITIAL_VELOCITY);
    for (int p = 0; p < POINT_COUNT; p++) {
        grid->addDensityToCircle(pointsX[p], pointsY[p], DENSITY_RADIUS, initialDensity[p]);
    }

    while (!glfwWindowShouldClose(window)) {
        glfwPollEvents();
        
        for (int p = 0; p < POINT_COUNT; p++) {
            grid->addDensityToCircle(pointsX[p], pointsY[p], DENSITY_RADIUS, frameDensity[p] * CUBE_TIMESTEP_DEFAULT);
        }
        addSwirlVelocity(grid, FRAME_VELOCITY * CUBE_TIMESTEP_DEFAULT);

        simulation->step();

        glClear(GL_COLOR_BUFFER_BIT);

        drawGridAsTexture(grid);

        glfwSwapBuffers(window);
    }

    glDeleteTextures(1, &texture);
}

void runWaterSimulation(GLFWwindow* window, fsim::WaterSimulation2d *simulation) {
    fsim::VoxelizedCube *grid = (fsim::VoxelizedCube*) simulation->getGrid();

    glfwMakeContextCurrent(window);
    glClear(GL_COLOR_BUFFER_BIT);

    GLuint texture;
    glGenTextures(1, &texture);
    glBindTexture(GL_TEXTURE_2D, texture);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP);

    // addSwirlVelocity(grid, INITIAL_VELOCITY);
    for (int p = 0; p < POINT_COUNT; p++) {
        grid->addDensityToCircle(pointsX[p], pointsY[p], DENSITY_RADIUS, initialDensity[p]);
    }

    while (!glfwWindowShouldClose(window)) {
        glfwPollEvents();
        
        for (int p = 0; p < POINT_COUNT; p++) {
            grid->addDensityToCircle(pointsX[p], pointsY[p], DENSITY_RADIUS, frameDensity[p] * CUBE_TIMESTEP_DEFAULT);
        }
        // addSwirlVelocity(grid, FRAME_VELOCITY * CUBE_TIMESTEP_DEFAULT);

        simulation->step();

        glClear(GL_COLOR_BUFFER_BIT);

        drawWater(grid);

        glfwSwapBuffers(window);
    }

    glDeleteTextures(1, &texture);
}

void runFlipSimulation(GLFWwindow* window, fsim::FlipSimulation2d *simulation) {
    fsim::MacGrid *grid = (fsim::MacGrid*) simulation->getGrid();
    fsim::ParticleCollection *particles = simulation->getParticles();

    glfwMakeContextCurrent(window);
    glClear(GL_COLOR_BUFFER_BIT);

    GLuint texture;
    glGenTextures(1, &texture);
    glBindTexture(GL_TEXTURE_2D, texture);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP);

    // addSwirlVelocity(grid, INITIAL_VELOCITY);
    // for (int p = 0; p < POINT_COUNT; p++) {
    //     grid->addDensityToCircle(pointsX[p], pointsY[p], DENSITY_RADIUS, initialDensity[p]);
    // }
    for (int p = 0; p < particles->size; p++) {
        particles->px[p] = 100 + (rand() % 1000) * 0.1f;
        particles->py[p] = 100 + (rand() % 1000) * 0.1f;
    }

    while (!glfwWindowShouldClose(window)) {
        glfwPollEvents();
        
        // for (int p = 0; p < POINT_COUNT; p++) {
        //     grid->addDensityToCircle(pointsX[p], pointsY[p], DENSITY_RADIUS, frameDensity[p] * CUBE_TIMESTEP_DEFAULT);
        // }
        addSwirlVelocity(grid, FRAME_VELOCITY * 1000);

        simulation->step();

        glClear(GL_COLOR_BUFFER_BIT);

        drawParticles(particles, simulation->getScene());

        glfwSwapBuffers(window);
    }

    glDeleteTextures(1, &texture);
}

int main() {

    if (!glfwInit()) {
        std::cerr << "Failed to initialize GLFW\n";
        return 1;
    }

    GLFWwindow* window = glfwCreateWindow(600, 600, "Fluid Simulation", nullptr, nullptr);

    if (!window) {
        std::cerr << "Failed to create GLFW window\n";
        glfwTerminate();
        return 1;
    }

    fsim::Scene scene = fsim::Scene();
    scene.type = fsim::FLIP_SIMULATION_2D;

    if (scene.type == fsim::WATER_SIMULATION_2D) {
        scene.iterations = WATER_PRESSURE_ITER;
    } else 
    if (scene.type == fsim::SMOKE_SIMULATION_2D) {
        scene.iterations = GRID_BASED_ITER;
        scene.diffusion = 0.0001f;
        scene.viscosity = 0.001f;
        scene.timestep = 0.01f;
        POINT_COUNT = 2;
    } else
    if (scene.type == fsim::FLIP_SIMULATION_2D) {
        scene.iterations = WATER_PRESSURE_ITER;
    }
    
    fsim::ISimulation *simulation = fsim::SimulationFactory().createSimulation(scene);

    if (scene.type == fsim::WATER_SIMULATION_2D) {
        runWaterSimulation(window, (fsim::WaterSimulation2d*) simulation);
    } else if (scene.type == fsim::SMOKE_SIMULATION_2D) {
        runSmokeSimulation(window, (fsim::SmokeSimulation2d*) simulation);
    } else if (scene.type == fsim::FLIP_SIMULATION_2D) {
        runFlipSimulation(window, (fsim::FlipSimulation2d*) simulation);
    }

    glfwDestroyWindow(window);
    glfwTerminate();

    delete simulation;

    return 0;
}

std::ostream& operator<<(std::ostream& os, fsim::FluidCube& cube) {
    for (int y = 0; y < cube.size; y++) {
        for (int x = 0; x < cube.size; x++)
            os << std::fixed << std::setprecision(2) << std::setw(7)
               << cube.density[cube.indexOf(x, y)];
        os << '\n';
    }
    
    int centerX = cube.size / 2, centerY = cube.size / 2;
    os << "\nCenter density: " << cube.density[cube.indexOf(centerX, centerY)] << '\n';
    os << "Center velocity: ("
       << cube.Vx[cube.indexOf(centerX, centerY)] << ", "
       << cube.Vy[cube.indexOf(centerX, centerY)] << ")\n";

    return os;
}
