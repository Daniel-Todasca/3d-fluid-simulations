#include <iostream>
#include <iomanip>
#include <algorithm>
#include <vector>
#include <cmath>

#include <GLFW/glfw3.h>

#include "types/Scene.hpp"
#include "types/SceneFactory.hpp"
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

void drawParticles(fsim::ParticleCollection *particles, const fsim::Scene &scene) {
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

void seedKelvinHelmholtz(fsim::FluidCube *grid, float shear) {
    int N = grid->size;
    int mid = N / 2;
    for (int y = 0; y < N; y++) {
        float vx = (y > mid) ? shear : -shear;
        for (int x = 0; x < N; x++) {
            int i = grid->indexOf(x, y);
            grid->Vx[i] = vx;
            grid->Vy[i] = 0.0f;
            grid->density[i] = (y <= mid) ? 1.0f : 0.0f;
        }
    }
    for (int x = 0; x < N; x++) {
        // superpose several wavelengths so the interface rolls up into many vortices
        float kick = shear * (
              0.20f * std::sin(6.2831853f *  6.0f * x / N)
            + 0.12f * std::sin(6.2831853f * 10.0f * x / N + 1.3f)
            + 0.08f * std::sin(6.2831853f *  3.0f * x / N + 2.7f));
        for (int y = mid - 1; y <= mid + 1; y++)
            grid->Vy[grid->indexOf(x, y)] = kick;
    }
}

void driveKelvinHelmholtz(fsim::FluidCube *grid, float shear, float strength) {
    int N = grid->size;
    int mid = N / 2;
    for (int y = 0; y < N; y++) {
        float target = (y > mid) ? shear : -shear;
        for (int x = 0; x < N; x++) {
            int i = grid->indexOf(x, y);
            grid->Vx[i] += (target - grid->Vx[i]) * strength;
        }
    }
    for (int x = 0; x < N; x++)
        grid->density[grid->indexOf(x, 1)] = 1.0f;
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

    glfwSwapInterval(1); // cap to monitor refresh so the roll-up is watchable

    const float shear = 1.4f;

    seedKelvinHelmholtz(grid, shear);

    while (!glfwWindowShouldClose(window)) {
        glfwPollEvents();

        driveKelvinHelmholtz(grid, shear, 0.008f);

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

    for (int p = 0; p < POINT_COUNT; p++) {
        grid->addDensityToCircle(pointsX[p], pointsY[p], DENSITY_RADIUS, initialDensity[p]);
    }

    while (!glfwWindowShouldClose(window)) {
        glfwPollEvents();
        
        for (int p = 0; p < POINT_COUNT; p++) {
            grid->addDensityToCircle(pointsX[p], pointsY[p], DENSITY_RADIUS, frameDensity[p] * CUBE_TIMESTEP_DEFAULT);
        }
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

    glfwSwapInterval(1); 

    srand(simulation->getScene().seed);
    for (int p = 0; p < particles->size / 4; p++) {
        float x = 100 + (rand() % 1000) * 0.1f;
        float y = 100 + (rand() % 1000) * 0.1f;
        particles->px[4*p] = x;
        particles->py[4*p] = y;
        particles->px[4*p+1] = x;
        particles->py[4*p+1] = y;
        particles->px[4*p+2] = x;
        particles->py[4*p+2] = y;
        particles->px[4*p+3] = x;
        particles->py[4*p+3] = y;
    }

    while (!glfwWindowShouldClose(window)) {
        glfwPollEvents();

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
        scene = fsim::SceneFactory().getWaterScene();
    } else 
    if (scene.type == fsim::SMOKE_SIMULATION_2D) {
        scene = fsim::SceneFactory().getSmokeScene();
        POINT_COUNT = 2;
    } else
    if (scene.type == fsim::FLIP_SIMULATION_2D) {
        scene = fsim::SceneFactory().getFlipScene();
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
