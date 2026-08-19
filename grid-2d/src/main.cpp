#include <iostream>
#include <iomanip>
#include <algorithm>
#include <vector>
#include <cmath>

#include <GLFW/glfw3.h>

#include "simulation.cpp"

#define INITIAL_VELOCITY 2.0f
#define FRAME_VELOCITY 0.1f
#define DENSITY_RADIUS 3.0f
#define POINT_COUNT 2

const int N = CUBE_SIZE_DEFAULT;

float pointsX[]         = { N/4, 3 * N/4 };
float pointsY[]         = { 3*N/4, N/4 };
float initialDensity[]  = { 2.0f, 2.0f };
float frameDensity[]    = { 0.5f, 0.5f };

std::ostream& operator<<(std::ostream& os, fsim::FluidCube& cube);

void drawGridAsTexture(fsim::FluidCube *grid) {
    int N = grid->size;

    std::vector<float> pixels(N * N);
    for (int y = 0; y < N; ++y)
        for (int x = 0; x < N; ++x)
            pixels[y * N + x] = std::clamp(grid->density[grid->indexOf(x, y)], 0.0f, 1.0f);

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

void addSwirlVelocity(fsim::FluidCube *grid, float velocity) {
    int distance = grid->size / 8;
    for (int i = 0; i < grid->size/2; i++) {
        grid->addVelocity(pointsX[0] - distance, pointsY[0] + i, 0, -velocity);
        grid->addVelocity(pointsX[0] + distance, pointsY[0] + i, 0, velocity);
        grid->addVelocity(pointsX[1] - distance, pointsY[1] + i, 0, velocity);
        grid->addVelocity(pointsX[1] + distance, pointsY[1] + i, 0, -velocity);
    }
}

void runSmokeSimulation(GLFWwindow* window, fsim::GridSimulation2d *simulation) {

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

int main() {

    fsim::FluidCube *cube = new fsim::FluidCube();
    fsim::GridSimulation2d *simulation = new fsim::GridSimulation2d(cube);

    std::cout << *cube;

    if (!glfwInit()) {
        std::cerr << "Failed to initialize GLFW\n";
        return 1;
    }

    GLFWwindow* window = glfwCreateWindow(300, 300, "Fluid Simulation", nullptr, nullptr);

    if (!window) {
        std::cerr << "Failed to create GLFW window\n";
        glfwTerminate();
        return 1;
    }

    runSmokeSimulation(window, simulation);

    glfwDestroyWindow(window);
    glfwTerminate();

    delete simulation;
    delete cube;

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
