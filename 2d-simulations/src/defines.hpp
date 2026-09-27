#pragma once

#define makeFloatArray(size)   new float[(size)]()
#define makeIntArray(size)     new int[(size)]()
#define deleteFloatArray(arr)  delete[] (arr)
#define deleteIntArray(arr)    delete[] (arr)
#define runForNSteps(n)        for (int _step = 0; _step < (n); _step++)
#define LOG_ENABLED            1
#define LOG(msg)               if (LOG_ENABLED) std::cout << msg << '\n'

#define CUBE_SIZE_DEFAULT       64
#define CUBE_DIFFUSION_DEFAULT  0.000f
#define CUBE_VISCOSITY_DEFAULT  0.01f
#define CUBE_TIMESTEP_DEFAULT   0.0005f

#define NUM_PARTICLES_DEFAULT    5000
#define FLIP_PARTICLE_RADIUS     1.5f
#define FLIP_RATIO_DEFAULT       0.8f
#define PARTICLE_COLLISION_DAMPING_DEFAULT 0.0f

#define GRID_BASED_ITER         4
#define WATER_PRESSURE_ITER     40

#define SCENE_WIDTH_2D          300
#define SCENE_HEIGHT_2D         300

#define MIN_DENSITY             0.1f
#define WATER_SHARPEN_THRESHOLD 0.8f

#define SCENE_SEED_DEFAULT      0

#define WATER_SIMULATION_2D_GRAVITY_DEFAULT 98.1
