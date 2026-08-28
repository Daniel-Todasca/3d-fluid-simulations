#pragma once

#define makeFloatArray(size)   new float[(size)]()
#define deleteFloatArray(arr)  delete[] (arr)
#define runForNSteps(n)        for (int _step = 0; _step < (n); _step++)
#define LOG_ENABLED            0
#define LOG(msg)               if (LOG_ENABLED) std::cout << msg << '\n'
