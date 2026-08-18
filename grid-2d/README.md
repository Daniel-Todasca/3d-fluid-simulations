..\w64devkit\bin\g++.exe -o out\main.exe src\main.cpp
g++.exe -o out\main.exe src\main.cpp -I"$env:GLFW\include" -L"$env:GLFW\lib-mingw-w64" -lglfw3 -lopengl32 -lgdi32

