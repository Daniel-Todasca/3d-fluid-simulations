..\w64devkit\bin\g++.exe -o out\main.exe src\main.cpp
$env:GLFW = "C:\Users\dan3lu\Documents\3d-fluid-simulations\do_not_push\glfw-3.5.1.bin.WIN64"
g++.exe -o out\main.exe src\main.cpp -I"$env:GLFW\include" -L"$env:GLFW\lib-mingw-w64" -lglfw3 -lopengl32 -lgdi32 -static
.\out\main.exe

1. Replace bounce collision on walls with a no-folow condition on solid faces: mark cells as solid, zero any velociy component, push particles a fixed margin back
2. Use Poisson instead of Gauss-Seidel
3. Integrate viscosity
4. Create IGrid
6. New images
7. Better debugging
