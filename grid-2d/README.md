..\w64devkit\bin\g++.exe -o out\main.exe src\main.cpp
$env:GLFW = "C:\Users\DanielTodasca\Documents\cpp\3d-fluid-simulations\do_not_push\glfw-3.5.1.bin.WIN64"
g++.exe -o out\main.exe src\main.cpp -I"$env:GLFW\include" -L"$env:GLFW\lib-mingw-w64" -lglfw3 -lopengl32 -lgdi32 -static
.\out\main.exe
