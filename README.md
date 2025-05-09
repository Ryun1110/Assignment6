# Assignment6
 Assignment Overview
This repository implements a simple software rasterization pipeline and three different shading techniques on a tessellated sphere:

Flat Shading
Shading computed once per triangle at its centroid (Q1.cpp).

Gouraud Shading
Per-vertex lighting in world space with linear interpolation across the triangle (Q2.cpp).

Phong Shading
Per-pixel lighting in world space using interpolated world positions and normals (Q3.cpp).

Material and light properties are identical for all techniques:

Ambient coefficient ka = (0, 1, 0)

Diffuse coefficient kd = (0, 0.5, 0)

Specular coefficient ks = (0.5, 0.5, 0.5)

Shininess p = 32

Ambient light intensity Ia = 0.2

Point light at (−4, 4, −3), unit white intensity

Gamma correction γ = 2.2

 Repository Structure
├── CG_hw6.pdf           # Assignment specification
├── Q1.cpp               # Flat shading implementation
├── Q2.cpp               # Gouraud shading implementation
├── Q3.cpp               # Phong shading implementation
├── Q?_output.ppm        # Generated PPM images (after running)
├── Q?_output.bmp        # Generated BMP images (Windows only)
└── README.md            # This file
 Prerequisites
A C++17-compatible compiler (e.g. g++, clang++)

Standard C++ library

Windows-only: <windows.h> and ShellApi.h for automatic BMP display (optional on non-Windows platforms)

 Build & Run
Compile


# Flat shading
g++ -std=c++17 Q1.cpp -o flat_shading

# Gouraud shading
g++ -std=c++17 Q2.cpp -o gouraud_shading

# Phong shading
g++ -std=c++17 Q3.cpp -o phong_shading
Execute
Each executable writes:

Q?_output.ppm (portable pixmap)

Q?_output.bmp (Windows-only BMP)

./flat_shading      # produces Q1_output.ppm & Q1_output.bmp
./gouraud_shading   # produces Q2_output.ppm & Q2_output.bmp
./phong_shading     # produces Q3_output.ppm & Q3_output.bmp
View the Results

PPM: open with any image viewer that supports PPM (e.g. GIMP, ImageMagick).

BMP (Windows): automatically opened via ShellExecuteA.

 Sample Outputs
Flat Shading	Gouraud Shading	Phong Shading

Replace the above with your actual rendered screenshots in docs/.
