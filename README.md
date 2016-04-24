# Lattice Boltzmann simulation for Computational Physics project #3


This folder contains assignment 3 from Computational transport phenomena.

Students:
Duco van Holte tot Echten
Jesse Slim
Sybren Zwetsloot

## File directory
The interesting files in this project are:
### Zwetsloot - Attempt 3.ipynb	
Simulates flow around a cylinder in a rectangular d2q9 domain
### Movie.mp4
Results of flow around the cylinder. 25000 Data points. In the last 20 seconds a vortex sheet can be seen behind the cylinder
### LBJesseInletOutletBC.ipynb
Vectorized approach to solving the lattice Boltzmann problem on a triangular lattice using a neighbour list
In this case the flow through a pipe is simulated
### comp-phys-lb/comp-phys-lb-worker/main.cpp
Translation of the Python algorithm in the file above to C for maximum speed (20-40x speedup). Input and output is handled using by numpy archives
Can be compiled using compile.sh on OS X with GCC 5.2 with OpenMP (or using the Xcode project, but then no OpenMP)
### C-*.ipynb
Simulations of different types of domains using the C program (especially the bent pipe one is nice)

