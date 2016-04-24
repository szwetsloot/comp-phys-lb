g++-5 -fopenmp -funroll-loops -std=gnu++11 -O3 -DUSE_OMP /usr/lib/libz.dylib cnpy.cpp main.cpp -o comp-phys-lb-worker-gcc
