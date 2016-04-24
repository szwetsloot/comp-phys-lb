//
//  main.cpp
//  comp-phys-lb
//
//  Created by Jesse Slim on 21/04/16.
//  Copyright © 2016 Jesse Slim. All rights reserved.
//

#include <iostream>
#include <cstring>
#include <cstdlib>
#include <thread>
#include <cmath>
#include <ctime>
#include <chrono>

#ifdef USE_OMP
    #include <omp.h>
#endif
#include "cnpy.h"

#define for_gridpoints_i for (long i = 0; i < N; i++)
#define for_directions_k for (long k = 0; k < q; k++)
#define for_dimensions_a for (long a = 0; a < dim; a++)
#define for_time_t for (long t = 0; t < num_iter; t++)
#define for_threads_n for (long n = 0; n < num_threads; n++)

enum gridpoint_identity_t: long {
    id_deleted = -1,
    id_interior = 0,
    id_bc_noslip,
    id_bc_inflow,
    id_bc_outflow_dir_0,
    id_bc_outflow_dir_1,
    id_bc_outflow_dir_2,
    id_bc_outflow_dir_3,
    id_bc_outflow_dir_4,
    id_bc_outflow_dir_5,
    id_bc_outflow_dir_6,
};

const long q = 7;
long num_iter = -1;
long N = -1;
const long dim = 2;
long num_threads = -1;

double *u_inlet;
double omega = 0.0;
double *rho_inlet;

double *e;
double *w_e;

long *bb_dir;

double *lattice;
double *f_prev;
double *f_next;
double *f_eq;

double *initial_f;
double *initial_rho;
double *initial_u;

double *rho;
double *u;

long *link_list;
gridpoint_identity_t *lattice_identities;

void *memcpy_double(double *out, void *in, size_t size) {
    return std::memcpy(out, in, sizeof(double) * size);
}

void *memcpy_long(long* out, void *in, size_t size) {
    return std::memcpy(out, in, sizeof(long) * size);
}

bool contains_nans(const double *data, size_t size) {
    for (size_t j = 0; j < size; j++) {
        if (std::isnan(data[j]))
        {
            std::cout << "NaN at: " << j << ", value: " << data[j] << std::endl;
            return true;
        }
    }
    return false;
}

void initialize_load_data(std::string inFile) {
    cnpy::NpyArray lattice_in       = cnpy::npz_load(inFile, "lattice");
    cnpy::NpyArray link_list_in     = cnpy::npz_load(inFile, "link_list");
    cnpy::NpyArray lattice_id_in    = cnpy::npz_load(inFile, "lattice_identities");
    
    cnpy::NpyArray q_in             = cnpy::npz_load(inFile, "q");
    cnpy::NpyArray bb_dir_in        = cnpy::npz_load(inFile, "bb_dir");
    
    cnpy::NpyArray u_inlet_in       = cnpy::npz_load(inFile, "u_inlet");
    cnpy::NpyArray omega_in         = cnpy::npz_load(inFile, "omega");
    cnpy::NpyArray rho_inlet_in     = cnpy::npz_load(inFile, "rho_inlet");
    
    cnpy::NpyArray e_in             = cnpy::npz_load(inFile, "e");
    cnpy::NpyArray w_e_in           = cnpy::npz_load(inFile, "w_e");
    
    cnpy::NpyArray initial_rho_in   = cnpy::npz_load(inFile, "initial_rho");
    cnpy::NpyArray initial_u_in     = cnpy::npz_load(inFile, "initial_u");
    
    cnpy::NpyArray num_iter_in      = cnpy::npz_load(inFile, "num_iter");
    cnpy::NpyArray num_threads_in   = cnpy::npz_load(inFile, "num_threads");
    
    N = lattice_in.shape[0];
    // q = *((long*) q_in.data);
    num_iter = *((long*) num_iter_in.data);
    num_threads = *((long*) num_threads_in.data);
    
#ifdef USE_OMP
    omp_set_num_threads(num_threads);
#endif
    
    // time to allocate arrays
    u_inlet             = new double[N * dim]();
    rho_inlet           = new double[N]();
    
    e                   = new double[q * dim]();
    w_e                 = new double[q]();
    
    bb_dir              = new long[q]();
    
    lattice             = new double[N * dim]();
    f_prev              = new double[N * q]();
    f_next              = new double[N * q]();
    f_eq                = new double[N * q]();
    
    initial_f           = new double[N * q]();
    initial_rho         = new double[N]();
    initial_u           = new double[N * dim]();
    
    rho                 = new double[N]();
    u                   = new double[N * dim]();
    
    link_list           = new long[N * q]();
    lattice_identities  = new gridpoint_identity_t[N]();
    
    // copy data
    
    memcpy_double(lattice, lattice_in.data, N * dim);
    memcpy_long(link_list, link_list_in.data, N * q);
    memcpy_long((long*) lattice_identities, lattice_id_in.data, N);
    
    memcpy_long(bb_dir, bb_dir_in.data, q);
    
    memcpy_double(u_inlet, u_inlet_in.data, N * dim);
    memcpy_double(&omega, omega_in.data, 1);
    memcpy_double(rho_inlet, rho_inlet_in.data, N);
    
    memcpy_double(e, e_in.data, q * dim);
    memcpy_double(w_e, w_e_in.data, q);
    
    memcpy_double(initial_rho, initial_rho_in.data, N);
    memcpy_double(initial_u, initial_u_in.data, N * dim);
    
    memcpy_long(&num_iter, num_iter_in.data, 1);
}

void calc_f_eq_single(double rhoIn, double *uIn, double *fOut){
    double uu = 0.0;
    
    for_dimensions_a {
        uu += uIn[a]*uIn[a];
    }
    
    for_directions_k {
        double ue = 0.0;
        
        for_dimensions_a {
            ue += uIn[a] * e[k * dim + a];
        }
        
        fOut[k] = w_e[k] * rhoIn * (1 + 4*ue - 2*uu + 8*ue*ue);
    }
}

void calc_f_eq(double *rhoIn, double *uIn, double *fOut)
{
    for_gridpoints_i {
        if (lattice_identities[i] == id_deleted)
            continue;
        calc_f_eq_single(rhoIn[i], &uIn[i * dim], &fOut[i * q]);
    }
}

void copy_f(double *fOut, double *fIn) {
    memcpy_double(fOut, fIn, N * q);
}

void clear_f(double *fOut) {
    std::memset(fOut, 0, sizeof(double) * N * q);
}

void initial_conditions() {
    // calculate initial distribution
    calc_f_eq(initial_rho, initial_u, initial_f);
    copy_f(f_prev, initial_f);
    
    if (contains_nans(initial_u, N * dim))
        std::cout << "NaNs in initial u" << std::endl;
    
    if (contains_nans(f_prev, N * q))
        std::cout << "NaNs in initial f" << std::endl;
}

void collision_loop() {
#ifdef USE_OMP
    #pragma omp parallel for
#endif
    for_gridpoints_i {
        if (lattice_identities[i] == id_deleted)
            continue;
        // apply right wall bc
        if (lattice_identities[i] >= id_bc_outflow_dir_0)
        {
            long outflow_dir = lattice_identities[i] - id_bc_outflow_dir_0;
            long link_gp = link_list[i * q + bb_dir[outflow_dir]];
            if (link_gp == -1) {
                std::cout << "Invalid outflow link encountered!" << std::endl;
                exit(-1);
            }
            for_directions_k {
                f_prev[i * q + k] = f_prev[link_gp * q + k];
            }
        }
        
        // apply noslip/bounce-back bc
        if (lattice_identities[i] == id_bc_noslip) {
            double reversed_f[q];
            
            for_directions_k {
                reversed_f[bb_dir[k]] = f_prev[i * q + k];
            }
            
            memcpy_double(&f_prev[i * q], reversed_f, q);
        }
        
        rho[i] = 0.0;
        for_directions_k {
            rho[i] += f_prev[i * q + k];
        }
        
        std::memset(&u[i * dim], 0, sizeof(double) * dim);
        
        if (rho[i] > 0.0 && lattice_identities[i] != id_bc_noslip) {
            for_directions_k {
                for_dimensions_a {
                    u[i * dim + a] += f_prev[i * q + k] * e[k * dim + a] / rho[i];
                }
            }
        }
        
        // left wall boundary condition
        if (lattice_identities[i] == id_bc_inflow) {
            rho[i] = u_inlet[i];
            
            memcpy_double(&u[i * dim], &u_inlet[i * dim], dim);
            memcpy_double(&f_prev[i * q], &initial_f[i * q], q);
        }
        
        calc_f_eq_single(rho[i], &u[i * dim], &f_eq[i * q]);
        
        // for interior points: do collision step
        if (lattice_identities[i] == id_interior) {
            for_directions_k {
                f_prev[i * q + k] -= omega * (f_prev[i * q + k] - f_eq[i * q + k]);
            }
        }
    }
}

void streaming_loop() {
#ifdef USE_OMP
    #pragma omp parallel for
#endif
    for_gridpoints_i {
        if (lattice_identities[i] == id_deleted)
            continue;
        for_directions_k {
            // check if the link in this direction is valid
            long link = link_list[i * q + k];
            if (link > -1) {
                f_next[link * q + k] = f_prev[i * q + k];
            }
        }
    }
}

void do_simulation(bool report = true) {
    for_time_t {
        
        collision_loop();
        streaming_loop();
        
        // finish iteration
        copy_f(f_prev, f_next);
        
        if (report && t % (num_iter / 10) == 0) {
            std::cout << "Simulation progress: " << ((t * 100.0) / num_iter) << "%" << std::endl;
        }
    }
}

int main(int argc, const char * argv[]) {
    
    
    initialize_load_data("in.npz");
    initial_conditions();
    
    std::cout << "Starting simulation with " << N << " grid points and " << num_iter << " iterations";
#ifdef USE_OMP
    std::cout << " with OpenMP";
#endif
    std::cout << std::endl;
    
    auto start_time = std::chrono::steady_clock::now();
    do_simulation();
    auto finish_time = std::chrono::steady_clock::now();
    
    std::chrono::duration<double> diff = finish_time - start_time;
    
    std::cout << "Finished in " << diff.count() << " s" << std::endl;
    
    unsigned int u_shape[] = {static_cast<unsigned int>(N), static_cast<unsigned int>(dim)};
    unsigned int rho_shape[] = {static_cast<unsigned int>(N)};
    unsigned int lattice_shape[] = {static_cast<unsigned int>(N), static_cast<unsigned int>(dim)};
    unsigned int element_shape[] = {1};
    unsigned int link_list_shape[] = {static_cast<unsigned int>(N), static_cast<unsigned int>(q)};
    unsigned int lattice_id_shape[] = {static_cast<unsigned int>(N)};
    unsigned int f_shape[] = {static_cast<unsigned int>(N), static_cast<unsigned int>(q)};
    std::string outFile = "out.npz";
    
    cnpy::npz_save(outFile, "N", &N, element_shape, 1, "w");
    cnpy::npz_save(outFile, "lattice", lattice, lattice_shape, 2, "a");
    cnpy::npz_save(outFile, "u", u, u_shape, 2, "a");
    cnpy::npz_save(outFile, "rho", rho, rho_shape, 1, "a");
    cnpy::npz_save(outFile, "link_list", link_list, link_list_shape, 2, "a");
    cnpy::npz_save(outFile, "lattice_identities", (long*) lattice_identities, lattice_id_shape, 1, "a");
    cnpy::npz_save(outFile, "f_prev", f_prev, f_shape, 2, "a");
    
    return 0;
}
