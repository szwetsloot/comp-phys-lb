//
//  main.cpp
//  comp-phys-lb
//
//  Created by Jesse Slim on 21/04/16.
//  Copyright © 2016 Jesse Slim. All rights reserved.
//

#include <iostream>
#include <cmath>
#include <ctime>
#include "cnpy.h"

#define for_gridpoints_i for (int i = 0; i < N; i++)
#define for_directions_k for (int k = 0; k < q; k++)
#define for_dimensions_a for (int a = 0; a < dim; a++)
#define for_time_t for (int t = 0; t < num_iter; t++)

enum gridpoint_identity_t {
    id_interior,
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

const int q = 7;
const int num_iter = 10000;
const int N_x = 200;
const int N_y = 100;
const int N = N_x * N_y;
const int dim = 2;

const double u_inlet[2] = {0.004, 0.0};
const double omega = 1.0 / 0.8;
const double rho_inlet = 1.0;

double e[q][dim] = {};
const double w_e[q] = {1./2., 1./12., 1./12., 1./12., 1./12., 1./12., 1./12.};

const int bb_dir[q] = {0, 4, 5, 6, 1, 2, 3};

double lattice[N][dim] = {};
double f_prev[N][q] = {};
double f_next[N][q] = {};
double f_eq[N][q] = {};

double initial_f[N][q] = {};
double initial_rho[N] = {};
double initial_u[N][dim] = {};

double rho[N] = {};
double u[N][dim] = {};

int link_matrix[N_x][N_y][q][dim] = {};
int link_list[N][q] = {};
gridpoint_identity_t lattice_identities[N] = {};

bool contains_nans(const double *data, unsigned int size) {
    for (int j = 0; j < size; j++) {
        if (std::isnan(data[j]))
        {
            std::cout << "NaN at: " << j << ", value: " << data[j] << std::endl;
            return true;
        }
    }
    return false;
}

void initialize_velocities() {
    double angle_diff = 2 * M_PI / double(q-1);
    for (int j = 1; j < q; j++) {
        e[j][0] = std::cos((j-1) * angle_diff);
        e[j][1] = std::sin((j-1) * angle_diff);
    }
}

void initialize_grid() {
    // set positions of lattice points
    for_gridpoints_i {
        int x = i / N_y;
        int y = i % N_y;
        
        lattice[i][0] = double(x);
        // offset even rows
        if (y % 2 == 0)
            lattice[i][0] += 0.5;
        
        lattice[i][1] = double(y) * 0.5 * std::sqrt(3.0);
    }
    
    // clear link matrix
    int *ptr_link_matrix = &link_matrix[0][0][0][0];
    for(int j = 0; j < N_x*N_y*q*dim; j++) {
        ptr_link_matrix[j] = -1;
    }
    
    // construct link matrix
    for_gridpoints_i {
        int x = i / N_y;
        int y = i % N_y;
        
        // stationary direction
        link_matrix[x][y][0][0] = x;
        link_matrix[x][y][0][1] = y;
        
        // +x direction
        link_matrix[x][y][1][0] = x + 1;
        link_matrix[x][y][1][1] = y;
        
        // up-right direction
        link_matrix[x][y][2][0] = x + 1 - (y % 2);
        link_matrix[x][y][2][1] = y + 1;
        
        // up-left direction
        link_matrix[x][y][3][0] = x - (y % 2);
        link_matrix[x][y][3][1] = y + 1;
        
        // -x direction
        link_matrix[x][y][4][0] = x - 1;
        link_matrix[x][y][4][1] = y;
        
        // down-left direction
        link_matrix[x][y][5][0] = x - (y % 2);
        link_matrix[x][y][5][1] = y - 1;
        
        // down-right direction
        link_matrix[x][y][6][0] = x + 1 - (y % 2);
        link_matrix[x][y][6][1] = y - 1;
        
        for_directions_k {
            bool x_high = link_matrix[x][y][k][0] >= N_x;
            bool x_low  = link_matrix[x][y][k][0] <  0;
            bool y_high = link_matrix[x][y][k][1] >= N_y;
            bool y_low  = link_matrix[x][y][k][1] <  0;
            if(x_high || x_low || y_high || y_low) {
                link_matrix[x][y][k][0] = -1;
                link_matrix[x][y][k][1] = -1;
            }
        }
    }
    
    // flatten link matrix
    for_gridpoints_i {
        int x = i / N_y;
        int y = i % N_y;
        
        for_directions_k {
            if (link_matrix[x][y][k][0] >= 0 && link_matrix[x][y][k][1] >= 0) {
                // we have a valid link
                link_list[i][k] = link_matrix[x][y][k][0] * N_y + link_matrix[x][y][k][1];
            } else {
                // link is invalid
                link_list[i][k] = -1;
            }
        }
    }
    
    // identify points
    for_gridpoints_i {
        // upper wall
        if (lattice[i][1] > (N_y - 1.5) * 0.5 * std::sqrt(3.0)) {
            lattice_identities[i] = id_bc_noslip;
            continue;
        }
        
        // lower wall
        if (lattice[i][1] < 0.5 * 0.5 * std::sqrt(3.0)) {
            lattice_identities[i] = id_bc_noslip;
            continue;
        }
        
        // left wall
        if (lattice[i][0] < 0.75) {
            lattice_identities[i] = id_bc_inflow;
            continue;
        }
        
        // right wall
        if (lattice[i][0] > (N_x - 1.0) - 0.25) {
            lattice_identities[i] = id_bc_outflow_dir_1;
            continue;
        }
        
        lattice_identities[i] = id_interior;
    }
}

void calc_f_eq_single(double rhoIn, double uIn[dim], double fOut[q]){
    double uu = 0.0;
    
    for_dimensions_a {
        uu += uIn[a]*uIn[a];
    }
    
    for_directions_k {
        double ue = 0.0;
        
        for_dimensions_a {
            ue += uIn[a] * e[k][a];
        }
        
        fOut[k] = w_e[k] * rhoIn * (1 + 4*ue - 2*uu + 8*ue*ue);
    }
}

void calc_f_eq(double rhoIn[N], double uIn[N][dim], double fOut[N][q])
{
    for_gridpoints_i {
        calc_f_eq_single(rhoIn[i], uIn[i], fOut[i]);
    }
}

void copy_f(double fOut[N][q], const double fIn[N][q]) {
    memcpy(&fOut[0][0], &fIn[0][0], sizeof(double) * N * q);
}

void clear_f(double fOut[N][q]) {
    memset(&fOut[0][0], 0, sizeof(double) * N * q);
}

void initial_conditions() {
    memset(&initial_f[0][0], 0, sizeof(double) * N * q);
    memset(&initial_u[0][0], 0, sizeof(double) * N * dim);
    
    if (contains_nans(&initial_u[0][0], N * dim))
        std::cout << "NaNs in initial u" << std::endl;
    
    memset(&initial_rho[0], 0, sizeof(double) * N);
    
    for_gridpoints_i {
        if (lattice_identities[i] == id_bc_inflow) {
            for_dimensions_a {
                initial_u[i][a] = u_inlet[a];
            }
        } else {
            for_dimensions_a {
                initial_u[i][a] = 0.0;
            }
        }
        
        if (lattice_identities[i] != id_bc_noslip) {
            initial_rho[i] = rho_inlet;
        } else {
            initial_rho[i] = 0.0;
        }
    }
    
    // calculate initial distribution
    calc_f_eq(initial_rho, initial_u, initial_f);
    copy_f(f_prev, initial_f);
    
    if (contains_nans(&initial_u[0][0], N * dim))
        std::cout << "NaNs in initial u" << std::endl;
    
    if (contains_nans(&f_prev[0][0], N * q))
        std::cout << "NaNs in initial f" << std::endl;
}

void do_simulation() {
    for_time_t {
        // calculation of quantities over the grid + collision step
        for_gridpoints_i {
            // apply right wall bc
            if (lattice_identities[i] >= id_bc_outflow_dir_0)
            {
                int outflow_dir = lattice_identities[i] - id_bc_outflow_dir_0;
                int link_gp = link_list[i][bb_dir[outflow_dir]];
                for_directions_k {
                    f_prev[i][k] = f_prev[link_gp][k];
                }
            }
            
            // apply noslip/bounce-back bc
            if (lattice_identities[i] == id_bc_noslip) {
                double reversed_f[q];
                
                for_directions_k {
                    reversed_f[bb_dir[k]] = f_prev[i][k];
                }
                
                memcpy(f_prev[i], reversed_f, sizeof(double) * q);
            }
            
            rho[i] = 0.0;
            for_directions_k {
                rho[i] += f_prev[i][k];
            }
            
            memset(u[i], 0, sizeof(double) * dim);
            
            if (rho[i] > 0.0 && lattice_identities[i] != id_bc_noslip) {
                for_directions_k {
                    for_dimensions_a {
                        u[i][a] += f_prev[i][k] * e[k][a] / rho[i];
                    }
                }
            }
            
            // left wall boundary condition
            if (lattice_identities[i] == id_bc_inflow) {
                rho[i] = initial_rho[i];
                
                memcpy(u[i], u_inlet, sizeof(double) * dim);
                memcpy(f_prev[i], initial_f[i], sizeof(double) * q);
            }
            
            calc_f_eq_single(rho[i], u[i], f_eq[i]);
            
            // for interior points: do collision step
            if (lattice_identities[i] == id_interior) {
                for_directions_k {
                    f_prev[i][k] -= omega * (f_prev[i][k] - f_eq[i][k]);
                }
            }
        }
        
        // streaming step
        for_gridpoints_i {
            for_directions_k {
                // check if the link in this direction is valid
                if (link_list[i][k] > -1) {
                    f_next[link_list[i][k]][k] = f_prev[i][k];
                }
            }
        }
        
        // finish iteration
        copy_f(f_prev, f_next);
    }
}

int main(int argc, const char * argv[]) {
    double start_time;
    double finish_time;
    
    
    initialize_velocities();
    initialize_grid();
    initial_conditions();
    
    start_time = (double) std::clock() / CLOCKS_PER_SEC;
    do_simulation();
    finish_time = (double) std::clock() / CLOCKS_PER_SEC;
    
    std::cout << "Finished in " << (finish_time - start_time) << " s" << std::endl;
    
    unsigned int u_shape[] = {N, dim};
    unsigned int rho_shape[] = {N};
    unsigned int lattice_shape[] = {N, dim};
    unsigned int element_shape[] = {1};
    unsigned int link_list_shape[] = {N, q};
    unsigned int lattice_id_shape[] = {N};
    std::string outFile = "/Users/jesse/Code/comp-phys-lb/out.npz";
    
    cnpy::npz_save(outFile, "N_x", &N_x, element_shape, 1, "w");
    cnpy::npz_save(outFile, "N_y", &N_y, element_shape, 1, "a");
    cnpy::npz_save(outFile, "lattice", &lattice[0][0], lattice_shape, 2, "a");
    cnpy::npz_save(outFile, "u", &u[0][0], u_shape, 2, "a");
    cnpy::npz_save(outFile, "rho", &rho[0], rho_shape, 1, "a");
    cnpy::npz_save(outFile, "link_list", &link_list[0][0], link_list_shape, 2, "a");
    cnpy::npz_save(outFile, "lattice_identities", (int*) &lattice_identities[0], lattice_id_shape, 1, "a");
    
    return 0;
}
