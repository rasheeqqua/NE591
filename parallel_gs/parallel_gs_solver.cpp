// parallel_gs/parallel_gs_solver.cpp
#include "parallel_gs_solver.h"
#include <cmath>
#include <algorithm>
#include <mpi.h>
#include <iostream>

// Include main.cpp structures
struct ProblemParameters;

void parallelGaussSeidelSolve(const ProblemParameters& params,
                              std::vector<std::vector<double>>& local_phi,
                              int& iterations,
                              double& final_error) {
    // Get MPI process info
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    // For square mesh, calculate process grid dimensions
    int sqrt_P = static_cast<int>(std::sqrt(size));

    // Check if sqrt_P is an integer
    if (sqrt_P * sqrt_P != size) {
        if (rank == 0) {
            std::cerr << "Error: Number of processes must be a perfect square for parallel GS" << std::endl;
        }
        return;
    }

    // Get row and column of this process in the process grid
    int proc_row = rank / sqrt_P;
    int proc_col = rank % sqrt_P;

    // Calculate physical parameters
    double delta = params.a / (params.m + 1);
    double center_coef = 2 * params.D / (delta * delta) +
                         2 * params.D / (delta * delta) +
                         params.sigma_a;
    double neighbor_coef = -params.D / (delta * delta);

    // Calculate subdomain size (assuming square grid with m=n)
    int n = params.n;
    int subdomain_size = n / sqrt_P;

    // Calculate global indices for this subdomain
    int global_i_start = proc_row * subdomain_size;
    int global_j_start = proc_col * subdomain_size;

    // Allocate local grid with halo cells (+2 for each dimension)
    int local_size = subdomain_size + 2;
    local_phi.resize(local_size, std::vector<double>(local_size, 0.0));

    // Initialize local source array
    std::vector<std::vector<double>> local_q(subdomain_size, std::vector<double>(subdomain_size, 0.0));
    for (int i = 0; i < subdomain_size; i++) {
        for (int j = 0; j < subdomain_size; j++) {
            int global_i = global_i_start + i;
            int global_j = global_j_start + j;
            if (global_i < params.m && global_j < params.n) {
                local_q[i][j] = params.q[global_i][global_j];
            }
        }
    }

    // Determine neighbor ranks
    int north_rank = (proc_row > 0) ? rank - sqrt_P : MPI_PROC_NULL;
    int south_rank = (proc_row < sqrt_P - 1) ? rank + sqrt_P : MPI_PROC_NULL;
    int west_rank = (proc_col > 0) ? rank - 1 : MPI_PROC_NULL;
    int east_rank = (proc_col < sqrt_P - 1) ? rank + 1 : MPI_PROC_NULL;

    // Allocate buffers for halo exchange
    std::vector<double> send_north(subdomain_size);
    std::vector<double> send_south(subdomain_size);
    std::vector<double> send_west(subdomain_size);
    std::vector<double> send_east(subdomain_size);

    std::vector<double> recv_north(subdomain_size);
    std::vector<double> recv_south(subdomain_size);
    std::vector<double> recv_west(subdomain_size);
    std::vector<double> recv_east(subdomain_size);

    // Iteration loop
    iterations = 0;
    double max_change = 0.0;
    double global_max_change = 0.0;

    do {
        iterations++;
        max_change = 0.0;

        // ---- RED CELLS UPDATE ----

        // Update red cells (i+j is even)
        for (int i = 1; i <= subdomain_size; i++) {
            for (int j = 1; j <= subdomain_size; j++) {
                if ((i + j) % 2 == 0) {  // Red cell
                    double old_value = local_phi[i][j];
                    double sum = local_q[i-1][j-1];

                    // Add contributions from neighbors (using old values)
                    sum += neighbor_coef * (local_phi[i+1][j] + local_phi[i-1][j] +
                                          local_phi[i][j+1] + local_phi[i][j-1]);

                    // Update red cell
                    double new_value = sum / center_coef;
                    local_phi[i][j] = new_value;

                    // Calculate change for convergence check
                    if (old_value != 0.0) {
                        double rel_change = std::abs((new_value - old_value) / old_value);
                        max_change = std::max(max_change, rel_change);
                    }
                }
            }
        }

        // Exchange red cell data with neighbors

        // Prepare red data to send
        for (int j = 1; j <= subdomain_size; j++) {
            // Only send red cells (i+j is even)
            send_north[j-1] = ((1 + j) % 2 == 0) ? local_phi[1][j] : 0.0;
            send_south[j-1] = ((subdomain_size + j) % 2 == 0) ? local_phi[subdomain_size][j] : 0.0;
        }

        for (int i = 1; i <= subdomain_size; i++) {
            // Only send red cells (i+j is even)
            send_west[i-1] = ((i + 1) % 2 == 0) ? local_phi[i][1] : 0.0;
            send_east[i-1] = ((i + subdomain_size) % 2 == 0) ? local_phi[i][subdomain_size] : 0.0;
        }

        // Exchange red halo data
        MPI_Request reqs[8];

        // Send red data to neighbors
        MPI_Isend(send_north.data(), subdomain_size, MPI_DOUBLE, north_rank, 0, MPI_COMM_WORLD, &reqs[0]);
        MPI_Isend(send_south.data(), subdomain_size, MPI_DOUBLE, south_rank, 0, MPI_COMM_WORLD, &reqs[1]);
        MPI_Isend(send_west.data(), subdomain_size, MPI_DOUBLE, west_rank, 0, MPI_COMM_WORLD, &reqs[2]);
        MPI_Isend(send_east.data(), subdomain_size, MPI_DOUBLE, east_rank, 0, MPI_COMM_WORLD, &reqs[3]);

        // Receive red data from neighbors
        MPI_Irecv(recv_south.data(), subdomain_size, MPI_DOUBLE, south_rank, 0, MPI_COMM_WORLD, &reqs[4]);
        MPI_Irecv(recv_north.data(), subdomain_size, MPI_DOUBLE, north_rank, 0, MPI_COMM_WORLD, &reqs[5]);
        MPI_Irecv(recv_east.data(), subdomain_size, MPI_DOUBLE, east_rank, 0, MPI_COMM_WORLD, &reqs[6]);
        MPI_Irecv(recv_west.data(), subdomain_size, MPI_DOUBLE, west_rank, 0, MPI_COMM_WORLD, &reqs[7]);

        // Wait for all red communications to complete
        MPI_Waitall(8, reqs, MPI_STATUSES_IGNORE);

        // Update halo cells with received red data
        if (north_rank != MPI_PROC_NULL) {
            for (int j = 1; j <= subdomain_size; j++) {
                if ((0 + j) % 2 == 0) {  // Red halo cell
                    local_phi[0][j] = recv_north[j-1];
                }
            }
        }

        if (south_rank != MPI_PROC_NULL) {
            for (int j = 1; j <= subdomain_size; j++) {
                if (((subdomain_size+1) + j) % 2 == 0) {  // Red halo cell
                    local_phi[subdomain_size+1][j] = recv_south[j-1];
                }
            }
        }

        if (west_rank != MPI_PROC_NULL) {
            for (int i = 1; i <= subdomain_size; i++) {
                if ((i + 0) % 2 == 0) {  // Red halo cell
                    local_phi[i][0] = recv_west[i-1];
                }
            }
        }

        if (east_rank != MPI_PROC_NULL) {
            for (int i = 1; i <= subdomain_size; i++) {
                if ((i + (subdomain_size+1)) % 2 == 0) {  // Red halo cell
                    local_phi[i][subdomain_size+1] = recv_east[i-1];
                }
            }
        }

        // Apply boundary conditions for red cells
        if (north_rank == MPI_PROC_NULL) {
            for (int j = 0; j <= subdomain_size+1; j++) {
                if ((0 + j) % 2 == 0) {  // Red boundary
                    local_phi[0][j] = 0.0;
                }
            }
        }

        if (south_rank == MPI_PROC_NULL) {
            for (int j = 0; j <= subdomain_size+1; j++) {
                if (((subdomain_size+1) + j) % 2 == 0) {  // Red boundary
                    local_phi[subdomain_size+1][j] = 0.0;
                }
            }
        }

        if (west_rank == MPI_PROC_NULL) {
            for (int i = 0; i <= subdomain_size+1; i++) {
                if ((i + 0) % 2 == 0) {  // Red boundary
                    local_phi[i][0] = 0.0;
                }
            }
        }

        if (east_rank == MPI_PROC_NULL) {
            for (int i = 0; i <= subdomain_size+1; i++) {
                if ((i + (subdomain_size+1)) % 2 == 0) {  // Red boundary
                    local_phi[i][subdomain_size+1] = 0.0;
                }
            }
        }

        // Synchronize all processes after red update
        MPI_Barrier(MPI_COMM_WORLD);

        // ---- BLACK CELLS UPDATE ----

        // Update black cells (i+j is odd)
        for (int i = 1; i <= subdomain_size; i++) {
            for (int j = 1; j <= subdomain_size; j++) {
                if ((i + j) % 2 == 1) {  // Black cell
                    double old_value = local_phi[i][j];
                    double sum = local_q[i-1][j-1];

                    // Add contributions from neighbors
                    // (using new values from red cells, old values from black cells)
                    sum += neighbor_coef * (local_phi[i+1][j] + local_phi[i-1][j] +
                                          local_phi[i][j+1] + local_phi[i][j-1]);

                    // Update black cell
                    double new_value = sum / center_coef;
                    local_phi[i][j] = new_value;

                    // Calculate change for convergence check
                    if (old_value != 0.0) {
                        double rel_change = std::abs((new_value - old_value) / old_value);
                        max_change = std::max(max_change, rel_change);
                    }
                }
            }
        }

        // Exchange black cell data with neighbors

        // Prepare black data to send
        for (int j = 1; j <= subdomain_size; j++) {
            // Only send black cells (i+j is odd)
            send_north[j-1] = ((1 + j) % 2 == 1) ? local_phi[1][j] : 0.0;
            send_south[j-1] = ((subdomain_size + j) % 2 == 1) ? local_phi[subdomain_size][j] : 0.0;
        }

        for (int i = 1; i <= subdomain_size; i++) {
            // Only send black cells (i+j is odd)
            send_west[i-1] = ((i + 1) % 2 == 1) ? local_phi[i][1] : 0.0;
            send_east[i-1] = ((i + subdomain_size) % 2 == 1) ? local_phi[i][subdomain_size] : 0.0;
        }

        // Exchange black halo data

        // Send black data to neighbors
        MPI_Isend(send_north.data(), subdomain_size, MPI_DOUBLE, north_rank, 1, MPI_COMM_WORLD, &reqs[0]);
        MPI_Isend(send_south.data(), subdomain_size, MPI_DOUBLE, south_rank, 1, MPI_COMM_WORLD, &reqs[1]);
        MPI_Isend(send_west.data(), subdomain_size, MPI_DOUBLE, west_rank, 1, MPI_COMM_WORLD, &reqs[2]);
        MPI_Isend(send_east.data(), subdomain_size, MPI_DOUBLE, east_rank, 1, MPI_COMM_WORLD, &reqs[3]);

        // Receive black data from neighbors
        MPI_Irecv(recv_south.data(), subdomain_size, MPI_DOUBLE, south_rank, 1, MPI_COMM_WORLD, &reqs[4]);
        MPI_Irecv(recv_north.data(), subdomain_size, MPI_DOUBLE, north_rank, 1, MPI_COMM_WORLD, &reqs[5]);
        MPI_Irecv(recv_east.data(), subdomain_size, MPI_DOUBLE, east_rank, 1, MPI_COMM_WORLD, &reqs[6]);
        MPI_Irecv(recv_west.data(), subdomain_size, MPI_DOUBLE, west_rank, 1, MPI_COMM_WORLD, &reqs[7]);

        // Wait for all black communications to complete
        MPI_Waitall(8, reqs, MPI_STATUSES_IGNORE);

        // Update halo cells with received black data
        if (north_rank != MPI_PROC_NULL) {
            for (int j = 1; j <= subdomain_size; j++) {
                if ((0 + j) % 2 == 1) {  // Black halo cell
                    local_phi[0][j] = recv_north[j-1];
                }
            }
        }

        if (south_rank != MPI_PROC_NULL) {
            for (int j = 1; j <= subdomain_size; j++) {
                if (((subdomain_size+1) + j) % 2 == 1) {  // Black halo cell
                    local_phi[subdomain_size+1][j] = recv_south[j-1];
                }
            }
        }

        if (west_rank != MPI_PROC_NULL) {
            for (int i = 1; i <= subdomain_size; i++) {
                if ((i + 0) % 2 == 1) {  // Black halo cell
                    local_phi[i][0] = recv_west[i-1];
                }
            }
        }

        if (east_rank != MPI_PROC_NULL) {
            for (int i = 1; i <= subdomain_size; i++) {
                if ((i + (subdomain_size+1)) % 2 == 1) {  // Black halo cell
                    local_phi[i][subdomain_size+1] = recv_east[i-1];
                }
            }
        }

        // Apply boundary conditions for black cells
        if (north_rank == MPI_PROC_NULL) {
            for (int j = 0; j <= subdomain_size+1; j++) {
                if ((0 + j) % 2 == 1) {  // Black boundary
                    local_phi[0][j] = 0.0;
                }
            }
        }

        if (south_rank == MPI_PROC_NULL) {
            for (int j = 0; j <= subdomain_size+1; j++) {
                if (((subdomain_size+1) + j) % 2 == 1) {  // Black boundary
                    local_phi[subdomain_size+1][j] = 0.0;
                }
            }
        }

        if (west_rank == MPI_PROC_NULL) {
            for (int i = 0; i <= subdomain_size+1; i++) {
                if ((i + 0) % 2 == 1) {  // Black boundary
                    local_phi[i][0] = 0.0;
                }
            }
        }

        if (east_rank == MPI_PROC_NULL) {
            for (int i = 0; i <= subdomain_size+1; i++) {
                if ((i + (subdomain_size+1)) % 2 == 1) {  // Black boundary
                    local_phi[i][subdomain_size+1] = 0.0;
                }
            }
        }

        // Find global maximum change for convergence check
        MPI_Allreduce(&max_change, &global_max_change, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);

        // Synchronize all processes after full iteration
        MPI_Barrier(MPI_COMM_WORLD);

    } while (global_max_change > params.tolerance && iterations < params.max_iterations);

    final_error = global_max_change;
}