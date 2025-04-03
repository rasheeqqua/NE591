#include "parallel_pj_solver.h"
#include <cmath>
#include <algorithm>
#include <mpi.h>
#include <iostream>

// Include main.cpp structures
struct ProblemParameters;

void parallelPointJacobiSolve(const ProblemParameters& params,
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
            std::cerr << "Error: Number of processes must be a perfect square for parallel PJ" << std::endl;
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
    std::vector<std::vector<double>> new_local_phi(local_size, std::vector<double>(local_size, 0.0));

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

        // Apply Point Jacobi iteration to interior points
        max_change = 0.0;
        for (int i = 1; i <= subdomain_size; i++) {
            for (int j = 1; j <= subdomain_size; j++) {
                // Compute new value using Point Jacobi formula
                double old_value = local_phi[i][j];
                double sum = local_q[i-1][j-1];

                // Add contributions from neighbors
                sum += neighbor_coef * (local_phi[i+1][j] + local_phi[i-1][j] +
                                       local_phi[i][j+1] + local_phi[i][j-1]);

                // New value = solution to diffusion equation
                double new_value = sum / center_coef;
                new_local_phi[i][j] = new_value;

                // Calculate change for convergence check
                if (old_value != 0.0) {
                    double rel_change = std::abs((new_value - old_value) / old_value);
                    max_change = std::max(max_change, rel_change);
                }
            }
        }

        // Exchange data with neighbors to update halos

        // Prepare data to send
        for (int j = 1; j <= subdomain_size; j++) {
            send_north[j-1] = new_local_phi[1][j];          // First row
            send_south[j-1] = new_local_phi[subdomain_size][j]; // Last row
        }

        for (int i = 1; i <= subdomain_size; i++) {
            send_west[i-1] = new_local_phi[i][1];           // First column
            send_east[i-1] = new_local_phi[i][subdomain_size];  // Last column
        }

        // Exchange data with neighbors (non-blocking sends, blocking receives)
        MPI_Request reqs[8];

        // Send data to neighbors
        MPI_Isend(send_north.data(), subdomain_size, MPI_DOUBLE, north_rank, 0, MPI_COMM_WORLD, &reqs[0]);
        MPI_Isend(send_south.data(), subdomain_size, MPI_DOUBLE, south_rank, 0, MPI_COMM_WORLD, &reqs[1]);
        MPI_Isend(send_west.data(), subdomain_size, MPI_DOUBLE, west_rank, 0, MPI_COMM_WORLD, &reqs[2]);
        MPI_Isend(send_east.data(), subdomain_size, MPI_DOUBLE, east_rank, 0, MPI_COMM_WORLD, &reqs[3]);

        // Receive data from neighbors
        MPI_Irecv(recv_south.data(), subdomain_size, MPI_DOUBLE, south_rank, 0, MPI_COMM_WORLD, &reqs[4]);
        MPI_Irecv(recv_north.data(), subdomain_size, MPI_DOUBLE, north_rank, 0, MPI_COMM_WORLD, &reqs[5]);
        MPI_Irecv(recv_east.data(), subdomain_size, MPI_DOUBLE, east_rank, 0, MPI_COMM_WORLD, &reqs[6]);
        MPI_Irecv(recv_west.data(), subdomain_size, MPI_DOUBLE, west_rank, 0, MPI_COMM_WORLD, &reqs[7]);

        // Wait for all communications to complete
        MPI_Waitall(8, reqs, MPI_STATUSES_IGNORE);

        // Update halo cells with received data
        if (north_rank != MPI_PROC_NULL) {
            for (int j = 1; j <= subdomain_size; j++) {
                new_local_phi[0][j] = recv_north[j-1];
            }
        }

        if (south_rank != MPI_PROC_NULL) {
            for (int j = 1; j <= subdomain_size; j++) {
                new_local_phi[subdomain_size+1][j] = recv_south[j-1];
            }
        }

        if (west_rank != MPI_PROC_NULL) {
            for (int i = 1; i <= subdomain_size; i++) {
                new_local_phi[i][0] = recv_west[i-1];
            }
        }

        if (east_rank != MPI_PROC_NULL) {
            for (int i = 1; i <= subdomain_size; i++) {
                new_local_phi[i][subdomain_size+1] = recv_east[i-1];
            }
        }

        // Apply external boundary conditions (if this process has a boundary)
        if (north_rank == MPI_PROC_NULL) {  // Top boundary
            for (int j = 0; j <= subdomain_size+1; j++) {
                new_local_phi[0][j] = 0.0;
            }
        }

        if (south_rank == MPI_PROC_NULL) {  // Bottom boundary
            for (int j = 0; j <= subdomain_size+1; j++) {
                new_local_phi[subdomain_size+1][j] = 0.0;
            }
        }

        if (west_rank == MPI_PROC_NULL) {  // Left boundary
            for (int i = 0; i <= subdomain_size+1; i++) {
                new_local_phi[i][0] = 0.0;
            }
        }

        if (east_rank == MPI_PROC_NULL) {  // Right boundary
            for (int i = 0; i <= subdomain_size+1; i++) {
                new_local_phi[i][subdomain_size+1] = 0.0;
            }
        }

        // Find global maximum change for convergence check
        MPI_Allreduce(&max_change, &global_max_change, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);

        // Copy new solution to old solution for next iteration
        local_phi = new_local_phi;

    } while (global_max_change > params.tolerance && iterations < params.max_iterations);

    final_error = global_max_change;
}