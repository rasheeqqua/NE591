#include <mpi.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <iomanip>
#include <algorithm>
#include <string>

// Function to check if input data is valid
bool validateInput(int N, int I, double sigma_T, double sigma_S, double q, double L,
                  double eps_bar, int k_bar, const std::vector<double>& mu,
                  const std::vector<double>& omega) {

    // Check basic requirements
    if (N <= 0 || I <= 0) {
        std::cerr << "Error: Number of angles (N) and cells (I) must be positive." << std::endl;
        return false;
    }

    if (sigma_T < 0 || sigma_S < 0) {
        std::cerr << "Error: Cross sections must be non-negative." << std::endl;
        return false;
    }

    if (L <= 0) {
        std::cerr << "Error: Slab width must be positive." << std::endl;
        return false;
    }

    if (eps_bar <= 0 || k_bar <= 0) {
        std::cerr << "Error: Stopping criterion and max iterations must be positive." << std::endl;
        return false;
    }

    // Check quadrature points and weights
    if (mu.size() != N || omega.size() != N) {
        std::cerr << "Error: Number of quadrature points/weights does not match N." << std::endl;
        return false;
    }

    return true;
}

// Function to implement parallel source iteration algorithm
bool sourceIteration(int N, int I, double sigma_T, double sigma_S, double q, double L,
                     double eps_bar, int k_bar, const std::vector<double>& mu,
                     const std::vector<double>& omega, std::vector<double>& flux,
                     int& iterations, double& final_error, int rank, int size) {

    double delta = L / I; // Cell width

    // Initialize the scalar flux to zero
    std::vector<double> phi_old(I, 0.0);
    std::vector<double> phi_new(I, 0.0);
    std::vector<double> phi_local(I, 0.0); // Local contribution to flux

    // Create the full quadrature set (including both positive and negative mu)
    std::vector<double> all_mu;
    std::vector<double> all_omega;

    for (int n = 0; n < N; n++) {
        all_mu.push_back(mu[n]);     // Positive mu
        all_mu.push_back(-mu[n]);    // Negative mu
        all_omega.push_back(omega[n]);  // Weight for positive mu
        all_omega.push_back(omega[n]);  // Same weight for negative mu
    }

    int total_angles = all_mu.size();

    // Calculate angles per process (assuming N is divisible by size)
    int total_angles_per_proc = total_angles / size;
    int start_angle = rank * total_angles_per_proc;
    int end_angle = (rank + 1) * total_angles_per_proc;

    bool converged = false;
    double epsilon = 0.0;

    // Iteration loop
    int k;
    for (k = 1; k <= k_bar; k++) {
        // Reset local flux
        std::fill(phi_local.begin(), phi_local.end(), 0.0);

        // Loop over assigned angles for this process
        for (int n = start_angle; n < end_angle; n++) {
            std::vector<double> psi_edge(I + 1, 0.0);
            std::vector<double> psi_avg(I, 0.0);

            // For positive mu values (left to right sweep)
            if (all_mu[n] > 0.0) {
                // Left boundary condition
                psi_edge[0] = 0.0;

                // Sweep from left to right
                for (int i = 0; i < I; i++) {
                    // Calculate cell-averaged angular flux
                    double source = sigma_S * phi_old[i] + q;
                    double numer = source + (2.0 * all_mu[n] / delta) * psi_edge[i];
                    double denom = (2.0 * all_mu[n] / delta) + sigma_T;

                    psi_avg[i] = numer / denom;

                    // Calculate edge value using diamond difference
                    psi_edge[i+1] = 2.0 * psi_avg[i] - psi_edge[i];

                    // Accumulate contribution to local scalar flux
                    phi_local[i] += all_omega[n] * psi_avg[i];
                }
            }
            // For negative mu values (right to left sweep)
            else {
                // Right boundary condition
                psi_edge[I] = 0.0;

                // Sweep from right to left
                for (int i = I-1; i >= 0; i--) {
                    // Calculate cell-averaged angular flux
                    double source = sigma_S * phi_old[i] + q;
                    double numer = source + (2.0 * std::abs(all_mu[n]) / delta) * psi_edge[i+1];
                    double denom = (2.0 * std::abs(all_mu[n]) / delta) + sigma_T;

                    psi_avg[i] = numer / denom;

                    // Calculate edge value using diamond difference
                    psi_edge[i] = 2.0 * psi_avg[i] - psi_edge[i+1];

                    // Accumulate contribution to local scalar flux
                    phi_local[i] += all_omega[n] * psi_avg[i];
                }
            }
        }

        // Sum the local contributions across all processes to get the global phi_new
        MPI_Allreduce(phi_local.data(), phi_new.data(), I, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

        // Calculate maximum relative error (done by all processes)
        epsilon = 0.0;
        for (int i = 0; i < I; i++) {
            if (std::abs(phi_old[i]) > 1e-10) {
                double rel_err = std::abs(phi_new[i] / phi_old[i] - 1.0);
                epsilon = std::max(epsilon, rel_err);
            } else if (std::abs(phi_new[i]) > 1e-10) {
                // Large error if old is zero but new isn't
                epsilon = 1.0;
            }
        }

        // Check for convergence
        if (epsilon < eps_bar) {
            converged = true;
            break;
        }

        // Update old flux for next iteration
        phi_old = phi_new;
    }

    // Store the final flux result and iteration info
    flux = phi_new;
    iterations = k;
    final_error = epsilon;

    return converged;
}

int main(int argc, char *argv[]) {
    // Initialize MPI environment
    MPI_Init(&argc, &argv);

    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    // Start timing
    double start_time = MPI_Wtime();

    // Input and output filenames
    std::string input_filename = "input.txt";
    std::string output_filename = "output.txt";

    // Variables to hold problem data
    int N, I;
    double sigma_T, sigma_S;
    double q, L;
    double eps_bar;
    int k_bar;
    std::vector<double> omega;
    std::vector<double> mu;
    bool valid_input = true;

    // Manager process (rank 0) reads the input file
    if (rank == 0) {
        // Open input file
        std::ifstream input_file(input_filename);
        if (!input_file.is_open()) {
            std::cerr << "Error: Unable to open input file " << input_filename << std::endl;
            valid_input = false;
        } else {
            // Read input data
            input_file >> N >> I;
            input_file >> sigma_T >> sigma_S;
            input_file >> q >> L;
            input_file >> eps_bar >> k_bar;

            // Check if N is divisible by the number of processes
            if (N * 2 % size != 0) {
                std::cerr << "Error: Number of angles (2*N = " << 2*N << ") must be divisible by the number of processes (" << size << ")." << std::endl;
                valid_input = false;
            } else {
                // Read quadrature points and weights
                omega.resize(N);
                mu.resize(N);

                for (int n = 0; n < N; n++) {
                    input_file >> omega[n] >> mu[n];
                }

                // Validate input data
                valid_input = validateInput(N, I, sigma_T, sigma_S, q, L, eps_bar, k_bar, mu, omega);
            }

            input_file.close();
        }
    }

    // Broadcast whether input is valid
    MPI_Bcast(&valid_input, 1, MPI_C_BOOL, 0, MPI_COMM_WORLD);

    // If input is invalid, terminate all processes
    if (!valid_input) {
        if (rank == 0) {
            std::ofstream output_file(output_filename);
            if (output_file.is_open()) {
                output_file << "Error: Invalid input data" << std::endl;
                output_file.close();
            }
        }
        MPI_Finalize();
        return 1;
    }

    // Broadcast input parameters to all processes
    MPI_Bcast(&N, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&I, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&sigma_T, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&sigma_S, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&q, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&L, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&eps_bar, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&k_bar, 1, MPI_INT, 0, MPI_COMM_WORLD);

    // Non-root processes need to allocate memory for mu and omega
    if (rank != 0) {
        omega.resize(N);
        mu.resize(N);
    }

    // Broadcast quadrature weights and points
    MPI_Bcast(omega.data(), N, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(mu.data(), N, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    // Manager process writes header and echoes input data
    if (rank == 0) {
        std::ofstream output_file(output_filename);
        if (output_file.is_open()) {
            // Write header
            output_file << "==========================================" << std::endl;
            output_file << "Parallel Neutron Transport Equation Solver" << std::endl;
            output_file << "One-dimensional Slab Geometry, Discrete Ordinates Method" << std::endl;
            output_file << "Author: Hasibul H. Rasheeq" << std::endl;
            output_file << "Affiliation: NC State University" << std::endl;
            output_file << "Date: " << __DATE__ << std::endl;
            output_file << "Number of MPI Processes: " << size << std::endl;
            output_file << "==========================================" << std::endl << std::endl;

            // Echo input data
            output_file << "Input Parameters:" << std::endl;
            output_file << "Number of angles (N): " << N << std::endl;
            output_file << "Number of cells (I): " << I << std::endl;
            output_file << "Total cross section (sigma_T): " << sigma_T << std::endl;
            output_file << "Scattering cross section (sigma_S): " << sigma_S << std::endl;
            output_file << "Fixed source strength (q): " << q << std::endl;
            output_file << "Slab width (L): " << L << std::endl;
            output_file << "Stopping criterion (eps_bar): " << eps_bar << std::endl;
            output_file << "Maximum iterations (k_bar): " << k_bar << std::endl << std::endl;

            output_file << "Quadrature Points and Weights:" << std::endl;
            output_file << "n\tomega_n\t\tmu_n" << std::endl;
            for (int n = 0; n < N; n++) {
                output_file << n+1 << "\t" << omega[n] << "\t" << mu[n] << std::endl;
            }
            output_file << std::endl;

            output_file.close();
        }
    }

    // Perform parallel source iteration
    std::vector<double> flux(I);
    int iterations;
    double final_error;

    bool converged = sourceIteration(N, I, sigma_T, sigma_S, q, L, eps_bar, k_bar,
                                   mu, omega, flux, iterations, final_error, rank, size);

    // End timing
    double end_time = MPI_Wtime();
    double elapsed_time = end_time - start_time;

    // Get the maximum elapsed time across all processes
    double global_elapsed_time;
    MPI_Reduce(&elapsed_time, &global_elapsed_time, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);

    // Manager writes results to output file
    if (rank == 0) {
        std::ofstream output_file(output_filename, std::ios::app);
        if (output_file.is_open()) {
            // Report convergence status
            output_file << "Source Iteration Results:" << std::endl;
            if (converged) {
                output_file << "Iterations converged successfully." << std::endl;
            } else {
                output_file << "Iterations did NOT converge within maximum iterations." << std::endl;
            }

            output_file << "Number of iterations: " << iterations << std::endl;
            output_file << "Final relative error: " << final_error << std::endl << std::endl;

            // Write flux values
            output_file << "Final Scalar Flux Values:" << std::endl;
            output_file << "i\tflux" << std::endl;

            for (int i = 0; i < I; i++) {
                output_file << i+1 << "\t" << std::scientific << std::setprecision(6) << std::setw(14) << flux[i] << std::endl;
            }

            // Report elapsed time
            output_file << std::endl;
            output_file << "Wall-clock time: " << global_elapsed_time << " seconds" << std::endl;

            output_file.close();
        }
    }

    // Finalize MPI environment
    MPI_Finalize();

    return 0;
}