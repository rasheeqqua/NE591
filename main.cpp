#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <iomanip>
#include <algorithm>
#include <chrono>
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

    // Check quadrature weights sum to 2
    double sum = 0.0;
    for (int n = 0; n < N; n++) {
        sum += omega[n];
    }

    if (std::abs(sum - 2.0) > 1e-10) {
        std::cerr << "Warning: Quadrature weights sum to " << sum << " instead of 2.0." << std::endl;
    }

    return true;
}

// Function to implement source iteration algorithm
bool sourceIteration(int N, int I, double sigma_T, double sigma_S, double q, double L,
                     double eps_bar, int k_bar, const std::vector<double>& mu,
                     const std::vector<double>& omega, std::vector<double>& flux,
                     int& iterations, double& final_error) {

    double delta = L / I; // Cell width

    // Initialize the scalar flux to zero
    std::vector<double> phi_old(I, 0.0);
    std::vector<double> phi_new(I, 0.0);

    // Angular flux at cell edges
    std::vector<double> psi_edge(I + 1, 0.0);

    // Angular flux average in cells
    std::vector<double> psi_bar(I, 0.0);

    bool converged = false;
    double epsilon = 0.0;

    // Iteration loop
    for (int k = 1; k <= k_bar; k++) {

        // Set new flux vector to zero
        std::fill(phi_new.begin(), phi_new.end(), 0.0);

        // Loop over all angles
        for (int n = 0; n < N; n++) {

            // Process positive mu (left to right sweep)
            if (mu[n] > 0) {
                // Left boundary condition
                psi_edge[0] = 0.0;

                // Sweep from left to right
                for (int i = 0; i < I; i++) {
                    double source = sigma_S * phi_old[i] + q;
                    double denom = (2.0 * mu[n] / delta) + sigma_T;
                    double numer = source + (2.0 * mu[n] / delta) * psi_edge[i];

                    psi_bar[i] = numer / denom;
                    psi_edge[i+1] = 2.0 * psi_bar[i] - psi_edge[i];

                    // Accumulate contribution to scalar flux
                    phi_new[i] += omega[n] * psi_bar[i];
                }
            }
            // Process negative mu (right to left sweep)
            else {
                // Right boundary condition
                psi_edge[I] = 0.0;

                // Sweep from right to left
                for (int i = I-1; i >= 0; i--) {
                    double source = sigma_S * phi_old[i] + q;
                    double denom = (2.0 * std::abs(mu[n]) / delta) + sigma_T;
                    double numer = source + (2.0 * std::abs(mu[n]) / delta) * psi_edge[i+1];

                    psi_bar[i] = numer / denom;
                    psi_edge[i] = 2.0 * psi_bar[i] - psi_edge[i+1];

                    // Accumulate contribution to scalar flux
                    phi_new[i] += omega[n] * psi_bar[i];
                }
            }
        }

        // Calculate relative error
        epsilon = 0.0;
        for (int i = 0; i < I; i++) {
            // Protect against division by zero
            if (std::abs(phi_old[i]) > 1e-10) {
                double rel_err = std::abs(phi_new[i] / phi_old[i] - 1.0);
                epsilon = std::max(epsilon, rel_err);
            } else {
                if (std::abs(phi_new[i]) > 1e-10) {
                    // If old flux is zero but new isn't, use a large relative error
                    epsilon = std::max(epsilon, 1.0);
                }
            }
        }

        // Check for convergence
        if (epsilon < eps_bar) {
            converged = true;
            iterations = k;
            break;
        }

        // Update old flux for next iteration
        phi_old = phi_new;
    }

    // Store final flux and error
    flux = phi_new;
    final_error = epsilon;
    iterations = converged ? iterations : k_bar;

    return converged;
}

int main() {
    // Start timing
    auto start_time = std::chrono::high_resolution_clock::now();

    std::string input_filename = "../input.txt";
    std::string output_filename = "./output.txt";

    // Open input file
    std::ifstream input_file(input_filename);
    if (!input_file.is_open()) {
        std::cerr << "Error: Unable to open input file " << input_filename << std::endl;
        return 1;
    }

    // Open output file
    std::ofstream output_file(output_filename);
    if (!output_file.is_open()) {
        std::cerr << "Error: Unable to open output file " << output_filename << std::endl;
        return 1;
    }

    // Write header to output file
    output_file << "==========================================" << std::endl;
    output_file << "Neutron Transport Equation Solver" << std::endl;
    output_file << "One-dimensional Slab Geometry, Discrete Ordinates Method" << std::endl;
    output_file << "Author: Hasibul H. Rasheeq" << std::endl;
    output_file << "Affiliation: NC State University" << std::endl;
    output_file << "Date: " << __DATE__ << std::endl;
    output_file << "==========================================" << std::endl << std::endl;

    // Read input data
    int N, I;
    double sigma_T, sigma_S;
    double q, L;
    double eps_bar;
    int k_bar;

    input_file >> N >> I;
    input_file >> sigma_T >> sigma_S;
    input_file >> q >> L;
    input_file >> eps_bar >> k_bar;

    // Read quadrature points and weights
    std::vector<double> mu(N);
    std::vector<double> omega(N);

    for (int n = 0; n < N; n++) {
        input_file >> omega[n] >> mu[n];
    }

    // Validate input data
    if (!validateInput(N, I, sigma_T, sigma_S, q, L, eps_bar, k_bar, mu, omega)) {
        output_file << "Error: Invalid input data" << std::endl;
        return 1;
    }

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

    // Perform source iteration
    std::vector<double> flux(I);
    int iterations;
    double final_error;

    bool converged = sourceIteration(N, I, sigma_T, sigma_S, q, L, eps_bar, k_bar,
                                   mu, omega, flux, iterations, final_error);

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
        output_file << i+1 << "\t" << std::scientific << std::setprecision(6) << flux[i] << std::endl;
    }

    // Calculate and report elapsed time
    auto end_time = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = end_time - start_time;

    output_file << std::endl;
    output_file << "Wall-clock time: " << elapsed.count() << " seconds" << std::endl;

    // Close files
    input_file.close();
    output_file.close();

    return 0;
}