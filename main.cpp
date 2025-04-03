#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <iomanip>
#include <cmath>
#include <chrono>
#include <algorithm>

// Function declarations for different solver methods
#include "lup/lup_solver.h"
#include "point_jacobi/point_jacobi_solver.h"
#include "gauss_seidel/gauss_seidel_solver.h"
#include "sor/sor_solver.h"
//#include "parallel_pj/parallel_pj_solver.h"
//#include "parallel_gs/parallel_gs_solver.h"

// Structure to hold problem parameters
struct ProblemParameters {
    int flag;                 // Solution method flag
    int max_iterations;       // Maximum number of iterations
    double tolerance;         // Convergence tolerance
    double sor_weight;        // SOR relaxation parameter
    double a, b;              // Rectangle dimensions
    int m, n;                 // Grid dimensions
    double D, sigma_a;        // Physical parameters
    std::vector<std::vector<double>> q; // Source distribution
};

// Structure to hold solution results
struct SolutionResults {
    std::vector<std::vector<double>> phi;  // Flux solution
    int iterations;                        // Number of iterations performed
    double final_error;                    // Final error
    double max_residual;                   // Maximum absolute residual
    double execution_time;                 // Execution time in ms
    double memory_usage;                   // Memory usage in MB
};

// Function to read input parameters
bool readInputParameters(const std::string& filename, ProblemParameters& params) {
    std::ifstream inputFile(filename);
    if (!inputFile.is_open()) {
        std::cerr << "Error: Unable to open input file " << filename << std::endl;
        return false;
    }

    // Read solution method flag
    inputFile >> params.flag;

    // Read iteration parameters
    inputFile >> params.max_iterations >> params.tolerance >> params.sor_weight;

    // Read rectangle dimensions
    inputFile >> params.a >> params.b;

    // Read grid dimensions
    inputFile >> params.m >> params.n;

    // Read physical parameters
    inputFile >> params.D >> params.sigma_a;

    // Read source distribution
    params.q.resize(params.m, std::vector<double>(params.n));
    for (int i = 0; i < params.m; i++) {
        for (int j = 0; j < params.n; j++) {
            inputFile >> params.q[i][j];
        }
    }

    inputFile.close();

    // Validate input parameters
    if (params.flag < 0 || params.flag > 5) {
        std::cerr << "Error: Invalid solution method flag" << std::endl;
        return false;
    }

    if (params.max_iterations <= 0 || params.tolerance <= 0) {
        std::cerr << "Error: Invalid iteration parameters" << std::endl;
        return false;
    }

    if (params.a <= 0 || params.b <= 0 || params.m <= 0 || params.n <= 0) {
        std::cerr << "Error: Invalid grid parameters" << std::endl;
        return false;
    }

    if (params.D <= 0 || params.sigma_a < 0) {
        std::cerr << "Error: Invalid physical parameters" << std::endl;
        return false;
    }

    return true;
}

// Function to construct the diffusion matrix
void constructDiffusionMatrix(const ProblemParameters& params,
                             std::vector<std::vector<double>>& A,
                             std::vector<double>& b) {
    int size = params.m * params.n;
    double delta = params.a / (params.m + 1);
    double gamma = params.b / (params.n + 1);

    // Initialize matrix and RHS vector
    A.resize(size, std::vector<double>(size, 0.0));
    b.resize(size, 0.0);

    // Compute matrix coefficients
    double center_coef = 2 * params.D / (delta * delta) +
                         2 * params.D / (gamma * gamma) +
                         params.sigma_a;
    double x_coef = -params.D / (delta * delta);
    double y_coef = -params.D / (gamma * gamma);

    // Fill matrix and RHS vector
    for (int i = 1; i <= params.m; i++) {
        for (int j = 1; j <= params.n; j++) {
            int row = (i - 1) * params.n + (j - 1);

            // Diagonal term
            A[row][row] = center_coef;

            // x-direction neighbors
            if (i > 1) {
                int col = (i - 2) * params.n + (j - 1);
                A[row][col] = x_coef;
            }

            if (i < params.m) {
                int col = i * params.n + (j - 1);
                A[row][col] = x_coef;
            }

            // y-direction neighbors
            if (j > 1) {
                int col = (i - 1) * params.n + (j - 2);
                A[row][col] = y_coef;
            }

            if (j < params.n) {
                int col = (i - 1) * params.n + j;
                A[row][col] = y_coef;
            }

            // RHS vector
            b[row] = params.q[i-1][j-1];
        }
    }
}

// Function to convert linear vector to 2D grid
std::vector<std::vector<double>> convertToGrid(const std::vector<double>& x, int m, int n) {
    std::vector<std::vector<double>> grid(m + 2, std::vector<double>(n + 2, 0.0));

    for (int i = 1; i <= m; i++) {
        for (int j = 1; j <= n; j++) {
            int idx = (i - 1) * n + (j - 1);
            grid[i][j] = x[idx];
        }
    }

    return grid;
}

// Function to calculate maximum residual
double calculateMaxResidual(const ProblemParameters& params,
                          const std::vector<std::vector<double>>& phi) {
    double delta = params.a / (params.m + 1);
    double gamma = params.b / (params.n + 1);
    double max_residual = 0.0;

    for (int i = 1; i <= params.m; i++) {
        for (int j = 1; j <= params.n; j++) {
            double laplacian = (phi[i+1][j] - 2 * phi[i][j] + phi[i-1][j]) / (delta * delta) +
                              (phi[i][j+1] - 2 * phi[i][j] + phi[i][j-1]) / (gamma * gamma);

            double residual = std::abs(-params.D * laplacian + params.sigma_a * phi[i][j] - params.q[i-1][j-1]);
            max_residual = std::max(max_residual, residual);
        }
    }

    return max_residual;
}

// Function to estimate memory usage
double estimateMemoryUsage(const ProblemParameters& params, int method_flag) {
    double memory = 0.0;
    int size = params.m * params.n;

    // Memory for problem parameters
    memory += sizeof(ProblemParameters);

    // Memory for source array
    memory += sizeof(double) * params.m * params.n;

    // Memory for solution array
    memory += sizeof(double) * (params.m + 2) * (params.n + 2);

    // Additional memory based on method
    if (method_flag == 0) {
        // LUP: Full matrix + permutation vector
        memory += sizeof(double) * size * size;
        memory += sizeof(int) * size;
    } else if (method_flag == 1 || method_flag == 4) {
        // Point Jacobi or Parallel PJ: Two solution vectors
        memory += sizeof(double) * size * 2;
    } else if (method_flag == 2 || method_flag == 3 || method_flag == 5) {
        // GS, SOR, or Parallel GS: One solution vector
        memory += sizeof(double) * size;
    }

    // Convert to MB
    return memory / (1024 * 1024);
}

// Function to write output
void writeOutput(const std::string& filename, const ProblemParameters& params, const SolutionResults& results) {
    std::ofstream outputFile(filename);
    if (!outputFile.is_open()) {
        std::cerr << "Error: Unable to open output file " << filename << std::endl;
        return;
    }

    // Write header
    outputFile << "NE 591 - Steady State One-Speed Diffusion Equation Solver" << std::endl;
    outputFile << "Version: 5.0" << std::endl;
    outputFile << "Author: Hasibul H. Rasheeq" << std::endl;
    outputFile << "Date: April 03, 2025" << std::endl;
    outputFile << "----------------------------------------" << std::endl << std::endl;

    // Write solution method
    outputFile << "Solution Method: ";
    switch (params.flag) {
        case 0: outputFile << "LUP Decomposition"; break;
        case 1: outputFile << "Point Jacobi Iteration"; break;
        case 2: outputFile << "Gauss-Seidel Iteration"; break;
        case 3: outputFile << "Successive Over-Relaxation"; break;
        case 4: outputFile << "Parallel Point Jacobi Iteration"; break;
        case 5: outputFile << "Parallel Gauss-Seidel Iteration"; break;
    }
    outputFile << std::endl << std::endl;

    // Write problem parameters
    outputFile << "Problem Parameters:" << std::endl;
    outputFile << "  Grid dimensions: " << params.m << " x " << params.n
              << " (" << params.m * params.n << " unknowns)" << std::endl;
    outputFile << "  Rectangle dimensions: " << params.a << " cm x " << params.b << " cm" << std::endl;

    double delta = params.a / (params.m + 1);
    double gamma = params.b / (params.n + 1);
    outputFile << "  Grid spacing: δ = " << delta << " cm, γ = " << gamma << " cm" << std::endl;

    outputFile << "  Physical parameters: D = " << params.D << " cm, Σa = "
              << params.sigma_a << " cm-1" << std::endl << std::endl;

    // Write source distribution
    outputFile << "Distributed fixed source:" << std::endl;
    for (int i = 0; i < params.m; i++) {
        outputFile << "  ";
        for (int j = 0; j < params.n; j++) {
            outputFile << std::scientific << std::setprecision(2) << params.q[i][j] << " ";
        }
        outputFile << std::endl;
    }
    outputFile << std::endl;

    // Write iteration information
    if (params.flag > 0) {
        outputFile << "Iteration Information:" << std::endl;
        outputFile << "  Maximum iterations allowed: " << params.max_iterations << std::endl;
        outputFile << "  Convergence tolerance: " << std::scientific << std::setprecision(2)
                  << params.tolerance << std::endl << std::endl;
    }

    // Write results
    outputFile << "Results:" << std::endl;
    if (params.flag > 0) {
        outputFile << "  Converged after: " << results.iterations << " iterations" << std::endl;
        outputFile << "  Final error: " << std::scientific << std::setprecision(2)
                  << results.final_error << std::endl;
    }
    outputFile << "  Maximum absolute residual: " << std::scientific << std::setprecision(2)
              << results.max_residual << std::endl << std::endl;

    // Write performance metrics
    outputFile << "Performance Metrics:" << std::endl;
    outputFile << "  Execution time: " << std::fixed << std::setprecision(6)
              << results.execution_time << " ms" << std::endl;
    outputFile << "  Memory usage: " << std::fixed << std::setprecision(6)
              << results.memory_usage << " Mb" << std::endl << std::endl;

    // Write solution
    outputFile << "Scalar Flux Solution:" << std::endl;
    for (std::size_t i = 0; i < results.phi.size(); i++) {
        outputFile << "  ";
        for (std::size_t j = 0; j < results.phi[i].size(); j++) {
            outputFile << std::scientific << std::setprecision(2) << results.phi[i][j] << " ";
        }
        outputFile << std::endl;
    }

    outputFile.close();
}

int main() {
    std::string inputFile = "input.txt";
    std::string outputFile = "output.txt";

    // Read problem parameters
    ProblemParameters params;
    if (!readInputParameters(inputFile, params)) {
        return 1;
    }

    // Construct diffusion matrix and RHS vector
    std::vector<std::vector<double>> A;
    std::vector<double> b;
    constructDiffusionMatrix(params, A, b);

    // Initialize solution results
    SolutionResults results;
    results.phi.resize(params.m + 2, std::vector<double>(params.n + 2, 0.0));

    // Solve the system using the specified method
    auto start_time = std::chrono::high_resolution_clock::now();

    std::vector<double> x(params.m * params.n, 0.0);

    switch (params.flag) {
        case 0:
            // LUP method
            lupSolve(A, b, x);
            results.iterations = 1;
            results.final_error = 0.0;
            break;
        case 1:
            // Point Jacobi method
            pointJacobiSolve(A, b, x, params.max_iterations, params.tolerance,
                             results.iterations, results.final_error);
            break;
        case 2:
            // Gauss-Seidel method
            gaussSeidelSolve(A, b, x, params.max_iterations, params.tolerance,
                            results.iterations, results.final_error);
            break;
        case 3:
            // SOR method
            sorSolve(A, b, x, params.max_iterations, params.tolerance, params.sor_weight,
                     results.iterations, results.final_error);
            break;
        /*case 4:
            // Parallel Point Jacobi method
            parallelPointJacobiSolve(A, b, x, params.max_iterations, params.tolerance,
                                    results.iterations, results.final_error);
            break;
        case 5:
            // Parallel Gauss-Seidel method
            parallelGaussSeidelSolve(A, b, x, params.max_iterations, params.tolerance,
                                    results.iterations, results.final_error);
            break;*/
    }

    auto end_time = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double, std::milli> duration = end_time - start_time;
    results.execution_time = duration.count();

    // Convert solution to grid format
    results.phi = convertToGrid(x, params.m, params.n);

    // Calculate maximum residual
    results.max_residual = calculateMaxResidual(params, results.phi);

    // Estimate memory usage
    results.memory_usage = estimateMemoryUsage(params, params.flag);

    // Write output
    writeOutput(outputFile, params, results);

    return 0;
}