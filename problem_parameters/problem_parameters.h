#ifndef PROBLEM_PARAMETERS_H
#define PROBLEM_PARAMETERS_H

#include <vector>

/**
 * Structure to hold problem parameters for diffusion equation solver
 */
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

#endif // PROBLEM_PARAMETERS_H