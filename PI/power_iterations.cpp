// power_iterations.cpp
#include "power_iterations.h"
#include "../matrix_modules/matrix_operations.h"
#include <cmath>
#include <iomanip>
#include <algorithm>

double calculate_Linf_norm(const std::vector<double>& v) {
    double max_abs = 0.0;
    for (const double& val : v) {
        double abs_val = std::abs(val);
        if (abs_val > max_abs) {
            max_abs = abs_val;
        }
    }
    return max_abs;
}

std::vector<double> normalize_Linf(const std::vector<double>& v) {
    double norm = calculate_Linf_norm(v);
    if (norm == 0.0) {
        return v; // Avoid division by zero
    }
    std::vector<double> result(v.size());
    for (std::size_t i = 0; i < v.size(); i++) {
        result[i] = v[i] / norm;
    }
    return result;
}

double computeRayleighQuotient(const std::vector<std::vector<double>>& A,
                               const std::vector<double>& x) {
    // Rayleigh quotient: λ = (x^T * A * x) / (x^T * x)
    std::vector<double> Ax = matrix_vector_product(A, x);
    double numerator = 0.0;
    double denominator = 0.0;

    for (std::size_t i = 0; i < x.size(); i++) {
        numerator += x[i] * Ax[i];
        denominator += x[i] * x[i];
    }

    if (denominator == 0.0) {
        return 0.0; // Avoid division by zero
    }

    return numerator / denominator;
}

bool solvePowerIterations(const std::vector<std::vector<double>>& A,
                          const std::vector<double>& x0,
                          double epsilon,
                          int maxIter,
                          std::vector<double>& x,
                          int& iterations,
                          double& error,
                          double& eigenvalue,
                          double& eigenvalueError,
                          bool useRayleighQuotient) {
    int n = A.size();

    // Initialize x with the initial guess x0
    x = x0;

    // Normalize the initial vector using L-infinity norm
    x = normalize_Linf(x);

    std::vector<double> x_prev(n, 0.0);
    double eigenvalue_prev = 0.0;
    eigenvalue = 0.0;
    eigenvalueError = 0.0;

    // Main Power Iterations loop
    for (iterations = 1; iterations <= maxIter; ++iterations) {
        // Save previous vector for convergence check
        x_prev = x;

        // Calculate x_k = A * x_{k-1}
        x = matrix_vector_product(A, x);

        // Compute eigenvalue before normalization
        eigenvalue_prev = eigenvalue;

        if (useRayleighQuotient) {
            // Compute eigenvalue using Rayleigh Quotient
            eigenvalue = computeRayleighQuotient(A, x);
        } else {
            // Compute eigenvalue using Power Iteration formula: λ = (A*x)_max / x_max
            // where _max is the position of the maximum absolute value in the vector
            int max_idx = 0;
            double max_abs = 0.0;

            for (int i = 0; i < n; i++) {
                double abs_val = std::abs(x_prev[i]);
                if (abs_val > max_abs) {
                    max_abs = abs_val;
                    max_idx = i;
                }
            }

            // Avoid division by zero
            if (std::abs(x_prev[max_idx]) > 1e-10) {
                eigenvalue = x[max_idx] / x_prev[max_idx];
            }
        }

        // Calculate eigenvalue error
        if (iterations > 1) {
            eigenvalueError = std::abs(eigenvalue - eigenvalue_prev);
        }

        // Normalize the vector using L-infinity norm
        x = normalize_Linf(x);

        // Calculate error as the L-infinity norm of the difference between consecutive normalized vectors
        std::vector<double> diff = vector_subtract(x, x_prev);
        error = calculate_Linf_norm(diff);

        // Check convergence - both eigenvector and eigenvalue must converge
        if (error <= epsilon && eigenvalueError <= epsilon) {
            return true;
        }
    }

    // If we reached max iterations without convergence
    iterations--;
    return false;
}

void writePowerIterationsResults(std::ofstream& outFile,
                                const std::vector<double>& x,
                                int iterations,
                                double error) {
    outFile << "Iterations converged in: " << iterations << " iterations" << std::endl;
    outFile << "Iterative error in eigenvector: " << std::scientific << std::setprecision(4) << error << std::endl << std::endl;

    outFile << "Last iterative vector, x:" << std::endl;
    for (std::size_t i = 0; i < x.size(); ++i) {
        outFile << std::fixed << std::setprecision(4) << x[i] << " ";
    }
    outFile << std::endl << std::endl;
}