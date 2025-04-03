#include "point_jacobi_solver.h"
#include <cmath>
#include <algorithm>

void pointJacobiSolve(const std::vector<std::vector<double>>& A,
                     const std::vector<double>& b,
                     std::vector<double>& x,
                     int max_iterations,
                     double tolerance,
                     int& iterations,
                     double& final_error) {

    int n = A.size();
    std::vector<double> x_new(n, 0.0);
    double error = 0.0;

    // Initialize iterations counter
    iterations = 0;

    // Point Jacobi iteration loop
    for (int iter = 0; iter < max_iterations; ++iter) {

        // Compute new solution
        for (int i = 0; i < n; ++i) {
            double sum = b[i];

            // Sum off-diagonal terms
            for (int j = 0; j < n; ++j) {
                if (i != j) {
                    sum -= A[i][j] * x[j];
                }
            }

            // Compute new value for x_i
            x_new[i] = sum / A[i][i];
        }

        // Compute error as infinity norm of solution difference
        error = 0.0;
        for (int i = 0; i < n; ++i) {
            error = std::max(error, std::abs(x_new[i] - x[i]));
        }

        // Update solution vector
        x = x_new;

        // Increment iteration counter
        iterations++;

        // Check convergence
        if (error < tolerance) {
            break;
        }
    }

    // Set final error
    final_error = error;
}