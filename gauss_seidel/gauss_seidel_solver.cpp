#include "gauss_seidel_solver.h"
#include <cmath>
#include <algorithm>

void gaussSeidelSolve(const std::vector<std::vector<double>>& A,
                      const std::vector<double>& b,
                      std::vector<double>& x,
                      int max_iterations,
                      double tolerance,
                      int& iterations,
                      double& final_error) {

    int n = A.size();
    std::vector<double> x_old(n);
    double error = 0.0;

    // Initialize iterations counter
    iterations = 0;

    // Gauss-Seidel iteration loop
    for (int iter = 0; iter < max_iterations; ++iter) {
        // Save current solution for error calculation
        x_old = x;

        // Compute new solution in-place
        for (int i = 0; i < n; ++i) {
            double sum = b[i];

            // Sum using updated values for indices < i
            for (int j = 0; j < i; ++j) {
                sum -= A[i][j] * x[j];
            }

            // Sum using previous iteration values for indices > i
            for (int j = i + 1; j < n; ++j) {
                sum -= A[i][j] * x_old[j];
            }

            // Compute new value for x_i
            x[i] = sum / A[i][i];
        }

        // Compute error as infinity norm of solution difference
        error = 0.0;
        for (int i = 0; i < n; ++i) {
            error = std::max(error, std::abs(x[i] - x_old[i]));
        }

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
