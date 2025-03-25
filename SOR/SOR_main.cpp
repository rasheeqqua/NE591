#include "SOR_solver.h"
#include <cmath>
#include <iomanip>

// Calculate L2 norm of a vector
double calculateL2Norm(const std::vector<double>& v) {
    double sum = 0.0;
    for (double val : v) {
        sum += val * val;
    }
    return std::sqrt(sum);
}

// Compute residual vector r = b - Ax
std::vector<double> computeResidualSOR(const std::vector<std::vector<double>>& A,
                                      const std::vector<double>& x,
                                      const std::vector<double>& b) {
    int n = A.size();
    std::vector<double> r(n, 0.0);

    // r = b - A*x
    for (int i = 0; i < n; ++i) {
        r[i] = b[i];
        for (int j = 0; j < n; ++j) {
            r[i] -= A[i][j] * x[j];
        }
    }

    return r;
}

// Solve linear system Ax = b using SOR method
bool solveSOR(const std::vector<std::vector<double>>& A,
             const std::vector<double>& b,
             std::vector<double>& x,
             double omega,
             double epsilon,
             int maxIter,
             int& iterations,
             double& residualNorm) {
    int n = A.size();

    // Initialize x vector with zeros if not initialized
    if (x.size() != n) {
        x.resize(n, 0.0);
    }

    // Create temporary vector for the next iteration
    std::vector<double> x_new(n, 0.0);

    // Main SOR iteration loop
    iterations = 0;
    bool converged = false;

    while (iterations < maxIter && !converged) {
        // One SOR iteration
        for (int i = 0; i < n; ++i) {
            double sigma = 0.0;

            // Sum for previously updated components (j < i)
            for (int j = 0; j < i; ++j) {
                sigma += A[i][j] * x_new[j];
            }

            // Sum for not yet updated components (j > i)
            for (int j = i + 1; j < n; ++j) {
                sigma += A[i][j] * x[j];
            }

            // SOR update formula
            x_new[i] = (1.0 - omega) * x[i] + (omega / A[i][i]) * (b[i] - sigma);
        }

        // Compute residual and check for convergence
        std::vector<double> residual = computeResidualSOR(A, x_new, b);
        residualNorm = calculateL2Norm(residual);

        if (residualNorm < epsilon) {
            converged = true;
        }

        // Update x for next iteration
        x = x_new;

        // Increment iteration counter
        iterations++;
    }

    return converged;
}

// Write SOR solution results to output file
void writeSORResults(std::ofstream& outFile,
                    const std::vector<double>& x,
                    int iterations,
                    double residualNorm,
                    double omega) {
    // Write specific SOR results
    outFile << "Iterative solution using SOR method" << std::endl;
    outFile << "With relaxation parameter omega = " << std::fixed << std::setprecision(2) << omega << std::endl << std::endl;

    outFile << "Iterations performed: " << iterations << std::endl << std::endl;

    outFile << "Iterative residual: " << std::scientific << std::setprecision(4) << residualNorm << std::endl << std::endl;

    outFile << "Solution vector, x:" << std::endl;
    for (std::size_t i = 0; i < x.size(); ++i) {
        outFile << std::scientific << std::setprecision(4) << x[i] << " ";
    }
    outFile << std::endl << std::endl;
}