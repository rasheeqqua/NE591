// II/inverse_iterations.cpp
#include "inverse_iterations.h"
#include "../matrix_modules/matrix_operations.h"
#include "../LUP/LUP_solver.h"
#include "../PI/power_iterations.h"
#include <cmath>
#include <iomanip>
#include <algorithm>

std::vector<std::vector<double>> createShiftedMatrix(const std::vector<std::vector<double>>& A,
                                                    double lambda) {
    int n = A.size();
    std::vector<std::vector<double>> shiftedA = A;

    // Subtract λI from A
    for (int i = 0; i < n; i++) {
        shiftedA[i][i] -= lambda;
    }

    return shiftedA;
}

bool solveInverseIterations(const std::vector<std::vector<double>>& A,
                           const std::vector<double>& x0,
                           double lambda,
                           double epsilon,
                           int maxIter,
                           std::vector<double>& x,
                           int& iterations,
                           double& error,
                           double& eigenvalue,
                           double& eigenvalueError) {
    int n = A.size();

    // Create the shifted matrix (A - λI)
    std::vector<std::vector<double>> shiftedA = createShiftedMatrix(A, lambda);

    // Perform LUP factorization once (reused in iterations)
    std::vector<std::vector<double>> L(n, std::vector<double>(n, 0.0));
    std::vector<std::vector<double>> U(n, std::vector<double>(n, 0.0));
    std::vector<std::vector<double>> P(n, std::vector<double>(n, 0.0));

    if (!lupFactorize(shiftedA, L, U, P)) {
        return false; // Factorization failed
    }

    // Initialize x with the initial guess x0
    x = x0;

    // Normalize the initial vector using L-infinity norm
    x = normalize_Linf(x);

    std::vector<double> x_prev(n, 0.0);
    eigenvalue = 0.0;
    double eigenvalue_prev = 0.0;
    eigenvalueError = 0.0;

    // Main Inverse Iteration loop
    for (iterations = 1; iterations <= maxIter; ++iterations) {
        // Save previous vector for convergence check
        x_prev = x;

        // Solve (A - λI)y = x_{k-1} using LUP decomposition
        // Step 1: Apply permutation to x_{k-1}
        std::vector<double> Px(n, 0.0);
        applyPermutationMatrix(P, x, Px);

        // Step 2: Solve Ly = Px
        std::vector<double> y(n, 0.0);
        forwardSubstitution(L, Px, y);

        // Step 3: Solve Uz = y
        std::vector<double> z(n, 0.0);
        backSubstitution(U, y, z);

        // Update x_k = z / ||z||_∞
        x = normalize_Linf(z);

        // Compute eigenvalue using Rayleigh Quotient
        eigenvalue_prev = eigenvalue;
        eigenvalue = computeRayleighQuotient(A, x);

        // Calculate eigenvalue error
        if (iterations > 1) {
            eigenvalueError = std::abs(eigenvalue - eigenvalue_prev);
        }

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

void writeInverseIterationsResults(std::ofstream& outFile,
                                 const std::vector<double>& x,
                                 int iterations,
                                 double error,
                                 double lambda) {
    outFile << "Approximate eigenvalue = " << std::fixed << std::setprecision(4) << lambda << std::endl << std::endl;

    outFile << "Iterations converged in: " << iterations << " iterations" << std::endl;
    outFile << "Iterative error in eigenvector: " << std::scientific << std::setprecision(4) << error << std::endl << std::endl;

    outFile << "Last iterate vector, x:" << std::endl;
    for (std::size_t i = 0; i < x.size(); ++i) {
        outFile << std::fixed << std::setprecision(4) << x[i] << " ";
    }
    outFile << std::endl << std::endl;
}