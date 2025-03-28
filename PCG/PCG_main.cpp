#include "PCG_solver.h"
#include "../matrix_modules/matrix_operations.h"
#include "jacobi_preconditioner.h"
#include <cmath>
#include <iomanip>

// Compute residual vector r = b - Ax
std::vector<double> computeResidualPCG(const std::vector<std::vector<double>>& A,
                                     const std::vector<double>& x,
                                     const std::vector<double>& b) {
    return vector_subtract(b, matrix_vector_product(A, x));
}

// Solve linear system using Preconditioned Conjugate Gradient method
bool solveJacobiPCG(const std::vector<std::vector<double>>& A,
                   const std::vector<double>& b,
                   std::vector<double>& x,
                   double epsilon,
                   int maxIter,
                   int& iterations,
                   double& residualNorm) {
    int n = A.size();

    // Initialize x vector with zeros if not initialized
    if (x.size() != n) {
        x.resize(n, 0.0);
    }

    // Calculate initial residual r_0 = b - A*x_0
    std::vector<double> r = computeResidualPCG(A, x, b);

    // Apply preconditioner to get z_0 = M^(-1) * r_0
    std::vector<double> z = applyJacobiPreconditionerEfficient(A, r);

    // Set initial search direction p_0 = z_0
    std::vector<double> p = z;

    // r_dot_z = (r_0, z_0)
    double r_dot_z = scalar_product(r, z);

    // Calculate initial residual norm
    double r_norm = calculate_L2_norm(r);
    residualNorm = r_norm;

    // Check if initial guess is close enough
    if (r_norm <= epsilon) {
        iterations = 0;
        return true;
    }

    // Main PCG iteration loop
    for (iterations = 1; iterations <= maxIter; ++iterations) {
        // Calculate A*p_k
        std::vector<double> Ap = matrix_vector_product(A, p);

        // Calculate step size alpha_k = (r_k, z_k) / (p_k, A*p_k)
        double p_dot_Ap = scalar_product(p, Ap);
        double alpha = r_dot_z / p_dot_Ap;

        // Update solution x_{k+1} = x_k + alpha_k * p_k
        x = vector_add(x, scalar_multiply(p, alpha));

        // Update residual r_{k+1} = r_k - alpha_k * A*p_k
        r = vector_subtract(r, scalar_multiply(Ap, alpha));

        // Calculate new residual norm
        r_norm = calculate_L2_norm(r);
        residualNorm = r_norm;

        // Check convergence
        if (r_norm <= epsilon) {
            return true;
        }

        // Apply preconditioner to get z_{k+1} = M^(-1) * r_{k+1}
        z = applyJacobiPreconditionerEfficient(A, r);

        // Calculate beta_k = (r_{k+1}, z_{k+1}) / (r_k, z_k)
        double r_dot_z_new = scalar_product(r, z);
        double beta = r_dot_z_new / r_dot_z;

        // Update search direction p_{k+1} = z_{k+1} + beta_k * p_k
        p = vector_add(z, scalar_multiply(p, beta));

        // Update r_dot_z for next iteration
        r_dot_z = r_dot_z_new;
    }

    // If we reached max iterations without convergence
    iterations--;
    return false;
}

// Write PCG solution results to output file
void writePCGResults(std::ofstream& outFile,
                    const std::vector<double>& x,
                    int iterations,
                    double residualNorm,
                    const std::vector<std::vector<double>>& A) {
    // Write specific PCG results
    outFile << "Iterative solution using Preconditioned Conjugate Gradient method" << std::endl;
    outFile << "With Jacobi preconditioner" << std::endl << std::endl;

    // Generate and display the Jacobi preconditioner matrix
    outFile << "Preconditioner matrix M (Jacobi):" << std::endl;
    std::vector<std::vector<double>> M = generateJacobiPreconditioner(A);
    for (std::size_t i = 0; i < M.size(); ++i) {
        for (std::size_t j = 0; j < M[i].size(); ++j) {
            outFile << std::scientific << std::setprecision(4) << M[i][j] << " ";
        }
        outFile << std::endl;
    }
    outFile << std::endl;

    outFile << "Iterations converged in: " << iterations << " iterations" << std::endl << std::endl;

    outFile << "Iterative residual: " << std::scientific << std::setprecision(4) << residualNorm << std::endl << std::endl;

    outFile << "Solution vector, x:" << std::endl;
    for (std::size_t i = 0; i < x.size(); ++i) {
        outFile << std::scientific << std::setprecision(4) << x[i] << " ";
    }
    outFile << std::endl << std::endl;
}