#include "CG_solver.h"
#include "../matrix_modules/matrix_operations.h"
#include <cmath>
#include <iomanip>

// Compute residual vector r = b - Ax
std::vector<double> computeResidualCG(const std::vector<std::vector<double>>& A,
                                     const std::vector<double>& x,
                                     const std::vector<double>& b) {
    return vector_subtract(b, matrix_vector_product(A, x));
}

// Solve linear system using Conjugate Gradient method
bool solveCG(const std::vector<std::vector<double>>& A,
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
    std::vector<double> r = computeResidualCG(A, x, b);

    // Set initial search direction p_0 = r_0
    std::vector<double> p = r;

    // Calculate ||b|| for normalized stopping criterion
    double b_norm = calculate_L2_norm(b);

    // Initial residual norm
    double r_norm = calculate_L2_norm(r);
    residualNorm = r_norm;

    // Check if initial guess is close enough
    if (r_norm / b_norm <= epsilon) {
        iterations = 0;
        return true;
    }

    // Main CG iteration loop
    for (iterations = 1; iterations <= maxIter; ++iterations) {
        // Calculate A*p_k
        std::vector<double> Ap = matrix_vector_product(A, p);

        // Calculate step size alpha_k = (r_k, r_k) / (p_k, A*p_k)
        double r_dot_r = scalar_product(r, r);
        double p_dot_Ap = scalar_product(p, Ap);
        double alpha = r_dot_r / p_dot_Ap;

        // Update solution x_{k+1} = x_k + alpha_k * p_k
        x = vector_add(x, scalar_multiply(p, alpha));

        // Update residual r_{k+1} = r_k - alpha_k * A*p_k
        std::vector<double> r_next = vector_subtract(r, scalar_multiply(Ap, alpha));

        // Calculate new residual norm
        double r_next_norm = calculate_L2_norm(r_next);
        residualNorm = r_next_norm;

        // Check convergence using normalized criterion: ||r_{k+1}|| / ||b|| <= epsilon
        if (r_next_norm / b_norm <= epsilon) {
            return true;
        }

        // Calculate beta_k = (r_{k+1}, r_{k+1}) / (r_k, r_k)
        double beta = scalar_product(r_next, r_next) / r_dot_r;

        // Update search direction p_{k+1} = r_{k+1} + beta_k * p_k
        p = vector_add(r_next, scalar_multiply(p, beta));

        // Update residual for next iteration
        r = r_next;
    }

    // If we reached max iterations without convergence
    iterations--;
    return false;
}

// Write CG solution results to output file
void writeCGResults(std::ofstream& outFile,
                   const std::vector<double>& x,
                   int iterations,
                   double residualNorm) {
    // Write specific CG results
    outFile << "Iterative solution using Conjugate Gradient method" << std::endl << std::endl;

    outFile << "Iterations converged in: " << iterations << " iterations" << std::endl << std::endl;

    outFile << "Iterative residual: " << std::scientific << std::setprecision(4) << residualNorm << std::endl << std::endl;

    outFile << "Last iterative vector, x:" << std::endl;
    for (std::size_t i = 0; i < x.size(); ++i) {
        outFile << std::scientific << std::setprecision(4) << x[i] << " ";
    }
    outFile << std::endl << std::endl;
}