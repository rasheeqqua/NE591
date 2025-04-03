#ifndef GAUSS_SEIDEL_SOLVER_H
#define GAUSS_SEIDEL_SOLVER_H

#include <vector>

/**
 * Solves a linear system Ax = b using the Gauss-Seidel iterative method
 *
 * @param A Matrix of coefficients
 * @param b Right-hand side vector
 * @param x Initial guess and final solution vector
 * @param max_iterations Maximum number of iterations allowed
 * @param tolerance Convergence tolerance
 * @param iterations Output parameter for actual number of iterations performed
 * @param final_error Output parameter for final error achieved
 */
void gaussSeidelSolve(const std::vector<std::vector<double>>& A,
                      const std::vector<double>& b,
                      std::vector<double>& x,
                      int max_iterations,
                      double tolerance,
                      int& iterations,
                      double& final_error);

#endif // GAUSS_SEIDEL_SOLVER_H