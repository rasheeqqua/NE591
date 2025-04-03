#ifndef SOR_SOLVER_H
#define SOR_SOLVER_H

#include <vector>

/**
 * Solves a linear system Ax = b using the Successive Over-Relaxation (SOR) method
 *
 * @param A Matrix of coefficients
 * @param b Right-hand side vector
 * @param x Initial guess and final solution vector
 * @param max_iterations Maximum number of iterations allowed
 * @param tolerance Convergence tolerance
 * @param omega Relaxation parameter (1 < omega < 2 for over-relaxation)
 * @param iterations Output parameter for actual number of iterations performed
 * @param final_error Output parameter for final error achieved
 */
void sorSolve(const std::vector<std::vector<double>>& A,
              const std::vector<double>& b,
              std::vector<double>& x,
              int max_iterations,
              double tolerance,
              double omega,
              int& iterations,
              double& final_error);

#endif // SOR_SOLVER_H