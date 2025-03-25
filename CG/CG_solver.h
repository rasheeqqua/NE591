#ifndef CG_SOLVER_H
#define CG_SOLVER_H

#include <vector>
#include <fstream>

/**
 * Compute residual vector r = b - Ax
 * @param A Coefficient matrix
 * @param x Solution vector
 * @param b Right-hand side vector
 * @return Residual vector
 */
std::vector<double> computeResidualCG(const std::vector<std::vector<double>>& A,
                                     const std::vector<double>& x,
                                     const std::vector<double>& b);

/**
 * Solve linear system Ax = b using Conjugate Gradient method
 * @param A Coefficient matrix (symmetric positive definite)
 * @param b Right-hand side vector
 * @param x Solution vector (output) - will be initialized with zeros if empty
 * @param epsilon Convergence criterion (stop if ||r||/||b|| <= epsilon)
 * @param maxIter Maximum number of iterations
 * @param iterations Actual number of iterations performed (output)
 * @param residualNorm Final residual norm (output)
 * @return true if converged, false if reached max iterations without convergence
 */
bool solveCG(const std::vector<std::vector<double>>& A,
            const std::vector<double>& b,
            std::vector<double>& x,
            double epsilon,
            int maxIter,
            int& iterations,
            double& residualNorm);

/**
 * Write CG solution results to output file
 * @param outFile Output file stream
 * @param x Solution vector
 * @param iterations Number of iterations performed
 * @param residualNorm Final residual norm
 */
void writeCGResults(std::ofstream& outFile,
                   const std::vector<double>& x,
                   int iterations,
                   double residualNorm);

#endif // CG_SOLVER_H