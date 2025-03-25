#ifndef SOR_SOLVER_H
#define SOR_SOLVER_H

#include <vector>
#include <fstream>

/**
 * Calculate L2 norm of a vector
 * @param v Input vector
 * @return L2 norm (Euclidean length)
 */
double calculateL2Norm(const std::vector<double>& v);

/**
 * Compute residual vector r = b - Ax
 * @param A Coefficient matrix
 * @param x Solution vector
 * @param b Right-hand side vector
 * @return Residual vector
 */
std::vector<double> computeResidualSOR(const std::vector<std::vector<double>>& A,
                                      const std::vector<double>& x,
                                      const std::vector<double>& b);

/**
 * Solve linear system Ax = b using SOR method
 * @param A Coefficient matrix
 * @param b Right-hand side vector
 * @param x Solution vector (output) - will be initialized with zeros if empty
 * @param omega Relaxation parameter (1 < omega < 2 for over-relaxation)
 * @param epsilon Convergence criterion for residual norm
 * @param maxIter Maximum number of iterations
 * @param iterations Actual number of iterations performed (output)
 * @param residualNorm Final residual norm (output)
 * @return true if converged, false if reached max iterations without convergence
 */
bool solveSOR(const std::vector<std::vector<double>>& A,
             const std::vector<double>& b,
             std::vector<double>& x,
             double omega,
             double epsilon,
             int maxIter,
             int& iterations,
             double& residualNorm);

/**
 * Write SOR solution results to output file
 * @param outFile Output file stream
 * @param x Solution vector
 * @param iterations Number of iterations performed
 * @param residualNorm Final residual norm
 * @param omega Relaxation parameter used
 */
void writeSORResults(std::ofstream& outFile,
                    const std::vector<double>& x,
                    int iterations,
                    double residualNorm,
                    double omega);

#endif // SOR_SOLVER_H