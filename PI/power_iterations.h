// power_iterations.h
#ifndef POWER_ITERATIONS_H
#define POWER_ITERATIONS_H

#include <vector>
#include <fstream>

/**
 * Solve for fundamental eigenvalue and eigenvector using Power Iterations method
 * @param A Coefficient matrix
 * @param x0 Initial guess vector (will be normalized)
 * @param epsilon Convergence criterion
 * @param maxIter Maximum number of iterations
 * @param x Resulting eigenvector (output)
 * @param iterations Actual number of iterations performed (output)
 * @param error Final error in eigenvector (output)
 * @param eigenvalue Computed eigenvalue (output)
 * @param eigenvalueError Final error in eigenvalue (output)
 * @param useRayleighQuotient Flag to select eigenvalue computation method
 * @return true if converged, false if reached max iterations without convergence
 */
bool solvePowerIterations(const std::vector<std::vector<double>>& A,
                          const std::vector<double>& x0,
                          double epsilon,
                          int maxIter,
                          std::vector<double>& x,
                          int& iterations,
                          double& error,
                          double& eigenvalue,
                          double& eigenvalueError,
                          bool useRayleighQuotient);

/**
 * Write Power Iterations results to output file
 * @param outFile Output file stream
 * @param x Eigenvector
 * @param iterations Number of iterations performed
 * @param error Final error
 */
void writePowerIterationsResults(std::ofstream& outFile,
                                const std::vector<double>& x,
                                int iterations,
                                double error);

/**
 * Calculate the L-infinity norm of a vector (maximum absolute value)
 * @param v Vector
 * @return L-infinity norm
 */
double calculate_Linf_norm(const std::vector<double>& v);

/**
 * Normalize a vector using L-infinity norm
 * @param v Vector to normalize
 * @return Normalized vector
 */
std::vector<double> normalize_Linf(const std::vector<double>& v);

/**
 * Compute eigenvalue using Rayleigh Quotient
 * @param A Matrix
 * @param x Eigenvector
 * @return Computed eigenvalue
 */
double computeRayleighQuotient(const std::vector<std::vector<double>>& A,
                               const std::vector<double>& x);

#endif // POWER_ITERATIONS_H