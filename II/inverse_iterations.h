// II/inverse_iterations.h
#ifndef INVERSE_ITERATIONS_H
#define INVERSE_ITERATIONS_H

#include <vector>
#include <fstream>

/**
 * Solve for eigenvalue and eigenvector using Inverse Iteration method
 * @param A Original coefficient matrix
 * @param x0 Initial guess vector (will be normalized)
 * @param lambda Approximate eigenvalue (shift value)
 * @param epsilon Convergence criterion
 * @param maxIter Maximum number of iterations
 * @param x Resulting eigenvector (output)
 * @param iterations Actual number of iterations performed (output)
 * @param error Final error in eigenvector (output)
 * @param eigenvalue Computed eigenvalue of the original matrix (output)
 * @param eigenvalueError Final error in eigenvalue (output)
 * @return true if converged, false if reached max iterations without convergence
 */
bool solveInverseIterations(const std::vector<std::vector<double>>& A,
                           const std::vector<double>& x0,
                           double lambda,
                           double epsilon,
                           int maxIter,
                           std::vector<double>& x,
                           int& iterations,
                           double& error,
                           double& eigenvalue,
                           double& eigenvalueError);

/**
 * Write Inverse Iterations results to output file
 * @param outFile Output file stream
 * @param x Eigenvector
 * @param iterations Number of iterations performed
 * @param error Final error
 * @param lambda Initial approximation of eigenvalue
 */
void writeInverseIterationsResults(std::ofstream& outFile,
                                 const std::vector<double>& x,
                                 int iterations,
                                 double error,
                                 double lambda);

/**
 * Create a matrix (A - λI)
 * @param A Original matrix
 * @param lambda Shift value
 * @return Shifted matrix (A - λI)
 */
std::vector<std::vector<double>> createShiftedMatrix(const std::vector<std::vector<double>>& A,
                                                    double lambda);

#endif // INVERSE_ITERATIONS_H