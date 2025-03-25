#ifndef LUP_SOLVER_H
#define LUP_SOLVER_H

#include <vector>
#include <fstream>

/**
 * Perform LUP factorization on matrix A: PA = LU
 * @param A Input matrix
 * @param L Lower triangular matrix with unit diagonal
 * @param U Upper triangular matrix
 * @param P Permutation matrix
 * @return true if factorization successful, false otherwise
 */
bool lupFactorize(const std::vector<std::vector<double>>& A,
                 std::vector<std::vector<double>>& L,
                 std::vector<std::vector<double>>& U,
                 std::vector<std::vector<double>>& P);

/**
 * Apply permutation matrix P to vector b
 * @param P Permutation matrix
 * @param b Input vector
 * @param Pb Output vector after permutation
 */
void applyPermutationMatrix(const std::vector<std::vector<double>>& P,
                             const std::vector<double>& b,
                             std::vector<double>& Pb);

/**
 * Solve Ly = b by forward substitution
 * @param L Lower triangular matrix with unit diagonal
 * @param b Right-hand side vector
 * @param y Solution vector
 */
void forwardSubstitution(const std::vector<std::vector<double>>& L,
                         const std::vector<double>& b,
                         std::vector<double>& y);

/**
 * Solve Ux = y by back substitution
 * @param U Upper triangular matrix
 * @param y Right-hand side vector
 * @param x Solution vector
 */
void backSubstitution(const std::vector<std::vector<double>>& U,
                      const std::vector<double>& y,
                      std::vector<double>& x);

/**
 * Compute residual vector r = b - Ax
 * @param A Coefficient matrix
 * @param x Solution vector
 * @param b Right-hand side vector
 * @return Residual vector
 */
std::vector<double> computeResidual(const std::vector<std::vector<double>>& A,
                                   const std::vector<double>& x,
                                   const std::vector<double>& b);

/**
 * Find maximum absolute value in a vector
 * @param v Input vector
 * @return Maximum absolute value
 */
double maxAbsoluteValue(const std::vector<double>& v);

/**
 * Solve linear system Ax = b using LUP decomposition
 * @param A Coefficient matrix
 * @param b Right-hand side vector
 * @param x Solution vector (output)
 * @param maxResidual Maximum absolute residual (output)
 * @return true if solution successful, false otherwise
 */
bool solveLUP(const std::vector<std::vector<double>>& A,
             const std::vector<double>& b,
             std::vector<double>& x,
             double& maxResidual);

/**
 * Write LUP solution results to output file
 * @param outFile Output file stream
 * @param x Solution vector
 * @param maxResidual Maximum absolute residual
 */
void writeLUPResults(std::ofstream& outFile,
                    const std::vector<double>& x,
                    double maxResidual);

#endif // LUP_SOLVER_H