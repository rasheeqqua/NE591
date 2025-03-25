#ifndef VERIFY_POSITIVE_DEFINITE_MATRIX_H
#define VERIFY_POSITIVE_DEFINITE_MATRIX_H

#include <vector>

/**
 * Check if a matrix is symmetric
 * @param A The matrix to check
 * @return true if symmetric, false otherwise
 */
bool isSymmetric(const std::vector<std::vector<double>>& A);

/**
 * Check if a matrix is diagonally dominant
 * A matrix is diagonally dominant if |a_ii| >= sum(|a_ij|) for all j != i
 * @param A The matrix to check
 * @return true if diagonally dominant, false otherwise
 */
bool isDiagonallyDominant(const std::vector<std::vector<double>>& A);

#endif // VERIFY_POSITIVE_DEFINITE_MATRIX_H