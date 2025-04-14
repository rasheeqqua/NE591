// lu_factorization.h
// LU Factorization solver for linear systems
// April 13, 2025

#ifndef LU_FACTORIZATION_H
#define LU_FACTORIZATION_H

/**
 * Solves a linear system Ax = b using LU factorization
 *
 * @param A Coefficient matrix (n x n)
 * @param b Right-hand side vector
 * @param x Solution vector (output)
 * @param n System size
 * @return true if successful, false if matrix is singular
 */
bool luSolve(double** A, double* b, double* x, int n);

/**
 * Performs LU factorization in-place on matrix A
 *
 * @param A Input matrix, overwritten with L and U factors
 * @param n Matrix size
 * @param perm Permutation array for partial pivoting
 * @return true if successful, false if matrix is singular
 */
bool luFactorize(double** A, int n, int* perm);

/**
 * Solves the system LUx = b given the LU factorization and permutation
 *
 * @param A Matrix containing the LU factorization
 * @param b Right-hand side vector
 * @param x Solution vector (output)
 * @param n System size
 * @param perm Permutation array
 */
void luBacksubstitute(double** A, double* b, double* x, int n, int* perm);

#endif // LU_FACTORIZATION_H