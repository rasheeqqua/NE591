#ifndef JACOBI_PRECONDITIONER_H
#define JACOBI_PRECONDITIONER_H

#include <vector>

/**
 * Generate the Jacobi preconditioner matrix (M = diag(A))
 * @param A Coefficient matrix
 * @return Jacobi preconditioner matrix
 */
std::vector<std::vector<double>> generateJacobiPreconditioner(const std::vector<std::vector<double>>& A);

/**
 * Generate the inverse of Jacobi preconditioner matrix (M⁻¹ = diag(1/A_11, 1/A_22, ...))
 * @param A Coefficient matrix
 * @return Inverse of Jacobi preconditioner matrix
 */
std::vector<std::vector<double>> generateJacobiPreconditionerInverse(const std::vector<std::vector<double>>& A);

/**
 * Apply preconditioner to a vector: M⁻¹ * r
 * @param M_inv Inverse of preconditioner matrix
 * @param r Vector to apply preconditioner to
 * @return Result of M⁻¹ * r
 */
std::vector<double> applyJacobiPreconditioner(const std::vector<std::vector<double>>& M_inv,
                                           const std::vector<double>& r);

/**
 * Apply Jacobi preconditioner efficiently without creating full matrix
 * @param A Original coefficient matrix (used for diagonal elements)
 * @param r Vector to apply preconditioner to
 * @return Result of M⁻¹ * r where M = diag(A)
 */
std::vector<double> applyJacobiPreconditionerEfficient(const std::vector<std::vector<double>>& A,
                                                    const std::vector<double>& r);

#endif // JACOBI_PRECONDITIONER_H