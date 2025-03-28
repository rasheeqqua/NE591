#include "jacobi_preconditioner.h"
#include <cmath>

std::vector<std::vector<double>> generateJacobiPreconditioner(const std::vector<std::vector<double>>& A) {
    int n = A.size();

    // Create the Jacobi preconditioner matrix M
    // The Jacobi preconditioner is simply the diagonal elements of A
    // M = diag(A) = [ [A_11, 0, 0, ...],
    //                 [0, A_22, 0, ...],
    //                 [0, 0, A_33, ...],
    //                 [...           ] ]
    std::vector<std::vector<double>> M(n, std::vector<double>(n, 0.0));

    // Extract the diagonal elements of A and create M
    for (int i = 0; i < n; i++) {
        M[i][i] = A[i][i];
    }

    return M;
}

std::vector<std::vector<double>> generateJacobiPreconditionerInverse(const std::vector<std::vector<double>>& A) {
    int n = A.size();

    // Create the inverse of Jacobi preconditioner matrix M⁻¹
    // For Jacobi, M⁻¹ = diag(1/A_11, 1/A_22, 1/A_33, ...)
    std::vector<std::vector<double>> M_inv(n, std::vector<double>(n, 0.0));

    // Calculate the inverse of each diagonal element
    for (int i = 0; i < n; i++) {
        // Check if diagonal element is not too close to zero
        if (std::abs(A[i][i]) < 1e-10) {
            // Handle potential division by zero
            M_inv[i][i] = 1.0;  // Use identity matrix element instead
        } else {
            M_inv[i][i] = 1.0 / A[i][i];
        }
    }

    return M_inv;
}

// Apply preconditioner to a vector: M⁻¹ * r
std::vector<double> applyJacobiPreconditioner(const std::vector<std::vector<double>>& M_inv,
                                           const std::vector<double>& r) {
    int n = r.size();
    std::vector<double> result(n, 0.0);

    // For Jacobi preconditioner, M⁻¹ is diagonal
    // So M⁻¹ * r is simply element-wise multiplication of diagonal M⁻¹ with r
    // result_i = M⁻¹_ii * r_i
    for (int i = 0; i < n; i++) {
        result[i] = M_inv[i][i] * r[i];
    }

    return result;
}

// Apply preconditioner efficiently (without creating full matrix)
std::vector<double> applyJacobiPreconditionerEfficient(const std::vector<std::vector<double>>& A,
                                                    const std::vector<double>& r) {
    int n = r.size();
    std::vector<double> result(n, 0.0);

    // Directly divide each element of r by corresponding diagonal element of A
    // This is equivalent to M⁻¹ * r, but without creating the M⁻¹ matrix
    for (int i = 0; i < n; i++) {
        if (std::abs(A[i][i]) < 1e-10) {
            result[i] = r[i];  // Avoid division by zero
        } else {
            result[i] = r[i] / A[i][i];
        }
    }

    return result;
}