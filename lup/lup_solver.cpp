#include "lup_solver.h"
#include <cmath>
#include <iostream>
#include <algorithm>

// LUP decomposition of matrix A into L, U, and P
void lupDecomposition(const std::vector<std::vector<double>>& A,
                      std::vector<std::vector<double>>& L,
                      std::vector<std::vector<double>>& U,
                      std::vector<int>& P) {
    int n = A.size();

    // Initialize L, U, and P
    L.resize(n, std::vector<double>(n, 0.0));
    U = A;  // Copy A to U
    P.resize(n);

    // Initialize permutation vector
    for (int i = 0; i < n; ++i) {
        P[i] = i;
    }

    // Perform LUP decomposition
    for (int k = 0; k < n; ++k) {
        // Find pivot
        double p_val = 0.0;
        int p_idx = k;

        for (int i = k; i < n; ++i) {
            if (std::abs(U[i][k]) > p_val) {
                p_val = std::abs(U[i][k]);
                p_idx = i;
            }
        }

        if (p_val < 1e-10) {
            std::cerr << "Error: Matrix is singular or nearly singular" << std::endl;
            // Continue with a small value to avoid division by zero
            U[k][k] = 1e-10;
        } else {
            // Swap rows in U and P
            if (p_idx != k) {
                std::swap(U[k], U[p_idx]);
                std::swap(P[k], P[p_idx]);

                // Swap rows in L up to column k-1
                for (int j = 0; j < k; ++j) {
                    std::swap(L[k][j], L[p_idx][j]);
                }
            }
        }

        // Compute L and U
        L[k][k] = 1.0;  // Diagonal of L is 1

        for (int i = k + 1; i < n; ++i) {
            L[i][k] = U[i][k] / U[k][k];

            for (int j = k; j < n; ++j) {
                U[i][j] -= L[i][k] * U[k][j];
            }
        }
    }
}

// Solve Ax = b using LUP decomposition (L, U, P)
std::vector<double> lupSolveSystem(const std::vector<std::vector<double>>& L,
                                  const std::vector<std::vector<double>>& U,
                                  const std::vector<int>& P,
                                  const std::vector<double>& b) {
    int n = L.size();
    std::vector<double> y(n, 0.0);  // Temporary solution for Ly = Pb
    std::vector<double> x(n, 0.0);  // Final solution for Ux = y

    // Forward substitution to solve Ly = Pb
    for (int i = 0; i < n; ++i) {
        y[i] = b[P[i]];  // Apply permutation to b

        for (int j = 0; j < i; ++j) {
            y[i] -= L[i][j] * y[j];
        }
        // No division needed since L[i][i] = 1.0
    }

    // Backward substitution to solve Ux = y
    for (int i = n - 1; i >= 0; --i) {
        x[i] = y[i];

        for (int j = i + 1; j < n; ++j) {
            x[i] -= U[i][j] * x[j];
        }

        x[i] /= U[i][i];
    }

    return x;
}

// Combined function to solve Ax = b using LUP decomposition
void lupSolve(const std::vector<std::vector<double>>& A,
             const std::vector<double>& b,
             std::vector<double>& x) {
    int n = A.size();

    // Check dimensions
    if (A[0].size() != n || b.size() != n) {
        std::cerr << "Error: Matrix dimensions mismatch" << std::endl;
        return;
    }

    // Perform LUP decomposition
    std::vector<std::vector<double>> L, U;
    std::vector<int> P;
    lupDecomposition(A, L, U, P);

    // Solve the system
    x = lupSolveSystem(L, U, P, b);
}