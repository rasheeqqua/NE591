//
// Created by Hasibul H. Rasheeq on 02/08/25.
//

#include <vector>
#include <cmath>

bool verifyFactorization(const std::vector<std::vector<double>>& A,
                        const std::vector<std::vector<double>>& L,
                        const std::vector<std::vector<double>>& U) {
    int n = A.size();
    double eps = 1e-10;  // Tolerance for floating-point comparison

    // Verify L*U = A
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            double sum = 0.0;
            for (int k = 0; k < n; ++k) {
                sum += L[i][k] * U[k][j];
            }
            if (std::abs(sum - A[i][j]) > eps) {
                return false;
            }
        }
    }

    // Verify L is lower triangular with ones on diagonal
    for (int i = 0; i < n; ++i) {
        if (std::abs(L[i][i] - 1.0) > eps) return false;
        for (int j = i + 1; j < n; ++j) {
            if (std::abs(L[i][j]) > eps) return false;
        }
    }

    // Verify U is upper triangular
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < i; ++j) {
            if (std::abs(U[i][j]) > eps) return false;
        }
    }

    return true;
}

bool luFactorize(const std::vector<std::vector<double>>& A,
                 std::vector<std::vector<double>>& L,
                 std::vector<std::vector<double>>& U) {
    int n = A.size();
    U = A;  // Copy A to U for in-place factorization
    L = std::vector<std::vector<double>>(n, std::vector<double>(n, 0.0));

    // Initialize L's diagonal to 1
    for (int i = 0; i < n; ++i) {
        L[i][i] = 1.0;
    }

    // Perform LU factorization
    for (int k = 0; k < n - 1; ++k) {
        if (std::abs(U[k][k]) < 1e-10) return false;  // Check for zero pivot

        for (int i = k + 1; i < n; ++i) {
            L[i][k] = U[i][k] / U[k][k];
            for (int j = k; j < n; ++j) {
                U[i][j] -= L[i][k] * U[k][j];
            }
        }
    }

    // Verify factorization
    if (!verifyFactorization(A, L, U)) {
        return false;
    }

    return true;
}
