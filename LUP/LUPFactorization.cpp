#include <vector>
#include <cmath>
#include <algorithm>

bool verifyLUPFactorization(const std::vector<std::vector<double>>& A,
                          const std::vector<std::vector<double>>& L,
                          const std::vector<std::vector<double>>& U,
                          const std::vector<std::vector<double>>& P) {
    int n = A.size();
    double eps = 1e-10;

    // Verify L*U = P*A
    std::vector<std::vector<double>> PA(n, std::vector<double>(n, 0.0));
    std::vector<std::vector<double>> LU(n, std::vector<double>(n, 0.0));

    // Compute P*A
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            for (int k = 0; k < n; ++k) {
                PA[i][j] += P[i][k] * A[k][j];
            }
        }
    }

    // Compute L*U
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            for (int k = 0; k < n; ++k) {
                LU[i][j] += L[i][k] * U[k][j];
            }
        }
    }

    // Compare PA and LU
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            if (std::abs(PA[i][j] - LU[i][j]) > eps) {
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

bool lupFactorize(const std::vector<std::vector<double>>& A,
                 std::vector<std::vector<double>>& L,
                 std::vector<std::vector<double>>& U,
                 std::vector<std::vector<double>>& P) {
    int n = A.size();

    // Create U as a copy of A
    U = A;

    // Initialize L's diagonal to 1 and P as identity matrix
    L = std::vector<std::vector<double>>(n, std::vector<double>(n, 0.0));
    P = std::vector<std::vector<double>>(n, std::vector<double>(n, 0.0));
    for (int i = 0; i < n; ++i) {
        L[i][i] = 1.0;
        P[i][i] = 1.0;
    }

    for (int k = 0; k < n - 1; ++k) {
        int maxRow = k;
        double maxVal = std::abs(U[k][k]);

        // Find the maximum element in current column
        for (int i = k + 1; i < n; ++i) {
            if (std::abs(U[i][k]) > maxVal) {
                maxVal = std::abs(U[i][k]);
                maxRow = i;
            }
        }

        // Swap rows if necessary
        if (maxRow != k) {
            // Swap rows in U
            std::swap(U[k], U[maxRow]);

            // Swap rows in L (only elements before k)
            for (int j = 0; j < k; ++j) {
                std::swap(L[k][j], L[maxRow][j]);
            }

            // Update P matrix
            std::swap(P[k], P[maxRow]);
        }

        // Check for zero pivot
        if (std::abs(U[k][k]) < 1e-10) {
            return false;
        }

        for (int i = k + 1; i < n; ++i) {
            L[i][k] = U[i][k] / U[k][k];
            for (int j = k; j < n; ++j) {
                U[i][j] -= L[i][k] * U[k][j];
            }
        }
    }

    // Verify the factorization
    return verifyLUPFactorization(A, L, U, P);
}
