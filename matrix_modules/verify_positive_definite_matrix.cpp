#include "verify_positive_definite_matrix.h"
#include <cmath>

bool isSymmetric(const std::vector<std::vector<double>>& A) {
    int n = A.size();
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            // Compare matrix elements with tolerance to account for floating-point errors
            if (std::abs(A[i][j] - A[j][i]) > 1e-10) {
                return false;
            }
        }
    }
    return true;
}

bool isDiagonallyDominant(const std::vector<std::vector<double>>& A) {
    int n = A.size();
    for (int i = 0; i < n; i++) {
        double diagonal = std::abs(A[i][i]);
        double sum = 0.0;

        // Sum of absolute values of all off-diagonal elements in row i
        for (int j = 0; j < n; j++) {
            if (i != j) {
                sum += std::abs(A[i][j]);
            }
        }

        // Check if diagonal element is less than sum of off-diagonal elements
        if (diagonal < sum) {
            return false;
        }
    }
    return true;
}