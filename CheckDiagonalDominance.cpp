//
// Created by Hasibul H. Rasheeq on 02/14/25.
//

#include <cmath>
#include <vector>

// Function to check if matrix is diagonally dominant
bool isDiagonallyDominant(const std::vector<std::vector<double>>& A) {
    int n = A.size();
    for (int i = 0; i < n; i++) {
        double diagElement = std::abs(A[i][i]);
        double sum = 0.0;
        for (int j = 0; j < n; j++) {
            if (i != j) sum += std::abs(A[i][j]);
        }
        if (diagElement <= sum) return false;
    }
    return true;
}
