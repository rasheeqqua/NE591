//
// Created by Hasibul H. Rasheeq on 02/13/25.
//

#include <vector>

// Matrix-vector multiplication
std::vector<double> matrixVectorProduct(const std::vector<std::vector<double>>& A,
                                      const std::vector<double>& x) {
    int n = A.size();
    std::vector<double> result(n, 0.0);

    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            result[i] += A[i][j] * x[j];
        }
    }
    return result;
}
