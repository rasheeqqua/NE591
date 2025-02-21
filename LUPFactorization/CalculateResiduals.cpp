//
// Created by Hasibul H. Rasheeq on 02/08/25.
//
#include <cmath>
#include "MatrixVectorProduct.cpp"

std::vector<double> calculateResiduals(const std::vector<std::vector<double>>& A,
                                     const std::vector<double>& x,
                                     const std::vector<double>& b) {
    std::vector<double> Ax = matrixVectorProduct(A, x);
    std::vector<double> residuals(b.size());
    for (int i = 0; i < b.size(); ++i) {
        residuals[i] = std::abs(Ax[i] - b[i]);
    }
    return residuals;
}
