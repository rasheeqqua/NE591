//
// Created by Hasibul H. Rasheeq on 02/20/25.
//
// Point-Jacobi method

#include "CalculateResidual.cpp"

bool pointJacobi(const std::vector<std::vector<double>>& A,
                 const std::vector<double>& b,
                 std::vector<double>& x,
                 double epsilon,
                 int maxIter,
                 int& actualIter,
                 double& finalError) {
    int n = A.size();
    std::vector<double> x_new(n, 0.0);

    for (actualIter = 0; actualIter < maxIter; actualIter++) {
        for (int i = 0; i < n; i++) {
            double sum = 0.0;
            for (int j = 0; j < n; j++) {
                if (i != j) sum += A[i][j] * x[j];
            }
            x_new[i] = (b[i] - sum) / A[i][i];
        }

        finalError = calculateError(x_new, x);
        if (finalError < epsilon) {
            x = x_new;
            return true;
        }
        x = x_new;
    }
    return false;
}
