#include "matrix_operations.h"

// a. Multiplication of a vector by a scalar: c.y
std::vector<double> scalar_multiply(const std::vector<double>& y, double c) {
    std::vector<double> result(y.size());
    for (std::size_t i = 0; i < y.size(); i++) {
        result[i] = c * y[i];
    }
    return result;
}

// b. Sum of two vectors: y+z
std::vector<double> vector_add(const std::vector<double>& y, const std::vector<double>& z) {
    std::vector<double> result(y.size());
    for (std::size_t i = 0; i < y.size(); i++) {
        result[i] = y[i] + z[i];
    }
    return result;
}

// c. Scalar product of two vectors: (y, z) = y^T.z
double scalar_product(const std::vector<double>& y, const std::vector<double>& z) {
    // Note: For vectors, the transpose operation doesn't change the computation
    // y^T.z for vectors is the dot product which is sum(y[i]*z[i])
    double result = 0.0;
    for (std::size_t i = 0; i < y.size(); i++) {
        result += y[i] * z[i];
    }
    return result;
}

// d. A-inner product: (y, z)_A = y^T.A.z
double a_inner_product(const std::vector<double>& y, const std::vector<double>& z,
                      const std::vector<std::vector<double>>& A) {
    // First compute A.z
    std::vector<double> Az = matrix_vector_product(A, z);

    // Then compute y^T.(A.z)
    // For vectors, this is just the dot product
    return scalar_product(y, Az);
}

// Matrix-vector product (from Outlab 1)
std::vector<double> matrix_vector_product(const std::vector<std::vector<double>>& A,
                                         const std::vector<double>& x) {
    std::size_t n = A.size();
    std::vector<double> result(n, 0.0);

    for (std::size_t i = 0; i < n; i++) {
        for (std::size_t j = 0; j < n; j++) {
            result[i] += A[i][j] * x[j];
        }
    }

    return result;
}