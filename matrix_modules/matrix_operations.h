#ifndef MATRIX_OPERATIONS_H
#define MATRIX_OPERATIONS_H

#include <vector>
#include <cmath>

// a. Multiplication of a vector by a scalar: c.y
std::vector<double> scalar_multiply(const std::vector<double>& y, double c);

// b. Sum of two vectors: y+z
std::vector<double> vector_add(const std::vector<double>& y, const std::vector<double>& z);

// c. Scalar product of two vectors: (y, z) = y^T.z
double scalar_product(const std::vector<double>& y, const std::vector<double>& z);

// d. A-inner product: (y, z)_A = y^T.A.z
double a_inner_product(const std::vector<double>& y, const std::vector<double>& z,
                      const std::vector<std::vector<double>>& A);

// Matrix-vector product (from Outlab 1)
std::vector<double> matrix_vector_product(const std::vector<std::vector<double>>& A,
                                         const std::vector<double>& x);

// Vector subtraction: y-z
std::vector<double> vector_subtract(const std::vector<double>& y, const std::vector<double>& z);

// Calculate L2 norm of a vector
double calculate_L2_norm(const std::vector<double>& v);

#endif // MATRIX_OPERATIONS_H