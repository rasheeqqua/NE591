#ifndef MATRIX_MODULES_H
#define MATRIX_MODULES_H

#include <vector>
#include <cmath>
#include <iostream>
#include <algorithm>

namespace matrix {
    // Vector operations
    double norm_inf(const std::vector<double>& v);
    double norm_2(const std::vector<double>& v);
    std::vector<double> subtract(const std::vector<double>& a, const std::vector<double>& b);

    // Matrix-vector operations
    std::vector<double> multiply(const std::vector<std::vector<double>>& A, const std::vector<double>& x);
    double residual_norm(const std::vector<std::vector<double>>& A, const std::vector<double>& x, const std::vector<double>& b);

    // LUP decomposition and solving
    void lup_decompose(std::vector<std::vector<double>> A, std::vector<std::vector<double>>& L,
                      std::vector<std::vector<double>>& U, std::vector<int>& P);
    std::vector<double> lup_solve(const std::vector<std::vector<double>>& L,
                                 const std::vector<std::vector<double>>& U,
                                 const std::vector<int>& P,
                                 const std::vector<double>& b);

    // Iterative method functions
    void jacobi_iteration(const std::vector<std::vector<double>>& A,
                         const std::vector<double>& b,
                         std::vector<double>& x,
                         std::vector<double>& x_new);

    void gauss_seidel_iteration(const std::vector<std::vector<double>>& A,
                              const std::vector<double>& b,
                              std::vector<double>& x);

    void sor_iteration(const std::vector<std::vector<double>>& A,
                      const std::vector<double>& b,
                      std::vector<double>& x,
                      double omega);

    // Matrix creation functions
    void extract_diag(const std::vector<std::vector<double>>& A, std::vector<double>& diag);
    void extract_lower(const std::vector<std::vector<double>>& A, std::vector<std::vector<double>>& L);
    void extract_upper(const std::vector<std::vector<double>>& A, std::vector<std::vector<double>>& U);

    // Sparse matrix operations (using CSR format)
    struct CSRMatrix {
        std::vector<double> values;
        std::vector<int> col_indices;
        std::vector<int> row_ptr;
        int rows;
        int cols;
    };

    CSRMatrix convert_to_csr(const std::vector<std::vector<double>>& A);
    std::vector<double> multiply_csr(const CSRMatrix& A, const std::vector<double>& x);
    void jacobi_iteration_csr(const CSRMatrix& A,
                             const std::vector<double>& b,
                             std::vector<double>& x,
                             std::vector<double>& x_new);
    void gauss_seidel_iteration_csr(const CSRMatrix& A,
                                  const std::vector<double>& b,
                                  std::vector<double>& x);

    // Print functions for debugging
    void print_matrix(const std::vector<std::vector<double>>& A);
    void print_vector(const std::vector<double>& v);
}

#endif // MATRIX_MODULES_H