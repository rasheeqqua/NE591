#include "matrix_modules.h"

namespace matrix {

    // Vector operations
    double norm_inf(const std::vector<double>& v) {
        double max_val = 0.0;
        for (const auto& val : v) {
            max_val = std::max(max_val, std::abs(val));
        }
        return max_val;
    }

    double norm_2(const std::vector<double>& v) {
        double sum = 0.0;
        for (const auto& val : v) {
            sum += val * val;
        }
        return std::sqrt(sum);
    }

    std::vector<double> subtract(const std::vector<double>& a, const std::vector<double>& b) {
        std::vector<double> result(a.size());
        for (std::size_t i = 0; i < a.size(); ++i) {
            result[i] = a[i] - b[i];
        }
        return result;
    }

    // Matrix-vector operations
    std::vector<double> multiply(const std::vector<std::vector<double>>& A, const std::vector<double>& x) {
        std::vector<double> result(A.size(), 0.0);
        for (std::size_t i = 0; i < A.size(); ++i) {
            for (std::size_t j = 0; j < A[i].size(); ++j) {
                result[i] += A[i][j] * x[j];
            }
        }
        return result;
    }

    double residual_norm(const std::vector<std::vector<double>>& A, const std::vector<double>& x, const std::vector<double>& b) {
        std::vector<double> Ax = multiply(A, x);
        std::vector<double> residual = subtract(b, Ax);
        return norm_inf(residual);
    }

    // LUP decomposition and solving
    void lup_decompose(std::vector<std::vector<double>> A, std::vector<std::vector<double>>& L,
                      std::vector<std::vector<double>>& U, std::vector<int>& P) {
        int n = A.size();
        L.resize(n, std::vector<double>(n, 0.0));
        U = A;
        P.resize(n);

        // Initialize permutation vector
        for (int i = 0; i < n; ++i) {
            P[i] = i;
        }

        for (int k = 0; k < n; ++k) {
            // Find pivot
            double p_val = 0.0;
            int p_idx = k;

            for (int i = k; i < n; ++i) {
                if (std::abs(U[i][k]) > p_val) {
                    p_val = std::abs(U[i][k]);
                    p_idx = i;
                }
            }

            if (p_val == 0.0) {
                std::cerr << "Error: Matrix is singular" << std::endl;
                return;
            }

            // Swap rows in U and P
            if (p_idx != k) {
                std::swap(U[k], U[p_idx]);
                std::swap(P[k], P[p_idx]);

                // Swap rows in L up to column k-1
                for (int j = 0; j < k; ++j) {
                    std::swap(L[k][j], L[p_idx][j]);
                }
            }

            // Compute L and U
            L[k][k] = 1.0;

            for (int i = k + 1; i < n; ++i) {
                L[i][k] = U[i][k] / U[k][k];
                for (int j = k; j < n; ++j) {
                    U[i][j] -= L[i][k] * U[k][j];
                }
            }
        }
    }

    std::vector<double> lup_solve(const std::vector<std::vector<double>>& L,
                                 const std::vector<std::vector<double>>& U,
                                 const std::vector<int>& P,
                                 const std::vector<double>& b) {
        int n = L.size();
        std::vector<double> y(n, 0.0);
        std::vector<double> x(n, 0.0);

        // Forward substitution (Ly = Pb)
        for (int i = 0; i < n; ++i) {
            y[i] = b[P[i]];
            for (int j = 0; j < i; ++j) {
                y[i] -= L[i][j] * y[j];
            }
        }

        // Backward substitution (Ux = y)
        for (int i = n - 1; i >= 0; --i) {
            x[i] = y[i];
            for (int j = i + 1; j < n; ++j) {
                x[i] -= U[i][j] * x[j];
            }
            x[i] /= U[i][i];
        }

        return x;
    }

    // Iterative method functions
    void jacobi_iteration(const std::vector<std::vector<double>>& A,
                         const std::vector<double>& b,
                         std::vector<double>& x,
                         std::vector<double>& x_new) {
        std::size_t n = A.size();

        for (std::size_t i = 0; i < n; ++i) {
            double sum = b[i];
            for (std::size_t j = 0; j < n; ++j) {
                if (i != j) {
                    sum -= A[i][j] * x[j];
                }
            }
            x_new[i] = sum / A[i][i];
        }

        // Copy x_new to x
        x = x_new;
    }

    void gauss_seidel_iteration(const std::vector<std::vector<double>>& A,
                              const std::vector<double>& b,
                              std::vector<double>& x) {
        std::size_t n = A.size();

        for (std::size_t i = 0; i < n; ++i) {
            double sum = b[i];

            // Sum using updated values (j < i)
            for (std::size_t j = 0; j < i; ++j) {
                sum -= A[i][j] * x[j];
            }

            // Sum using previous iteration values (j > i)
            for (std::size_t j = i + 1; j < n; ++j) {
                sum -= A[i][j] * x[j];
            }

            x[i] = sum / A[i][i];
        }
    }

    void sor_iteration(const std::vector<std::vector<double>>& A,
                      const std::vector<double>& b,
                      std::vector<double>& x,
                      double omega) {
        std::size_t n = A.size();
        std::vector<double> x_old = x;

        for (std::size_t i = 0; i < n; ++i) {
            double sum = b[i];

            // Sum using updated values (j < i)
            for (std::size_t j = 0; j < i; ++j) {
                sum -= A[i][j] * x[j];
            }

            // Sum using previous iteration values (j > i)
            for (std::size_t j = i + 1; j < n; ++j) {
                sum -= A[i][j] * x_old[j];
            }

            x[i] = (1.0 - omega) * x_old[i] + omega * (sum / A[i][i]);
        }
    }

    // Matrix creation functions
    void extract_diag(const std::vector<std::vector<double>>& A, std::vector<double>& diag) {
        std::size_t n = A.size();
        diag.resize(n);

        for (std::size_t i = 0; i < n; ++i) {
            diag[i] = A[i][i];
        }
    }

    void extract_lower(const std::vector<std::vector<double>>& A, std::vector<std::vector<double>>& L) {
        std::size_t n = A.size();
        L.resize(n, std::vector<double>(n, 0.0));

        for (std::size_t i = 0; i < n; ++i) {
            for (std::size_t j = 0; j < i; ++j) {
                L[i][j] = A[i][j];
            }
        }
    }

    void extract_upper(const std::vector<std::vector<double>>& A, std::vector<std::vector<double>>& U) {
        std::size_t n = A.size();
        U.resize(n, std::vector<double>(n, 0.0));

        for (std::size_t i = 0; i < n; ++i) {
            for (std::size_t j = i + 1; j < n; ++j) {
                U[i][j] = A[i][j];
            }
        }
    }

    // Sparse matrix operations
    CSRMatrix convert_to_csr(const std::vector<std::vector<double>>& A) {
        CSRMatrix csr;
        csr.rows = A.size();
        csr.cols = A[0].size();
        csr.row_ptr.push_back(0);

        for (std::size_t i = 0; i < A.size(); ++i) {
            for (std::size_t j = 0; j < A[i].size(); ++j) {
                if (std::abs(A[i][j]) > 1e-10) {
                    csr.values.push_back(A[i][j]);
                    csr.col_indices.push_back(j);
                }
            }
            csr.row_ptr.push_back(csr.values.size());
        }

        return csr;
    }

    std::vector<double> multiply_csr(const CSRMatrix& A, const std::vector<double>& x) {
        std::vector<double> result(A.rows, 0.0);

        for (int i = 0; i < A.rows; ++i) {
            for (int j = A.row_ptr[i]; j < A.row_ptr[i + 1]; ++j) {
                result[i] += A.values[j] * x[A.col_indices[j]];
            }
        }

        return result;
    }

    void jacobi_iteration_csr(const CSRMatrix& A,
                             const std::vector<double>& b,
                             std::vector<double>& x,
                             std::vector<double>& x_new) {
        // Extract diagonal elements
        std::vector<double> diag(A.rows);
        for (int i = 0; i < A.rows; ++i) {
            for (int j = A.row_ptr[i]; j < A.row_ptr[i + 1]; ++j) {
                if (A.col_indices[j] == i) {
                    diag[i] = A.values[j];
                    break;
                }
            }
        }

        // Perform Jacobi iteration
        for (int i = 0; i < A.rows; ++i) {
            double sum = b[i];
            for (int j = A.row_ptr[i]; j < A.row_ptr[i + 1]; ++j) {
                if (A.col_indices[j] != i) {
                    sum -= A.values[j] * x[A.col_indices[j]];
                }
            }
            x_new[i] = sum / diag[i];
        }

        // Copy x_new to x
        x = x_new;
    }

    void gauss_seidel_iteration_csr(const CSRMatrix& A,
                                  const std::vector<double>& b,
                                  std::vector<double>& x) {
        // Extract diagonal elements
        std::vector<double> diag(A.rows);
        for (int i = 0; i < A.rows; ++i) {
            for (int j = A.row_ptr[i]; j < A.row_ptr[i + 1]; ++j) {
                if (A.col_indices[j] == i) {
                    diag[i] = A.values[j];
                    break;
                }
            }
        }

        // Perform Gauss-Seidel iteration
        for (int i = 0; i < A.rows; ++i) {
            double sum = b[i];
            for (int j = A.row_ptr[i]; j < A.row_ptr[i + 1]; ++j) {
                if (A.col_indices[j] != i) {
                    sum -= A.values[j] * x[A.col_indices[j]];
                }
            }
            x[i] = sum / diag[i];
        }
    }

    // Print functions for debugging
    void print_matrix(const std::vector<std::vector<double>>& A) {
        for (const auto& row : A) {
            for (const auto& val : row) {
                std::cout << val << " ";
            }
            std::cout << std::endl;
        }
    }

    void print_vector(const std::vector<double>& v) {
        for (const auto& val : v) {
            std::cout << val << " ";
        }
        std::cout << std::endl;
    }
}