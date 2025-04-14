// lu_factorization.cpp
// Implementation of LU factorization solver for linear systems
// April 13, 2025

#include "lu_factorization.h"
#include <cmath>
#include <algorithm>

bool luFactorize(double** A, int n, int* perm) {
    // Initialize permutation array
    for (int i = 0; i < n; i++) {
        perm[i] = i;
    }

    // LU factorization with partial pivoting
    for (int k = 0; k < n; k++) {
        // Find pivot
        double maxVal = 0.0;
        int pivotRow = k;

        for (int i = k; i < n; i++) {
            double absVal = std::abs(A[perm[i]][k]);
            if (absVal > maxVal) {
                maxVal = absVal;
                pivotRow = i;
            }
        }

        // Check for singularity
        if (maxVal < 1e-10) {
            return false; // Matrix is singular
        }

        // Swap rows in permutation
        if (pivotRow != k) {
            std::swap(perm[k], perm[pivotRow]);
        }

        // Compute multipliers and eliminate
        for (int i = k + 1; i < n; i++) {
            int permI = perm[i];
            int permK = perm[k];

            A[permI][k] /= A[permK][k];

            for (int j = k + 1; j < n; j++) {
                A[permI][j] -= A[permI][k] * A[permK][j];
            }
        }
    }

    return true;
}

void luBacksubstitute(double** A, double* b, double* x, int n, int* perm) {
    // Allocate temporary array for permuted right-hand side
    double* y = new double[n];

    // Forward substitution (Ly = Pb)
    for (int i = 0; i < n; i++) {
        y[i] = b[perm[i]];
        for (int j = 0; j < i; j++) {
            y[i] -= A[perm[i]][j] * y[j];
        }
    }

    // Backward substitution (Ux = y)
    for (int i = n - 1; i >= 0; i--) {
        x[i] = y[i];
        for (int j = i + 1; j < n; j++) {
            x[i] -= A[perm[i]][j] * x[j];
        }
        x[i] /= A[perm[i]][i];
    }

    // Free memory
    delete[] y;
}

bool luSolve(double** A, double* b, double* x, int n) {
    // Allocate permutation array
    int* perm = new int[n];

    // Perform LU factorization
    bool success = luFactorize(A, n, perm);

    if (success) {
        // Solve the system using the factorization
        luBacksubstitute(A, b, x, n, perm);
    }

    // Free memory
    delete[] perm;

    return success;
}