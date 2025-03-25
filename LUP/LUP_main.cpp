#include <vector>
#include <cmath>
#include <fstream>
#include <iomanip>
#include "LUP_solver.h"

// Function to compute residual vector and its norm
std::vector<double> computeResidual(const std::vector<std::vector<double>>& A,
                                   const std::vector<double>& x,
                                   const std::vector<double>& b) {
    int n = A.size();
    std::vector<double> r(n, 0.0);

    // r = b - A*x
    for (int i = 0; i < n; ++i) {
        r[i] = b[i];
        for (int j = 0; j < n; ++j) {
            r[i] -= A[i][j] * x[j];
        }
    }

    return r;
}

// Function to compute maximum absolute value in a vector
double maxAbsoluteValue(const std::vector<double>& v) {
    double maxVal = 0.0;
    for (double value : v) {
        maxVal = std::max(maxVal, std::abs(value));
    }
    return maxVal;
}

// Main LUP solver function
bool solveLUP(const std::vector<std::vector<double>>& A,
             const std::vector<double>& b,
             std::vector<double>& x,
             double& maxResidual) {
    int n = A.size();

    // Initialize L, U, P matrices
    std::vector<std::vector<double>> L(n, std::vector<double>(n, 0.0));
    std::vector<std::vector<double>> U(n, std::vector<double>(n, 0.0));
    std::vector<std::vector<double>> P(n, std::vector<double>(n, 0.0));

    // Perform LUP factorization
    if (!lupFactorize(A, L, U, P)) {
        return false; // Factorization failed
    }

    // Apply permutation to b
    std::vector<double> Pb(n, 0.0);
    applyPermutationMatrix(P, b, Pb);

    // Solve Ly = Pb
    std::vector<double> y(n, 0.0);
    forwardSubstitution(L, Pb, y);

    // Solve Ux = y
    x.resize(n, 0.0);
    backSubstitution(U, y, x);

    // Compute residual
    std::vector<double> residual = computeResidual(A, x, b);
    maxResidual = maxAbsoluteValue(residual);

    return true;
}

// Function to write LUP results to output file
void writeLUPResults(std::ofstream& outFile,
                    const std::vector<double>& x,
                    double maxResidual) {
    // Write specific LUP results
    outFile << "Direct solution using LUP decomposition" << std::endl << std::endl;

    outFile << "Maximum absolute residual: " << std::scientific << std::setprecision(4) << maxResidual << std::endl << std::endl;

    outFile << "Solution vector, x:" << std::endl;
    for (std::size_t i = 0; i < x.size(); ++i) {
        outFile << std::scientific << std::setprecision(4) << x[i] << " ";
    }
    outFile << std::endl << std::endl;
}