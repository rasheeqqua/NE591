#include <fstream>
#include <vector>
#include <iomanip>
#include <chrono>
#include <cmath>
#include "PointJacobi.cpp"
#include "GaussSeidel.cpp"
#include "SOR.cpp"
#include "LUPFactorization/LUPFactorization.cpp"
#include "LUPFactorization/ApplyPermutationMatrix.cpp"
#include "LUPFactorization/Substitution.cpp"
#include "LUPFactorization/CalculateResiduals.cpp"

// Function to calculate the largest relative difference between two solution vectors
// Formula: δ^μ = max(1≤i≤n) |x_i^μ/x_i^LUP - 1|
double calculateMaxRelativeDifference(const std::vector<double>& x_reference,
                                    const std::vector<double>& x_iterative) {
    double max_diff = 0.0;
    for (size_t i = 0; i < x_reference.size(); i++) {
        // Avoid division by zero
        if (std::abs(x_reference[i]) > 1e-15) {
            double relative_diff = std::abs(x_iterative[i] / x_reference[i] - 1.0);
            max_diff = std::max(max_diff, relative_diff);
        }
    }
    return max_diff;
}

// Function to generate matrix A according to the specified requirements
std::vector<std::vector<double>> generateMatrixA(int n) {
    std::vector<std::vector<double>> A(n, std::vector<double>(n));

    for (int i = 0; i < n; i++) {
        double rowSum = 0.0;
        for (int j = 0; j < n; j++) {
            if (i != j) {
                A[i][j] = -1.0 / ((i + 1) + (j + 1));
                rowSum += std::abs(A[i][j]);
            }
        }
        A[i][i] = 1.0/n + rowSum;
    }
    return A;
}

// Function to generate vector b
std::vector<double> generateVectorB(int n) {
    return std::vector<double>(n, 1.0);
}

// Function to solve and print results for a given matrix size
void solveAndPrintResults(int n, std::ofstream& output, double epsilon = 1e-4, double omega = 1.5) {
    output << "\nResults for n = " << n << "\n";
    output << "----------------------------------------\n";

    // Generate matrix A and vector b
    auto startGen = std::chrono::high_resolution_clock::now();
    std::vector<std::vector<double>> A = generateMatrixA(n);
    std::vector<double> b = generateVectorB(n);
    auto endGen = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> genTime = endGen - startGen;

    // Print matrix A and vector b for small matrices
    if (n <= 10) {
        output << "\nMatrix A:\n";
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                output << std::scientific << std::setprecision(4) << A[i][j] << " ";
            }
            output << "\n";
        }

        output << "\nVector b:\n";
        for (int i = 0; i < n; i++) {
            output << std::scientific << std::setprecision(4) << b[i] << " ";
        }
        output << "\n";
    }

    // Solve using LUP factorization (reference solution)
    std::vector<std::vector<double>> L(n, std::vector<double>(n));
    std::vector<std::vector<double>> U(n, std::vector<double>(n));
    std::vector<std::vector<double>> P(n, std::vector<double>(n));
    std::vector<double> y(n), x_lup(n), Pb(n);

    auto startLUP = std::chrono::high_resolution_clock::now();
    bool lupSuccess = lupFactorize(A, L, U, P);
    if (lupSuccess) {
        applyPermutationMatrix(P, b, Pb);
        forwardSubstitution(L, Pb, y);
        backSubstitution(U, y, x_lup);
    }
    auto endLUP = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> lupTime = endLUP - startLUP;

    // Calculate LUP residual
    double lupResidual = calculateResidual(A, x_lup, b);

    // Print LUP results
    output << "\nLUP Factorization Results:\n";
    output << "Execution time: " << std::scientific << std::setprecision(4) << lupTime.count() << " seconds\n";
    output << "Maximum residual: " << std::scientific << std::setprecision(4) << lupResidual << "\n";

    if (n <= 10) {
        output << "Solution vector x_lup:\n";
        for (int i = 0; i < n; i++) {
            output << std::scientific << std::setprecision(8) << x_lup[i] << "\n";
        }
    }

    // Solve using relaxation methods
    const int maxIter = 10000;

    // Point Jacobi
    std::vector<double> x_jacobi(n, 0.0);
    int jacobi_iter;
    double jacobi_error;
    auto startJacobi = std::chrono::high_resolution_clock::now();
    bool jacobi_converged = pointJacobi(A, b, x_jacobi, epsilon, maxIter, jacobi_iter, jacobi_error);
    auto endJacobi = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> jacobiTime = endJacobi - startJacobi;

    // Gauss-Seidel
    std::vector<double> x_gs(n, 0.0);
    int gs_iter;
    double gs_error;
    auto startGS = std::chrono::high_resolution_clock::now();
    bool gs_converged = gaussSeidel(A, b, x_gs, epsilon, maxIter, gs_iter, gs_error);
    auto endGS = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> gsTime = endGS - startGS;

    // SOR
    std::vector<double> x_sor(n, 0.0);
    int sor_iter;
    double sor_error;
    auto startSOR = std::chrono::high_resolution_clock::now();
    bool sor_converged = sor(A, b, x_sor, omega, epsilon, maxIter, sor_iter, sor_error);
    auto endSOR = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> sorTime = endSOR - startSOR;

    // Calculate relative differences with respect to LUP solution
    output << "\nRelative Differences from LUP Solution:\n";
    output << "----------------------------------------\n";

    if (jacobi_converged) {
        double jacobi_rel_diff = calculateMaxRelativeDifference(x_lup, x_jacobi);
        output << "Point Jacobi max relative difference: " << std::scientific
               << std::setprecision(4) << jacobi_rel_diff << "\n";
    }

    if (gs_converged) {
        double gs_rel_diff = calculateMaxRelativeDifference(x_lup, x_gs);
        output << "Gauss-Seidel max relative difference: " << std::scientific
               << std::setprecision(4) << gs_rel_diff << "\n";
    }

    if (sor_converged) {
        double sor_rel_diff = calculateMaxRelativeDifference(x_lup, x_sor);
        output << "SOR max relative difference: " << std::scientific
               << std::setprecision(4) << sor_rel_diff << "\n";
    }

    // Print relaxation method results
    output << "\nRelaxation Methods Results (ε = " << epsilon << "):\n";
    output << "\nPoint Jacobi:\n";
    output << "Converged: " << (jacobi_converged ? "Yes" : "No") << "\n";
    if (jacobi_converged) {
        output << "Iterations: " << jacobi_iter << "\n";
        output << "Final error: " << std::scientific << std::setprecision(4) << jacobi_error << "\n";
        output << "Execution time: " << std::scientific << std::setprecision(4) << jacobiTime.count() << " seconds\n";
        output << "Maximum residual: " << std::scientific << std::setprecision(4)
               << calculateResidual(A, x_jacobi, b) << "\n";
    }

    output << "\nGauss-Seidel:\n";
    output << "Converged: " << (gs_converged ? "Yes" : "No") << "\n";
    if (gs_converged) {
        output << "Iterations: " << gs_iter << "\n";
        output << "Final error: " << std::scientific << std::setprecision(4) << gs_error << "\n";
        output << "Execution time: " << std::scientific << std::setprecision(4) << gsTime.count() << " seconds\n";
        output << "Maximum residual: " << std::scientific << std::setprecision(4)
               << calculateResidual(A, x_gs, b) << "\n";
    }

    output << "\nSOR (ω = " << omega << "):\n";
    output << "Converged: " << (sor_converged ? "Yes" : "No") << "\n";
    if (sor_converged) {
        output << "Iterations: " << sor_iter << "\n";
        output << "Final error: " << std::scientific << std::setprecision(4) << sor_error << "\n";
        output << "Execution time: " << std::scientific << std::setprecision(4) << sorTime.count() << " seconds\n";
        output << "Maximum residual: " << std::scientific << std::setprecision(4)
               << calculateResidual(A, x_sor, b) << "\n";
    }

    output << "\n";
}

int main() {
    std::ofstream output("matrix_results.txt");

    // Write header
    output << "NE 591 Matrix Generation and Solution Program\n";
    output << "Implemented by Hasibul H. Rasheeq, 02/20/25\n";
    output << "-------------------------------------------\n\n";

    // Define parameters
    const double epsilon = 1e-4;  // Stopping criterion for relaxation methods
    const double omega = 1.5;     // SOR relaxation parameter

    // Vector of matrix sizes to test
    std::vector<int> sizes = {32, 64, 128, 512, 1024};

    // Solve for each matrix size
    for (int n : sizes) {
        solveAndPrintResults(n, output, epsilon, omega);
    }

    output.close();
    return 0;
}