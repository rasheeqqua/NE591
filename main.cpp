#include <fstream>
#include <vector>
#include <iomanip>
#include <chrono>
#include "PointJacobi.cpp"
#include "GaussSeidel.cpp"
#include "SOR.cpp"
#include "LUPFactorization/LUPFactorization.cpp"
#include "LUPFactorization/ApplyPermutationMatrix.cpp"
#include "LUPFactorization/Substitution.cpp"
#include "LUPFactorization/CalculateResiduals.cpp"

// Function to validate input parameters
bool validateInputs(int n, int maxIter, double epsilon, std::ofstream& output) {
    bool isValid = true;

    if (n <= 0) {
        output << "Error: Matrix order must be positive.\n";
        isValid = false;
    }

    if (maxIter <= 0) {
        output << "Error: Maximum iterations must be positive.\n";
        isValid = false;
    }

    if (epsilon <= 0) {
        output << "Error: Stopping criterion must be positive.\n";
        isValid = false;
    }

    return isValid;
}

int main() {
    // Solve system based on method selection
    auto IOStartTime = std::chrono::high_resolution_clock::now();

    std::ifstream input("../input.txt");
    std::ofstream output("output.txt");

    // Write header
    output << "NE 591 Outlab 06 Code\n";
    output << "Implemented by Hasibul H. Rasheeq, 02/20/25\n";
    output << "-------------------------------------------\n\n";
    output << "Solve Matrix Equation with Relaxation Iterations\n";
    output << "------------------------------------------------\n\n";

    // Read input data
    int method;
    double omega;
    input >> method >> omega;

    double epsilon;
    int maxIter;
    input >> epsilon >> maxIter;

    int n;
    input >> n;

    // Validate input parameters before proceeding
    if (!validateInputs(n, maxIter, epsilon, output)) {
        input.close();
        output.close();
        return 1;
    }

    std::vector<std::vector<double>> A(n, std::vector<double>(n));
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            input >> A[i][j];
        }
    }

    std::vector<double> b(n);
    for (int i = 0; i < n; i++) {
        input >> b[i];
    }

    // Print method selection
    switch (method) {
        case 0: output << "selected Point-Jacobi Iterations\n"; break;
        case 1: output << "selected Gauss-Seidel Iterations\n"; break;
        case 2: output << "selected SOR Iterations\n"; break;
    }

    output << "Stopping criterion = " << std::scientific << std::setprecision(2) << epsilon << "\n";
    output << "Iteration limit = " << maxIter << "\n";
    output << "Matrix is of order = " << n << "\n";

    // Print matrix A
    output << "Matrix A:\n";
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            output << std::scientific << std::setprecision(1) << A[i][j] << " ";
        }
        output << "\n";
    }

    // Print vector b
    output << "\nRHS Vector b:\n";
    for (int i = 0; i < n; i++) {
        output << std::scientific << std::setprecision(1) << b[i] << " ";
    }
    output << "\n\n";

    std::vector<double> x(n, 0.0);
    int actualIter;
    double finalError;
    bool converged;

    std::vector<std::vector<double>> L(n, std::vector<double>(n));
    std::vector<std::vector<double>> U(n, std::vector<double>(n));
    std::vector<std::vector<double>> P(n, std::vector<double>(n));
    std::vector<double> y(n), x_lup(n), Pb(n);

    // Print execution time
    auto IOEndTime = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> IODuration = IOEndTime - IOStartTime;

    auto RelaxationStartTime = std::chrono::high_resolution_clock::now();

    switch (method) {
        case 0:
            converged = pointJacobi(A, b, x, epsilon, maxIter, actualIter, finalError);
            if (converged) {
                output << "Iterations converged in: " << actualIter << " iterations\n";
                output << "Iterative error = " << std::scientific << std::setprecision(6)
                       << finalError << "\n\n";
            } else {
                output << "Failed to converge within " << maxIter << " iterations\n";
            }
            break;
        case 1:
            converged = gaussSeidel(A, b, x, epsilon, maxIter, actualIter, finalError);
            if (converged) {
                output << "Iterations converged in: " << actualIter << " iterations\n";
                output << "Iterative error = " << std::scientific << std::setprecision(6)
                       << finalError << "\n\n";
            } else {
                output << "Failed to converge within " << maxIter << " iterations\n";
            }
            break;
        case 2:
            // Validate omega for SOR method
                if (omega <= 0.0 || omega >= 2.0) {
                    output << "Error: SOR relaxation parameter must be in range (0,2)\n";
                    return 1;
                }
            converged = sor(A, b, x, omega, epsilon, maxIter, actualIter, finalError);
            if (converged) {
                output << "Iterations converged in: " << actualIter << " iterations\n";
                output << "Iterative error = " << std::scientific << std::setprecision(6)
                       << finalError << "\n\n";
            } else {
                output << "Failed to converge within " << maxIter << " iterations\n";
            }
            break;
    }

    auto RelaxationEndTime = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> RelaxationDuration = RelaxationEndTime - RelaxationStartTime;

    auto LUPStartTime = std::chrono::high_resolution_clock::now();

    // Perform LUP factorization
    if (!lupFactorize(A, L, U, P)) {
        output << "Error: LUP Factorization failed.\n";
        return 1;
    }

    applyPermutationMatrix(P, b, Pb);
    forwardSubstitution(L, Pb, y);
    backSubstitution(U, y, x_lup);

    auto LUPEndTime = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> LUPDuration = LUPEndTime - LUPStartTime;
    auto IOStartTime2 = std::chrono::high_resolution_clock::now();

    // Print Relaxation solution vector and reference LUP solution vector
    output << "Solution Vector, x, Obtained using Relaxation Method:\n";
    for (int i = 0; i < n; i++) {
        output<< "x[" << i + 1 << "] = "  << std::scientific << std::setprecision(8) << x[i] << "\n";
    }
    output << "\n";

    // Print solution vector
    output << "Solution Vector, x_lup, Obtained using LUP Method:\n";
    for (int i = 0; i < n; ++i) {
        output << "x[" << i + 1 << "] = " << std::setprecision(8) << x_lup[i] << "\n";
    }
    output << "\n";

    // Calculate and print residual
    double maxResidual = calculateResidual(A, x, b);
    output << "Max residual for Relaxation Method = " << std::scientific << std::setprecision(8) << maxResidual << "\n";

    // Calculate and print maximum residual for LUP
    std::vector<double> residuals = calculateResiduals(A, x_lup, b);
    double maxLUPResidual = *std::max_element(residuals.begin(), residuals.end());
    output << "Max residual for LUP Method = : " << std::scientific << std::setprecision(8) << maxLUPResidual << "\n\n";

    auto IOEndTime2 = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> IODuration2 = IOEndTime2 - IOStartTime2;

    output << "Execution time for Relaxation Method(sec) = " << std::scientific << std::setprecision(8) << RelaxationDuration.count() << "\n";
    output << "Execution time for LUP Method(sec) = " << std::scientific << std::setprecision(8) << LUPDuration.count() << "\n";

    input.close();
    output.close();
    return 0;
}
