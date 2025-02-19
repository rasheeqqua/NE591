#include <fstream>
#include <vector>
#include <iomanip>
#include <chrono>
#include "PointJacobi.cpp"
#include "GaussSeidel.cpp"
#include "SOR.cpp"

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
    std::ifstream input("../input.txt");
    std::ofstream output("output.txt");

    // Write header
    output << "NE 591 Inlab 06 Code\n";
    output << "Implemented by Hasibul H. Rasheeq, 02/14/25\n";
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

    // Solve system based on method selection
    auto start = std::chrono::high_resolution_clock::now();

    std::vector<double> x(n, 0.0);
    int actualIter;
    double finalError;
    bool converged;

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
            gaussSeidel(output);
            break;
        case 2:
            sor(output);
            break;
    }

    // Print solution vector
    output << "Last iterate of vector x:\n";
    for (int i = 0; i < n; i++) {
        output << std::scientific << std::setprecision(6) << x[i] << "\n";
    }
    output << "\n";

    // Calculate and print residual
    double maxResidual = calculateResidual(A, x, b);
    output << "Max residual = " << std::scientific << std::setprecision(6) << maxResidual << "\n\n";

    // Print execution time
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> duration = end - start;
    output << "Execution time (sec) = " << std::fixed << std::setprecision(6) << duration.count() << "\n";

    input.close();
    output.close();
    return 0;
}
