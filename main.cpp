//
// Steady State One-Speed Diffusion Equation Solver
// Version: 1.0
// Author: Hasibul H. Rasheeq
// Date: February 27, 2025
//

#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>
#include <chrono>
#include <ctime>
#include "DiffusionClass.cpp"
#include "LUP/MatrixVectorProduct.cpp"

// Function to get current date and time as a string
std::string getCurrentDateTime() {
    auto now = std::chrono::system_clock::now();
    std::time_t now_time = std::chrono::system_clock::to_time_t(now);

    char buf[100];
    std::strftime(buf, sizeof(buf), "%b %d %Y %H:%M:%S", std::localtime(&now_time));
    return std::string(buf);
}

// Function to calculate residuals
std::vector<double> calculateResiduals(const std::vector<std::vector<double>>& A,
                                    const std::vector<double>& x,
                                    const std::vector<double>& b) {
    // Calculate Ax
    std::vector<double> Ax = matrixVectorProduct(A, x);

    // Calculate residual r = Ax - b
    std::vector<double> residuals(b.size());
    for (size_t i = 0; i < b.size(); ++i) {
        residuals[i] = Ax[i] - b[i];
    }

    return residuals;
}

// Helper function to flatten 2D solution to 1D vector
std::vector<double> flattenSolution(const std::vector<std::vector<double>>& phi, int m, int n) {
    std::vector<double> flattened;
    flattened.reserve(m * n);

    for (int i = 1; i <= m; ++i) {
        for (int j = 1; j <= n; ++j) {
            flattened.push_back(phi[i][j]);
        }
    }
    return flattened;
}

// Function to estimate memory usage
size_t estimateMemoryUsage(int m, int n) {
    size_t totalSize = 0;

    // Memory for coefficient matrix (m*n × m*n)
    totalSize += sizeof(double) * m * n * m * n;

    // Memory for L, U, and P matrices
    totalSize += 3 * sizeof(double) * m * n * m * n;

    // Memory for vectors (solution, RHS, intermediate vectors)
    totalSize += 5 * sizeof(double) * m * n;

    // Memory for source term matrix
    totalSize += sizeof(double) * m * n;

    return totalSize;
}

int main() {
    DiffusionSolver solver;

    // Define the input file we'll process
    const std::string filepath = "../input.txt";

    // Define corresponding output file
    const std::string outputfile = "./output.txt";

    // Measure execution time
    auto startTime = std::chrono::high_resolution_clock::now();

    // Read input file
    if (!solver.readInput(filepath)) {
        return 1;
    }

    // Solve the system
    auto solution = solver.solve();

    // End timing
    auto endTime = std::chrono::high_resolution_clock::now();
    double executionTime = std::chrono::duration<double>(endTime - startTime).count();

    // Calculate residuals
    auto A = solver.buildMatrix();
    auto b = solver.buildRHS();
    auto x = flattenSolution(solution, solver.m, solver.n);
    auto residuals = calculateResiduals(A, x, b);

    // Find maximum residual
    double maxResidual = 0.0;
    for (const auto& res : residuals) {
        maxResidual = std::max(maxResidual, std::abs(res));
    }

    // Estimate memory usage
    size_t memoryUsage = estimateMemoryUsage(solver.m, solver.n);

    // Write output with all the required information
    solver.writeOutputWithDetails(solution, outputfile, executionTime, memoryUsage, residuals, maxResidual);

    // Output to console as well
    std::cout << "Solved system with grid size " << solver.m << "x" << solver.n << "\n";
    std::cout << "Execution time: " << executionTime << " seconds\n";
    std::cout << "Estimated memory usage: " << (memoryUsage / 1024.0 / 1024.0) << " MB\n";
    std::cout << "Maximum residual: " << std::scientific << maxResidual << std::endl;

    return 0;
}