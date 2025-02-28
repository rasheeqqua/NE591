//
// Created by Hasibul H. Rasheeq on 02/13/25.
//

#include <chrono>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <memory>
#include "../DiffusionSolver.cpp"
#include "MatrixVectorProduct.cpp"

class PerformanceAnalyzer {
private:
    // Structure to hold performance metrics
    struct PerformanceMetrics {
        double executionTime;      // in seconds
        size_t memoryUsage;        // in bytes
        double maxResidual;        // largest magnitude residual
        int gridSize;             // m*n (total number of unknowns)
    };

    // Helper function to estimate memory usage for a given problem size
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

    // Helper function to calculate residuals
    std::vector<double> calculateResiduals(const std::vector<std::vector<double>>& A,
                                         const std::vector<double>& x,
                                         const std::vector<double>& b) {
        // Calculate Ax
        std::vector<double> Ax = matrixVectorProduct(A, x);

        // Calculate residual r = Ax - b
        std::vector<double> residuals(b.size());
        for (size_t i = 0; i < b.size(); ++i) {
            residuals[i] = std::abs(Ax[i] - b[i]);
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

public:
    // Function to run performance tests for different grid sizes
    std::vector<PerformanceMetrics> runPerformanceTests(
            const std::vector<std::pair<int, int>>& gridSizes,
            double a = 10.0, double b = 10.0,
            double D = 0.5, double sigma_a = 0.1) {

        std::vector<PerformanceMetrics> results;

        for (const auto& [m, n] : gridSizes) {
            std::cout << "Testing grid size " << m << "x" << n << "...\n";

            // Create input file for current grid size
            createInputFile(m, n, a, b, D, sigma_a);

            // Create solver instance
            DiffusionSolver solver;

            // Measure execution time
            auto startTime = std::chrono::high_resolution_clock::now();

            // Run solver
            solver.readInput("temp_input.txt");
            auto solution = solver.solve();

            auto endTime = std::chrono::high_resolution_clock::now();
            double executionTime = std::chrono::duration<double>(endTime - startTime).count();

            // Calculate memory usage
            size_t memoryUsage = estimateMemoryUsage(m, n);

            // Calculate maximum residual
            auto A = solver.buildMatrix();  // Assuming we make this method public
            auto b = solver.buildRHS();     // Assuming we make this method public
            auto x = flattenSolution(solution, m, n);
            auto residuals = calculateResiduals(A, x, b);
            double maxResidual = *std::max_element(residuals.begin(), residuals.end());

            // Store results
            results.push_back({
                executionTime,
                memoryUsage,
                maxResidual,
                m * n
            });
        }

        return results;
    }

    // Function to create input file for a given grid size
    void createInputFile(int m, int n, double a, double b, double D, double sigma_a) {
        std::ofstream inputFile("temp_input.txt");

        // Write header information
        inputFile << "0\n";           // flag for direct solution
        inputFile << a << " " << b << "\n";
        inputFile << m << " " << n << "\n";
        inputFile << D << " " << sigma_a << "\n";

        // Create a simple source distribution (Gaussian source at center)
        double centerX = a / 2.0;
        double centerY = b / 2.0;
        double sigma = std::min(a, b) / 8.0;  // width of Gaussian

        for (int i = 0; i < m; ++i) {
            for (int j = 0; j < n; ++j) {
                double x = a * (i + 1.0) / (m + 1.0);
                double y = b * (j + 1.0) / (n + 1.0);

                // Gaussian source term
                double q = 10.0 * std::exp(
                    -((x - centerX) * (x - centerX) + (y - centerY) * (y - centerY))
                    / (2.0 * sigma * sigma));

                inputFile << q << " ";
            }
            inputFile << "\n";
        }

        inputFile.close();
    }

    // Function to write results to file
    void writeResults(const std::vector<PerformanceMetrics>& results,
                     const std::string& filename) {
        std::ofstream outFile(filename);

        outFile << "Performance Analysis Results\n";
        outFile << "----------------------------\n\n";
        outFile << "Grid Size (m*n) | Execution Time (s) | Memory Usage (MB) | Max Residual\n";
        outFile << "------------------------------------------------------------\n";

        for (const auto& result : results) {
            outFile << std::setw(13) << result.gridSize << " | "
                   << std::setw(16) << std::fixed << std::setprecision(4)
                   << result.executionTime << " | "
                   << std::setw(15) << std::fixed << std::setprecision(2)
                   << (result.memoryUsage / 1024.0 / 1024.0) << " | "
                   << std::setw(11) << std::scientific << std::setprecision(3)
                   << result.maxResidual << "\n";
        }

        outFile.close();
    }
};

// Example usage:
int main() {
    PerformanceAnalyzer analyzer;

    // Define sequence of grid sizes to test
    std::vector<std::pair<int, int>> gridSizes = {
        {5, 5},      // 25 unknowns
        {10, 10},    // 100 unknowns
        {20, 20},    // 400 unknowns
        {40, 40},    // 1,600 unknowns
    };

    // Run performance tests
    auto results = analyzer.runPerformanceTests(gridSizes);

    // Write results to file
    analyzer.writeResults(results, "performance_results.txt");

    return 0;
}
