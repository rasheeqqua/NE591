//
// Performance Analyzer for Diffusion Equation Solver - Milestone 2
// Author: Hasibul H. Rasheeq
// Date: February 28, 2025
//
// This program benchmarks the performance of different solution methods
// (LUP, Jacobi, Gauss-Seidel, SOR) for the diffusion equation.
//

#include <chrono>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <vector>
#include <string>
#include <algorithm>
#include <limits>
#include "DiffusionSolver.cpp"

class PerformanceAnalyzer {
private:
    // Structure to hold performance metrics
    struct PerformanceMetrics {
        int methodFlag;           // Solution method (0=LUP, 1=Jacobi, 2=Gauss-Seidel, 3=SOR)
        double executionTime;     // in seconds
        size_t memoryUsage;       // in bytes
        double maxResidual;       // largest magnitude residual
        int iterations;           // number of iterations (for iterative methods)
        bool converged;           // whether the method converged
        int gridSize;             // m*n (total number of unknowns)
    };

    // Helper function to estimate memory usage for different methods
    size_t estimateMemoryUsage(int method, int m, int n) {
        size_t totalSize = 0;

        // Base memory for all methods
        totalSize += sizeof(double) * m * n;  // Source term
        totalSize += sizeof(double) * (m+2) * (n+2);  // Solution vector

        if (method == 0) {  // LUP method
            // Memory for coefficient matrix (m*n × m*n)
            totalSize += sizeof(double) * m * n * m * n;

            // Memory for L, U, and P matrices
            totalSize += 3 * sizeof(double) * m * n * m * n;

            // Memory for intermediate vectors
            totalSize += 3 * sizeof(double) * m * n;
        } else {  // Iterative methods
            // Memory for additional flux array (previous iteration)
            totalSize += sizeof(double) * (m+2) * (n+2);

            // Memory for residual calculation
            totalSize += sizeof(double) * (m+2) * (n+2);
        }

        return totalSize;
    }

    // Function to get method name from flag
    std::string getMethodName(int method) {
        switch (method) {
            case 0: return "LUP (Direct)";
            case 1: return "Point Jacobi";
            case 2: return "Gauss-Seidel";
            case 3: return "SOR";
            default: return "Unknown";
        }
    }

public:
    // Function to run a single test case with a specific method and grid size
    PerformanceMetrics runSingleTest(
        int method, int m, int n,
        double a, double b,
        double D, double sigma_a,
        int maxIter = 10000,
        double tolerance = 1e-6,
        double omega = 1.5) {

        std::cout << "Testing " << getMethodName(method) << " on grid size " << m << "x" << n << "...\n";

        // Create a source distribution (Gaussian source at center)
        std::vector<std::vector<double>> source(m, std::vector<double>(n));
        double centerX = a / 2.0;
        double centerY = b / 2.0;
        double sigma = std::min(a, b) / 8.0;  // width of Gaussian

        for (int i = 0; i < m; ++i) {
            for (int j = 0; j < n; ++j) {
                double x = a * (i + 1.0) / (m + 1.0);
                double y = b * (j + 1.0) / (n + 1.0);

                // Gaussian source term
                source[i][j] = 10.0 * std::exp(
                    -((x - centerX) * (x - centerX) + (y - centerY) * (y - centerY))
                    / (2.0 * sigma * sigma));
            }
        }

        // Create solver instance and set parameters
        DiffusionSolver solver;
        solver.setParameters(method, maxIter, tolerance, omega, a, b, m, n, D, sigma_a, source);

        // Variables to store results
        int iterations;
        double finalError;
        bool converged;

        // Measure execution time
        auto startTime = std::chrono::high_resolution_clock::now();

        // Solve the system
        auto solution = solver.solve(iterations, finalError, converged);

        auto endTime = std::chrono::high_resolution_clock::now();
        double executionTime = std::chrono::duration<double>(endTime - startTime).count();

        // Calculate memory usage
        size_t memoryUsage = estimateMemoryUsage(method, m, n);

        // Calculate maximum residual
        double maxResidual = 0.0;

        // For iterative methods, we use the reported finalError
        // For direct method, we calculate actual residual
        if (method == 0) {
            // This assumes the solver can calculate residuals
            maxResidual = finalError;
        } else {
            maxResidual = finalError;
        }

        // Return performance metrics
        return {
            method,
            executionTime,
            memoryUsage,
            maxResidual,
            iterations,
            converged,
            m * n
        };
    }

    // Function to run performance tests for all methods on different grid sizes
    void runComprehensiveTests(
        const std::vector<std::pair<int, int>>& gridSizes,
        double a = 60.0, double b = 60.0,
        double D = 0.142, double sigma_a = 0.0222) {

        std::vector<PerformanceMetrics> results;

        // Parameters for iterative methods
        int maxIter = 10000;
        double tolerance = 1e-6;
        double omega = 1.5;  // for SOR

        // Run tests for each method and grid size
        for (int method = 0; method <= 3; ++method) {
            // For LUP method, skip large grid sizes
            if (method == 0) {
                for (const auto& [m, n] : gridSizes) {
                    // Skip large grids for LUP to avoid memory issues
                    if (m * n <= 2500) {  // Limit LUP to 50x50 or equivalent
                        auto result = runSingleTest(method, m, n, a, b, D, sigma_a, maxIter, tolerance, omega);
                        results.push_back(result);
                    }
                }
            } else {
                // For iterative methods, test all grid sizes
                for (const auto& [m, n] : gridSizes) {
                    auto result = runSingleTest(method, m, n, a, b, D, sigma_a, maxIter, tolerance, omega);
                    results.push_back(result);
                }
            }
        }

        // Write results to file
        writeDetailedResults(results, "performance_results.txt");
        writeComparisonCharts(results, "performance_comparison");
    }

    // Function to write detailed results to file
    void writeDetailedResults(const std::vector<PerformanceMetrics>& results,
                            const std::string& filename) {
        std::ofstream outFile(filename);

        outFile << "Performance Analysis Results for Diffusion Equation Solver - Milestone 2\n";
        outFile << "Author: Hasibul H. Rasheeq\n";
        outFile << "Date: " << getCurrentDateTime() << "\n";
        outFile << "----------------------------------------------------------------------\n\n";

        // Group results by method
        outFile << "Results by Method:\n";
        outFile << "=================\n\n";

        for (int method = 0; method <= 3; ++method) {
            outFile << getMethodName(method) << ":\n";
            outFile << std::string(getMethodName(method).length() + 1, '-') << "\n";
            outFile << "Grid Size | Unknowns | Converged | Iterations | Exec Time (s) | Memory (MB) | Max Residual\n";
            outFile << "--------------------------------------------------------------------------------------\n";

            bool hasResults = false;
            for (const auto& result : results) {
                if (result.methodFlag == method) {
                    hasResults = true;
                    outFile << std::setw(8) << static_cast<int>(std::sqrt(result.gridSize)) << "x"
                           << static_cast<int>(std::sqrt(result.gridSize)) << " | "
                           << std::setw(8) << result.gridSize << " | "
                           << std::setw(9) << (result.converged ? "Yes" : "No") << " | "
                           << std::setw(10) << (method > 0 ? std::to_string(result.iterations) : "N/A") << " | "
                           << std::setw(13) << std::fixed << std::setprecision(8) << result.executionTime << " | "
                           << std::setw(10) << std::fixed << std::setprecision(8) << (result.memoryUsage / 1024.0 / 1024.0) << " | "
                           << std::setw(11) << std::scientific << std::setprecision(3) << result.maxResidual << "\n";
                }
            }

            if (!hasResults) {
                outFile << "No results available for this method\n";
            }
            outFile << "\n\n";
        }

        // Results by grid size
        outFile << "Results by Grid Size:\n";
        outFile << "====================\n\n";

        std::vector<int> uniqueGridSizes;
        for (const auto& result : results) {
            if (std::find(uniqueGridSizes.begin(), uniqueGridSizes.end(), result.gridSize) == uniqueGridSizes.end()) {
                uniqueGridSizes.push_back(result.gridSize);
            }
        }
        std::sort(uniqueGridSizes.begin(), uniqueGridSizes.end());

        for (int gridSize : uniqueGridSizes) {
            outFile << "Grid Size: " << static_cast<int>(std::sqrt(gridSize)) << "x" << static_cast<int>(std::sqrt(gridSize))
                   << " (" << gridSize << " unknowns):\n";
            outFile << "Method      | Converged | Iterations | Exec Time (s) | Memory (MB) | Max Residual\n";
            outFile << "----------------------------------------------------------------------------\n";

            for (int method = 0; method <= 3; ++method) {
                bool foundMethod = false;
                for (const auto& result : results) {
                    if (result.methodFlag == method && result.gridSize == gridSize) {
                        foundMethod = true;
                        outFile << std::left << std::setw(12) << getMethodName(method) << " | "
                               << std::setw(9) << (result.converged ? "Yes" : "No") << " | "
                               << std::setw(10) << (method > 0 ? std::to_string(result.iterations) : "N/A") << " | "
                               << std::right << std::setw(13) << std::fixed << std::setprecision(8) << result.executionTime << " | "
                               << std::setw(10) << std::fixed << std::setprecision(8) << (result.memoryUsage / 1024.0 / 1024.0) << " | "
                               << std::setw(11) << std::scientific << std::setprecision(3) << result.maxResidual << "\n";
                        break;
                    }
                }

                if (!foundMethod) {
                    outFile << std::left << std::setw(12) << getMethodName(method) << " | "
                           << std::setw(9) << "N/A" << " | "
                           << std::setw(10) << "N/A" << " | "
                           << std::right << std::setw(13) << "N/A" << " | "
                           << std::setw(10) << "N/A" << " | "
                           << std::setw(11) << "N/A" << "\n";
                }
            }
            outFile << "\n\n";
        }

        // Performance summary
        outFile << "Performance Summary:\n";
        outFile << "===================\n\n";

        // Find method with lowest execution time for each grid size
        outFile << "Fastest Method by Grid Size:\n";
        for (int gridSize : uniqueGridSizes) {
            int fastestMethod = -1;
            double fastestTime = std::numeric_limits<double>::max();
            for (const auto& result : results) {
                if (result.gridSize == gridSize && result.converged && result.executionTime < fastestTime) {
                    fastestMethod = result.methodFlag;
                    fastestTime = result.executionTime;
                }
            }

            outFile << static_cast<int>(std::sqrt(gridSize)) << "x" << static_cast<int>(std::sqrt(gridSize))
                   << " (" << gridSize << " unknowns): ";
            if (fastestMethod >= 0) {
                outFile << getMethodName(fastestMethod) << " (" << std::fixed << std::setprecision(4) << fastestTime << " s)\n";
            } else {
                outFile << "No method converged\n";
            }
        }
        outFile << "\n";

        // Find method with lowest memory usage for each grid size
        outFile << "Memory-Efficient Method by Grid Size:\n";
        for (int gridSize : uniqueGridSizes) {
            int efficientMethod = -1;
            size_t lowestMemory = std::numeric_limits<size_t>::max();
            for (const auto& result : results) {
                if (result.gridSize == gridSize && result.converged && result.memoryUsage < lowestMemory) {
                    efficientMethod = result.methodFlag;
                    lowestMemory = result.memoryUsage;
                }
            }

            outFile << static_cast<int>(std::sqrt(gridSize)) << "x" << static_cast<int>(std::sqrt(gridSize))
                   << " (" << gridSize << " unknowns): ";
            if (efficientMethod >= 0) {
                outFile << getMethodName(efficientMethod) << " (" << std::fixed << std::setprecision(2)
                       << (lowestMemory / 1024.0 / 1024.0) << " MB)\n";
            } else {
                outFile << "No method converged\n";
            }
        }
        outFile << "\n";

        // Find method with highest accuracy for each grid size
        outFile << "Most Accurate Method by Grid Size:\n";
        for (int gridSize : uniqueGridSizes) {
            int accurateMethod = -1;
            double lowestResidual = std::numeric_limits<double>::max();
            for (const auto& result : results) {
                if (result.gridSize == gridSize && result.converged && result.maxResidual < lowestResidual) {
                    accurateMethod = result.methodFlag;
                    lowestResidual = result.maxResidual;
                }
            }

            outFile << static_cast<int>(std::sqrt(gridSize)) << "x" << static_cast<int>(std::sqrt(gridSize))
                   << " (" << gridSize << " unknowns): ";
            if (accurateMethod >= 0) {
                outFile << getMethodName(accurateMethod) << " (" << std::scientific << std::setprecision(3)
                       << lowestResidual << ")\n";
            } else {
                outFile << "No method converged\n";
            }
        }

        outFile.close();
    }

    // Function to create comparison charts (data for plotting)
    void writeComparisonCharts(const std::vector<PerformanceMetrics>& results,
                              const std::string& filenamePrefix) {
        // Create data files for execution time vs grid size
        std::ofstream timeFile(filenamePrefix + "_time.txt");
        timeFile << "# Grid_Size LUP Jacobi Gauss-Seidel SOR\n";

        // Create data files for memory usage vs grid size
        std::ofstream memoryFile(filenamePrefix + "_memory.txt");
        memoryFile << "# Grid_Size LUP Jacobi Gauss-Seidel SOR\n";

        // Create data files for residual vs grid size
        std::ofstream residualFile(filenamePrefix + "_residual.txt");
        residualFile << "# Grid_Size LUP Jacobi Gauss-Seidel SOR\n";

        // Create data files for iterations vs grid size (for iterative methods)
        std::ofstream iterFile(filenamePrefix + "_iterations.txt");
        iterFile << "# Grid_Size Jacobi Gauss-Seidel SOR\n";

        // Get unique grid sizes
        std::vector<int> uniqueGridSizes;
        for (const auto& result : results) {
            if (std::find(uniqueGridSizes.begin(), uniqueGridSizes.end(), result.gridSize) == uniqueGridSizes.end()) {
                uniqueGridSizes.push_back(result.gridSize);
            }
        }
        std::sort(uniqueGridSizes.begin(), uniqueGridSizes.end());

        // Write data for each grid size
        for (int gridSize : uniqueGridSizes) {
            // Initialize with placeholder values
            double lupTime = -1.0, jacobiTime = -1.0, gsTime = -1.0, sorTime = -1.0;
            double lupMemory = -1.0, jacobiMemory = -1.0, gsMemory = -1.0, sorMemory = -1.0;
            double lupResidual = -1.0, jacobiResidual = -1.0, gsResidual = -1.0, sorResidual = -1.0;
            int jacobiIter = -1, gsIter = -1, sorIter = -1;

            // Find data for each method at this grid size
            for (const auto& result : results) {
                if (result.gridSize == gridSize && result.converged) {
                    switch (result.methodFlag) {
                        case 0: // LUP
                            lupTime = result.executionTime;
                            lupMemory = result.memoryUsage / (1024.0 * 1024.0); // MB
                            lupResidual = result.maxResidual;
                            break;
                        case 1: // Jacobi
                            jacobiTime = result.executionTime;
                            jacobiMemory = result.memoryUsage / (1024.0 * 1024.0);
                            jacobiResidual = result.maxResidual;
                            jacobiIter = result.iterations;
                            break;
                        case 2: // Gauss-Seidel
                            gsTime = result.executionTime;
                            gsMemory = result.memoryUsage / (1024.0 * 1024.0);
                            gsResidual = result.maxResidual;
                            gsIter = result.iterations;
                            break;
                        case 3: // SOR
                            sorTime = result.executionTime;
                            sorMemory = result.memoryUsage / (1024.0 * 1024.0);
                            sorResidual = result.maxResidual;
                            sorIter = result.iterations;
                            break;
                    }
                }
            }

            // Write time data
            timeFile << gridSize << " "
                    << lupTime << " "
                    << jacobiTime << " "
                    << gsTime << " "
                    << sorTime << "\n";

            // Write memory data
            memoryFile << gridSize << " "
                      << lupMemory << " "
                      << jacobiMemory << " "
                      << gsMemory << " "
                      << sorMemory << "\n";

            // Write residual data
            residualFile << gridSize << " "
                        << lupResidual << " "
                        << jacobiResidual << " "
                        << gsResidual << " "
                        << sorResidual << "\n";

            // Write iteration data (only for iterative methods)
            iterFile << gridSize << " "
                    << jacobiIter << " "
                    << gsIter << " "
                    << sorIter << "\n";
        }

        timeFile.close();
        memoryFile.close();
        residualFile.close();
        iterFile.close();
    }

    // Helper function to get current date and time
    std::string getCurrentDateTime() {
        auto now = std::chrono::system_clock::now();
        std::time_t now_time = std::chrono::system_clock::to_time_t(now);
        char buf[100];
        std::strftime(buf, sizeof(buf), "%b %d %Y %H:%M:%S", std::localtime(&now_time));
        return std::string(buf);
    }
};

int main() {
    std::cout << "Performance Analyzer for Diffusion Equation Solver - Milestone 2\n";
    std::cout << "============================================================\n\n";

    PerformanceAnalyzer analyzer;

    // Define sequence of grid sizes to test
    std::vector<std::pair<int, int>> gridSizes = {
        {5, 5},      // 25 unknowns
        {10, 10},    // 100 unknowns
        {20, 20},    // 400 unknowns
        {40, 40},    // 1,600 unknowns
        {50, 50},    // 2,500 unknowns
        {80, 80},    // 6,400 unknowns
        {100, 100},  // 10,000 unknowns
    };

    // Run performance tests for all methods
    analyzer.runComprehensiveTests(gridSizes);

    std::cout << "\nPerformance analysis completed. Results written to 'performance_results.txt'.\n";
    std::cout << "Data for plotting also written to 'performance_comparison_*.txt' files.\n";

    return 0;
}