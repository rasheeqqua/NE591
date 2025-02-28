//
// Simplified Verification Test for Diffusion Equation Solver - Milestone 2
// Author: Hasibul H. Rasheeq
// Date: February 28, 2025
//
// This program tests the diffusion equation solver by verifying
// reflective and rotational symmetry properties.
//

#include <iostream>
#include <fstream>
#include <iomanip>
#include <vector>
#include <string>
#include <cmath>
#include <chrono>
#include <string>
#include "DiffusionSolver.cpp"

// Class to handle the verification tests with hard-coded values
class VerificationTest {
private:
    // Test parameters - hard-coded based on the examples
    const double a = 60.0;
    const double b = 60.0;
    const int m = 4;
    const int n = 4;
    const double D = 0.142;
    const double sigma_a = 0.02220;
    const int maxIter = 100;
    const double tolerance = 1e-4;
    const double omega = 1.15;

public:
    // Run a single test with the given method and source distribution
    std::vector<std::vector<double>> runTest(int method, const std::vector<std::vector<double>>& source) {
        // Create solver instance
        DiffusionSolver solver;

        // Set parameters directly in the solver
        solver.setParameters(method, maxIter, tolerance, omega, a, b, m, n, D, sigma_a, source);

        // Variables to store results
        int iterations;
        double finalError;
        bool converged;

        // Solve the system
        auto solution = solver.solve(iterations, finalError, converged);

        if (!converged) {
            std::cerr << "Warning: Method " << method << " did not converge!" << std::endl;
        } else {
            std::cout << "Method " << method << " converged in " << iterations << " iterations." << std::endl;
        }

        return solution;
    }

    // Test reflective symmetry
    void testReflectiveSymmetry(int method) {
        std::cout << "\nRunning Reflective Symmetry Test for method " << method << std::endl;

        // Create source distributions for reflective symmetry test
        std::vector<std::vector<double>> source1 = {
            {0.0, 0.0, 0.0, 0.0},
            {0.0, 10.0, 0.0, 0.0},
            {0.0, 10.0, 0.0, 0.0},
            {0.0, 0.0, 0.0, 0.0}
        };

        std::vector<std::vector<double>> source2 = {
            {0.0, 0.0, 0.0, 0.0},
            {0.0, 0.0, 10.0, 0.0},
            {0.0, 0.0, 10.0, 0.0},
            {0.0, 0.0, 0.0, 0.0}
        };

        // Run tests
        auto solution1 = runTest(method, source1);
        auto solution2 = runTest(method, source2);

        // Write results to file
        std::string filename = "reflective_test" + std::to_string(method) + "_results.txt";
        writeComparisonResults(filename, "Reflective Symmetry Test", method, source1, source2, solution1, solution2);
    }

    // Test rotational symmetry
    void testRotationalSymmetry(int method) {
        std::cout << "\nRunning Rotational Symmetry Test for method " << method << std::endl;

        // Create source distributions for rotational symmetry test
        std::vector<std::vector<double>> source1 = {
            {0.0, 0.0, 0.0, 0.0},
            {0.0, 10.0, 10.0, 0.0},
            {0.0, 0.0, 0.0, 0.0},
            {0.0, 0.0, 0.0, 0.0}
        };

        std::vector<std::vector<double>> source2 = {
            {0.0, 0.0, 0.0, 0.0},
            {0.0, 10.0, 0.0, 0.0},
            {0.0, 10.0, 0.0, 0.0},
            {0.0, 0.0, 0.0, 0.0}
        };

        // Run tests
        auto solution1 = runTest(method, source1);
        auto solution2 = runTest(method, source2);

        // Write results to file
        std::string filename = "rotational_test" + std::to_string(method) + "_results.txt";
        writeComparisonResults(filename, "Rotational Symmetry Test", method, source1, source2, solution1, solution2);
    }

    // Write comparison results in a juxtaposed format
    void writeComparisonResults(const std::string& filename,
                               const std::string& testName,
                               int method,
                               const std::vector<std::vector<double>>& source1,
                               const std::vector<std::vector<double>>& source2,
                               const std::vector<std::vector<double>>& solution1,
                               const std::vector<std::vector<double>>& solution2) {
        std::ofstream file(filename);
        if (!file.is_open()) {
            std::cerr << "Error opening file: " << filename << std::endl;
            return;
        }

        // Method name
        std::string methodName;
        switch (method) {
            case 0: methodName = "LUP (Direct)"; break;
            case 1: methodName = "Point Jacobi"; break;
            case 2: methodName = "Gauss-Seidel"; break;
            case 3: methodName = "SOR"; break;
            default: methodName = "Unknown";
        }

        // Write header
        file << "Diffusion Equation Solver - " << testName << std::endl;
        file << "Method: " << methodName << std::endl;
        file << "Rectangle: " << a << " x " << b << " cm, Grid: " << m << " x " << n << std::endl;
        file << "Parameters: D = " << D << " cm, Σₐ = " << sigma_a << " cm⁻¹" << std::endl;
        file << "Iterative parameters: max iterations = " << maxIter << ", tolerance = " << tolerance;
        if (method == 3) {
            file << ", ω = " << omega;
        }
        file << std::endl << std::endl;

        // Write source distributions
        file << "Source Distribution 1:" << std::setw(36) << "Source Distribution 2:" << std::endl;
        for (int i = 0; i < m; i++) {
            for (int j = 0; j < n; j++) {
                file << std::setw(10) << source1[i][j];
            }
            file << std::setw(10) << " | ";
            for (int j = 0; j < n; j++) {
                file << std::setw(10) << source2[i][j];
            }
            file << std::endl;
        }
        file << std::endl;

        // Write solutions side by side
        file << "Solution 1 (including boundaries):" << std::setw(30) << "Solution 2 (including boundaries):" << std::endl;
        for (int i = 0; i < solution1.size(); i++) {
            for (int j = 0; j < solution1[i].size(); j++) {
                file << std::scientific << std::setw(15) << std::setprecision(6) << solution1[i][j];
            }
            file << " | ";
            for (int j = 0; j < solution2[i].size(); j++) {
                file << std::scientific << std::setw(15) << std::setprecision(6) << solution2[i][j];
            }
            file << std::endl;
        }

        file.close();
        std::cout << "Comparison results written to " << filename << std::endl;
    }
};

int main() {
    std::cout << "Simplified Verification Test for Diffusion Equation Solver - Milestone 2\n";
    std::cout << "================================================================\n\n";

    VerificationTest tester;

    // Test reflective symmetry for all methods
    for (int method = 0; method <= 3; ++method) {
        tester.testReflectiveSymmetry(method);
    }

    // Test rotational symmetry for all methods
    for (int method = 0; method <= 3; ++method) {
        tester.testRotationalSymmetry(method);
    }

    std::cout << "\nAll tests completed. Results are written to separate files for each test and method.\n";
    std::cout << "Please examine the output files to verify the symmetry properties.\n";

    return 0;
}