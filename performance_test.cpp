#include <iostream>
#include <fstream>
#include <vector>
#include <iomanip>
#include <chrono>
#include <cmath>
#include <random>
#include "matrix_modules/matrix_operations.h"
#include "matrix_modules/verify_positive_definite_matrix.h"
#include "LUP/LUP_solver.h"
#include "SOR/SOR_solver.h"
#include "CG/CG_solver.h"

// Function to generate a symmetric positive definite matrix of size n
std::vector<std::vector<double>> generateSPDMatrix(int n) {
    // Create a random matrix
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> dis(-1.0, 1.0);

    std::vector<std::vector<double>> A(n, std::vector<double>(n, 0.0));

    // Generate a random matrix
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            A[i][j] = dis(gen);
        }
    }

    // Make it symmetric: A = (A + A^T)/2
    for (int i = 0; i < n; i++) {
        for (int j = i+1; j < n; j++) {
            double avg = (A[i][j] + A[j][i]) / 2.0;
            A[i][j] = A[j][i] = avg;
        }
    }

    // Make it positive definite: A = A + n*I
    // This ensures diagonal dominance, which guarantees positive definiteness
    for (int i = 0; i < n; i++) {
        A[i][i] += n;
    }

    return A;
}

// Function to generate a random right-hand side vector of size n
std::vector<double> generateRandomVector(int n) {
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> dis(-10.0, 10.0);

    std::vector<double> b(n);
    for (int i = 0; i < n; i++) {
        b[i] = dis(gen);
    }

    return b;
}

// Function to write matrix to input file
void writeInputFile(const std::vector<std::vector<double>>& A, const std::vector<double>& b, int n, int method, double sorWeight, double epsilon, int maxIter) {
    std::string filename = "input_n" + std::to_string(n) + ".txt";
    std::ofstream inputFile(filename);

    if (!inputFile) {
        std::cerr << "Error opening input file: " << filename << std::endl;
        return;
    }

    // Write method flag and SOR weight
    inputFile << method << " " << sorWeight << std::endl;

    // Write stopping criterion and max iterations
    inputFile << epsilon << " " << maxIter << std::endl;

    // Write matrix order
    inputFile << n << std::endl;

    // Write matrix A
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            inputFile << A[i][j] << " ";
        }
        inputFile << std::endl;
    }

    // Write vector b
    for (int i = 0; i < n; i++) {
        inputFile << b[i] << " ";
    }
    inputFile << std::endl;

    inputFile.close();
}

int main() {
    std::ofstream resultsFile("performance_results.txt");

    if (!resultsFile) {
        std::cerr << "Error opening results file!" << std::endl;
        return 1;
    }

    // Parameters for solver
    double epsilon = 1e-4;
    int maxIter = 10000;
    double sorWeight = 1.5; // Optimized SOR weight

    // Matrix sizes to test
    std::vector<int> sizes = {32, 64, 128, 512, 1024};

    // Results storage
    std::vector<int> sorIterations;
    std::vector<int> cgIterations;
    std::vector<double> lupTimes;
    std::vector<double> sorTimes;
    std::vector<double> cgTimes;

    resultsFile << "Performance Comparison of LUP, SOR, and CG Methods" << std::endl;
    resultsFile << "==================================================" << std::endl << std::endl;
    resultsFile << "Stopping criterion: " << epsilon << std::endl;
    resultsFile << "Maximum iterations: " << maxIter << std::endl;
    resultsFile << "SOR weight: " << sorWeight << std::endl << std::endl;

    resultsFile << std::setw(10) << "Size" << std::setw(15) << "SOR Iter"
                << std::setw(15) << "CG Iter" << std::setw(20) << "LUP Time(ms)"
                << std::setw(20) << "SOR Time(ms)" << std::setw(20) << "CG Time(ms)" << std::endl;
    resultsFile << std::string(100, '-') << std::endl;

    for (int n : sizes) {
        std::cout << "Testing matrix size: " << n << std::endl;

        // Generate SPD matrix and RHS vector
        auto A = generateSPDMatrix(n);
        auto b = generateRandomVector(n);

        // Write input to file for each method
        writeInputFile(A, b, n, 0, sorWeight, epsilon, maxIter); // For LUP

        // Run LUP method
        std::vector<double> x_lup(n, 0.0);
        double maxResidual;

        auto lup_start = std::chrono::high_resolution_clock::now();
        bool lup_success = solveLUP(A, b, x_lup, maxResidual);
        auto lup_end = std::chrono::high_resolution_clock::now();

        std::chrono::duration<double> lup_elapsed_sec = lup_end - lup_start;
        double lup_elapsed_ms = lup_elapsed_sec.count() * 1000.0;
        lupTimes.push_back(lup_elapsed_ms);

        // Run SOR method
        std::vector<double> x_sor(n, 0.0);
        int sor_iterations;
        double sor_residual;

        auto sor_start = std::chrono::high_resolution_clock::now();
        bool sor_converged = solveSOR(A, b, x_sor, sorWeight, epsilon, maxIter, sor_iterations, sor_residual);
        auto sor_end = std::chrono::high_resolution_clock::now();

        std::chrono::duration<double> sor_elapsed_sec = sor_end - sor_start;
        double sor_elapsed_ms = sor_elapsed_sec.count() * 1000.0;
        sorTimes.push_back(sor_elapsed_ms);
        sorIterations.push_back(sor_iterations);

        // Run CG method
        std::vector<double> x_cg(n, 0.0);
        int cg_iterations;
        double cg_residual;

        auto cg_start = std::chrono::high_resolution_clock::now();
        bool cg_converged = solveCG(A, b, x_cg, epsilon, maxIter, cg_iterations, cg_residual);
        auto cg_end = std::chrono::high_resolution_clock::now();

        std::chrono::duration<double> cg_elapsed_sec = cg_end - cg_start;
        double cg_elapsed_ms = cg_elapsed_sec.count() * 1000.0;
        cgTimes.push_back(cg_elapsed_ms);
        cgIterations.push_back(cg_iterations);

        // Write result for this size
        resultsFile << std::setw(10) << n
                    << std::setw(15) << sor_iterations
                    << std::setw(15) << cg_iterations
                    << std::fixed << std::setprecision(8)
                    << std::setw(20) << lup_elapsed_ms
                    << std::setw(20) << sor_elapsed_ms
                    << std::setw(20) << cg_elapsed_ms
                    << std::endl;

        // Record detailed results for each size
        resultsFile << std::endl << "Detailed results for n = " << n << ":" << std::endl;
        resultsFile << "----------------------------------------" << std::endl;

        // LUP details
        resultsFile << "LUP Method:" << std::endl;
        resultsFile << "  - Maximum absolute residual: " << std::scientific << std::setprecision(8) << maxResidual << std::endl;
        resultsFile << "  - Direct solution (no iterations)" << std::endl;
        resultsFile << "  - Execution time: " << std::fixed << std::setprecision(8) << lup_elapsed_ms << " milliseconds" << std::endl << std::endl;

        // SOR details
        resultsFile << "SOR Method:" << std::endl;
        resultsFile << "  - Relaxation parameter (omega): " << sorWeight << std::endl;
        resultsFile << "  - Iterations performed: " << sor_iterations << std::endl;
        resultsFile << "  - Final residual norm: " << std::scientific << std::setprecision(8) << sor_residual << std::endl;
        resultsFile << "  - Execution time: " << std::fixed << std::setprecision(8) << sor_elapsed_ms << " milliseconds" << std::endl << std::endl;

        // CG details
        resultsFile << "CG Method:" << std::endl;
        resultsFile << "  - Iterations performed: " << cg_iterations << std::endl;
        resultsFile << "  - Final residual norm: " << std::scientific << std::setprecision(8) << cg_residual << std::endl;
        resultsFile << "  - Execution time: " << std::fixed << std::setprecision(8) << cg_elapsed_ms << " milliseconds" << std::endl << std::endl;

        // Solution comparison
        double solution_diff_lup_cg = 0.0;
        double solution_diff_lup_sor = 0.0;

        for (int i = 0; i < n; i++) {
            solution_diff_lup_cg = std::max(solution_diff_lup_cg, std::abs(x_lup[i] - x_cg[i]));
            solution_diff_lup_sor = std::max(solution_diff_lup_sor, std::abs(x_lup[i] - x_sor[i]));
        }

        resultsFile << "Solution Comparison:" << std::endl;
        resultsFile << "  - Maximum difference between LUP and CG solutions: "
                   << std::scientific << std::setprecision(8) << solution_diff_lup_cg << std::endl;
        resultsFile << "  - Maximum difference between LUP and SOR solutions: "
                   << std::scientific << std::setprecision(8) << solution_diff_lup_sor << std::endl;
        resultsFile << std::endl << std::string(50, '-') << std::endl << std::endl;
    }

    resultsFile.close();

    std::cout << "Performance testing complete. Results saved to performance_results.txt" << std::endl;
    std::cout << "Input matrices saved to input_n*.txt files" << std::endl;

    return 0;
}