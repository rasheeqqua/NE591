#include <iostream>
#include <fstream>
#include <vector>
#include <iomanip>
#include <chrono>
#include <cmath>
#include "matrix_modules/matrix_operations.h"
#include "matrix_modules/verify_positive_definite_matrix.h"
#include "LUP/LUP_solver.h"
#include "SOR/SOR_solver.h"
#include "CG/CG_solver.h"
#include "PCG/PCG_solver.h"
#include "PI/power_iterations.h"

int main() {
    std::ifstream inFile("input.txt");
    std::ofstream outFile("output.txt");

    if (!inFile || !outFile) {
        std::cerr << "Error opening files!" << std::endl;
        return 1;
    }

    // Read method flag and SOR weight
    int methodFlag;
    double sorWeight;
    inFile >> methodFlag >> sorWeight;

    // Read stopping criterion and maximum iterations
    double epsilon;
    int maxIter;
    inFile >> epsilon >> maxIter;

    // Read matrix order
    int n;
    inFile >> n;

    // Read matrix A
    std::vector<std::vector<double>> A(n, std::vector<double>(n));
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            inFile >> A[i][j];
        }
    }

    // Read vector b (or initial guess for Power Iterations)
    std::vector<double> b(n);
    for (int i = 0; i < n; i++) {
        inFile >> b[i];
    }

    // Write header to output file
    outFile << "NE 591 - Outlab 11 Code" << std::endl;
    outFile << "Implemented by Hasibul Hossain Rasheeq, March 31, 2025" << std::endl;
    outFile << "------------------------------------------------------" << std::endl << std::endl;

    if (methodFlag == 4) {
        outFile << "Compute fundamental eigenvector with Power Iterations" << std::endl << std::endl;
    } else {
        // Check matrix properties
        bool symmetric = isSymmetric(A);
        bool diagonallyDominant = isDiagonallyDominant(A);

        if (!symmetric) {
            outFile << "Error: Matrix is not symmetric!" << std::endl;
            return 1;
        }

        outFile << "Solve Symmetric Positive Definite Matrix" << std::endl;
        outFile << "Equation with Linear Solvers" << std::endl << std::endl;

        outFile << "Matrix symmetry checked" << std::endl;
        if (diagonallyDominant) {
            outFile << "Matrix is diagonally dominant" << std::endl;
        } else {
            outFile << "Warning: Matrix is not diagonally dominant" << std::endl;
        }
        outFile << "User must ensure it is positive definite" << std::endl;
        outFile << "-----------------------------------------" << std::endl << std::endl;
    }

    // Echo input data to output file
    outFile << "stopping criterion on residual norm = " << std::scientific << std::setprecision(2) << epsilon << std::endl;
    outFile << "maximum iterations = " << maxIter << std::endl;
    outFile << "matrix is of order: " << n << std::endl;
    outFile << "Matrix A:" << std::endl;

    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            outFile << std::scientific << std::setprecision(2) << A[i][j] << " ";
        }
        outFile << std::endl;
    }

    outFile << std::endl;

    if (methodFlag == 4) {
        // For Power Iterations, b is the initial guess
        outFile << "Initial guess:" << std::endl;
    } else {
        outFile << "RHS vector b:" << std::endl;
    }

    for (int i = 0; i < n; i++) {
        outFile << std::scientific << std::setprecision(2) << b[i] << " ";
    }
    outFile << std::endl << std::endl;

    // Solve the system based on method flag
    std::vector<double> x(n, 0.0);

    // Start the timer
    std::chrono::duration<double> elapsed;
    auto start = std::chrono::high_resolution_clock::now();

    if (methodFlag == 0) {
        // LUP method
        double maxResidual;
        bool success = solveLUP(A, b, x, maxResidual);

        // Record execution time
        auto end = std::chrono::high_resolution_clock::now();
        elapsed = end - start;

        if (success) {
            writeLUPResults(outFile, x, maxResidual);
        } else {
            outFile << "Error: LUP factorization failed!" << std::endl;
        }
    }
    else if (methodFlag == 1) {
        // SOR method
        int iterations;
        double residualNorm;
        bool converged = solveSOR(A, b, x, sorWeight, epsilon, maxIter, iterations, residualNorm);

        // Record execution time
        auto end = std::chrono::high_resolution_clock::now();
        elapsed = end - start;

        writeSORResults(outFile, x, iterations, residualNorm, sorWeight);

        if (!converged) {
            outFile << "Warning: SOR method did not converge within maximum iterations!" << std::endl;
        }
    }
    else if (methodFlag == 2) {
        // CG method
        int iterations;
        double residualNorm;
        bool converged = solveCG(A, b, x, epsilon, maxIter, iterations, residualNorm);

        // Record execution time
        auto end = std::chrono::high_resolution_clock::now();
        elapsed = end - start;

        writeCGResults(outFile, x, iterations, residualNorm);

        if (!converged) {
            outFile << "Warning: CG method did not converge within maximum iterations!" << std::endl;
        }
    }
    else if (methodFlag == 3) {
        // PCG method with Jacobi preconditioner
        int iterations;
        double residualNorm;
        bool converged = solveJacobiPCG(A, b, x, epsilon, maxIter, iterations, residualNorm);

        // Record execution time
        auto end = std::chrono::high_resolution_clock::now();
        elapsed = end - start;

        writePCGResults(outFile, x, iterations, residualNorm, A);

        if (!converged) {
            outFile << "Warning: PCG method did not converge within maximum iterations!" << std::endl;
        }
    }
    else if (methodFlag == 4) {
        // Power Iterations method
        int iterations;
        double error;
        bool converged = solvePowerIterations(A, b, epsilon, maxIter, x, iterations, error);

        // Record execution time
        auto end = std::chrono::high_resolution_clock::now();
        elapsed = end - start;

        writePowerIterationsResults(outFile, x, iterations, error);

        if (!converged) {
            outFile << "Warning: Power Iterations method did not converge within maximum iterations!" << std::endl;
        }
    }
    else {
        outFile << "Error: Invalid method flag!" << std::endl;
        return 1;
    }

    outFile << "Execution time (ms) = " << std::fixed << std::setprecision(8) << elapsed.count() * 1000.0 << std::endl;

    inFile.close();
    outFile.close();

    return 0;
}
