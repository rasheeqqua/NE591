#include <iostream>
#include <fstream>
#include <vector>
#include <iomanip>
#include <chrono>
#include <cmath>
#include "matrix_operations.h"

bool isSymmetric(const std::vector<std::vector<double>>& A) {
    int n = A.size();
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            if (std::abs(A[i][j] - A[j][i]) > 1e-10) {
                return false;
            }
        }
    }
    return true;
}

int main() {
    auto start = std::chrono::high_resolution_clock::now();

    std::ifstream inFile("input.txt");
    std::ofstream outFile("output.txt");

    if (!inFile || !outFile) {
        std::cerr << "Error opening files!" << std::endl;
        return 1;
    }

    // Task 1: Write header to output file
    outFile << "NE 591 - Inlab 10 Code" << std::endl;
    outFile << "Implemented by Hasibul Hossain Rasheeq, March 21, 2025" << std::endl;
    outFile << "------------------------------------------------------" << std::endl << std::endl;
    outFile << "Solve Symmetric Positive Definite Matrix" << std::endl;
    outFile << "Equation with Conjugate Gradient Method" << std::endl << std::endl;

    // Task 2: Read input data
    double epsilon;
    int maxIter, n;

    inFile >> epsilon >> maxIter >> n;

    std::vector<std::vector<double>> A(n, std::vector<double>(n));
    std::vector<double> b(n);

    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            inFile >> A[i][j];
        }
    }

    for (int i = 0; i < n; i++) {
        inFile >> b[i];
    }

    // Task 3: Check correctness of input data
    bool symmetric = isSymmetric(A);

    if (symmetric) {
        outFile << "Matrix symmetry checked" << std::endl;
        outFile << "User must ensure it is positive definite" << std::endl;
        outFile << "-----------------------------------------" << std::endl << std::endl;

        // Echo input data to output file
        outFile << "stopping criterion on residual norm = " << std::scientific << std::setprecision(2) << epsilon << std::endl;
        outFile << "matrix is of order: " << n << std::endl;
        outFile << "Matrix A:" << std::endl;

        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                outFile << std::scientific << std::setprecision(2) << A[i][j] << " ";
            }
            outFile << std::endl;
        }

        outFile << std::endl << "RHS vector b:" << std::endl;
        for (int i = 0; i < n; i++) {
            outFile << std::scientific << std::setprecision(2) << b[i] << " ";
        }
        outFile << std::endl << std::endl;
    } else {
        outFile << "Error: Matrix is not symmetric!" << std::endl;
        return 1;
    }

    // Task 4: Record execution time
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = end - start;

    outFile << "Execution time (sec) = " << std::fixed << std::setprecision(4) << elapsed.count() << std::endl;

    inFile.close();
    outFile.close();

    return 0;
}