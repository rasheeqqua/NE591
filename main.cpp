// main.cpp for NE 591 Inlab 13
// Nonlinear Neutron Diffusion Equation solver using Fixed-Point Iterations
// April 11, 2025

#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <chrono>
#include "fpi/fixed_point_iteration.h"

int main() {
    // Input parameters
    double epsilon;       // Stopping criterion
    int maxIterations;    // Maximum number of iterations
    int n;                // Number of nodes per dimension (excluding boundary)
    double rho0, beta;    // Linear and nonlinear components of removal term
    double L;             // Domain side length

    // Calculation variables
    double h;             // Mesh size
    double** flux;        // 2D array for flux values
    int iterations = 0;   // Actual iterations used
    double error = 0.0;   // Final error achieved
    double avgFlux = 0.0; // Average flux over domain

    // Read input file
    std::ifstream inputFile("input.txt");
    if (!inputFile.is_open()) {
        std::cerr << "Error: Cannot open input file!" << std::endl;
        return 1;
    }

    // Read parameters
    inputFile >> epsilon >> maxIterations;
    inputFile >> n;
    inputFile >> rho0 >> beta;
    inputFile >> L;
    inputFile.close();

    // Verify input parameters
    if (epsilon <= 0 || maxIterations <= 0 || n <= 0 || L <= 0) {
        std::cerr << "Error: Invalid input parameters!" << std::endl;
        return 1;
    }

    // Calculate mesh size
    h = L / (n + 1);

    // Allocate memory for flux array (including boundary nodes)
    flux = new double*[n+2];
    for (int i = 0; i < n+2; i++) {
        flux[i] = new double[n+2];
        // Initialize with zeros (including boundary)
        for (int j = 0; j < n+2; j++) {
            flux[i][j] = 0.0;
        }
    }

    // Initialize interior nodes with initial guess = 1.0
    for (int i = 1; i <= n; i++) {
        for (int j = 1; j <= n; j++) {
            flux[i][j] = 1.0;
        }
    }

    // Start timer
    auto startTime = std::chrono::high_resolution_clock::now();

    // Call fixed-point iteration solver
    fixedPointIteration(flux, n, h, rho0, beta, epsilon, maxIterations, iterations, error);

    // End timer
    auto endTime = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(endTime - startTime);
    double execTime = duration.count() / 1000.0; // Convert to milliseconds

    // Calculate average flux
    double sum = 0.0;
    for (int i = 1; i <= n; i++) {
        for (int j = 1; j <= n; j++) {
            sum += flux[i][j];
        }
    }
    avgFlux = sum / (n * n);

    // Write output file
    std::ofstream outputFile("output.txt");
    if (!outputFile.is_open()) {
        std::cerr << "Error: Cannot open output file!" << std::endl;
        return 1;
    }

    // Header information
    outputFile << "NE 591 Inlab 13 Code" << std::endl;
    outputFile << "Implemented by Hasibul H. Rasheeq, April 11, 2025" << std::endl;
    outputFile << "-------------------------------------------------" << std::endl << std::endl;

    outputFile << "Solve 2D Neutron Diffusion Nonlinear Equation with Fixed-Point Iterations" << std::endl;
    outputFile << "-------------------------------------------------------------------------" << std::endl << std::endl;

    // Echo input parameters
    outputFile << "Number of nodes/dimension = " << n << std::endl;
    outputFile << "Linear removal term in diffusion equation = " << std::scientific << std::setprecision(2) << rho0 << std::endl;
    outputFile << "Nonlinear removal term = " << std::scientific << std::setprecision(2) << beta << std::endl;
    outputFile << "Domain-side length = " << std::scientific << std::setprecision(2) << L << std::endl;
    outputFile << "Maximum number of iterations = " << maxIterations << std::endl;
    outputFile << "Stopping criterion = " << std::scientific << std::setprecision(2) << epsilon << std::endl;
    outputFile << "-----------------------------------------------------" << std::endl << std::endl;

    // Results of iterations
    outputFile << "Iterations converged in = " << iterations << " iterations" << std::endl;
    outputFile << "Iterative error = " << std::scientific << std::setprecision(4) << error << std::endl << std::endl;

    outputFile << "The scalar flux obtained from the last iterate: in file `Flux`" << std::endl;
    outputFile << "--------------------------------------------------------------" << std::endl << std::endl;

    // Summary
    outputFile << "Summary of calculation:" << std::endl;
    outputFile << "nodes, removal, beta = " << n << ", " << std::scientific << std::setprecision(2) << rho0 << " " << beta << std::endl;
    outputFile << "Domain-side length = " << std::scientific << std::setprecision(2) << L << std::endl;
    outputFile << "Iters, Average flux = " << iterations << ", " << std::scientific << std::setprecision(4) << avgFlux << std::endl << std::endl;

    // Execution time
    outputFile << "Execution time (ms) = " << std::fixed << std::setprecision(8) << execTime << std::endl;

    outputFile.close();

    // Write flux file
    std::ofstream fluxFile("Flux");
    if (!fluxFile.is_open()) {
        std::cerr << "Error: Cannot open Flux file!" << std::endl;
        return 1;
    }

    // Write flux values in specified format
    for (int i = 1; i <= n; i++) {
        for (int j = 1; j <= n; j++) {
            fluxFile << i << " " << j << " " << std::scientific << std::setprecision(6) << flux[i][j] << std::endl;
        }
    }

    fluxFile.close();

    // Free allocated memory
    for (int i = 0; i < n+2; i++) {
        delete[] flux[i];
    }
    delete[] flux;

    return 0;
}