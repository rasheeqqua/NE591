// newton_iteration.cpp
// Implementation of Newton's Iteration solver for the nonlinear neutron diffusion equation
// April 13, 2025

#include "newton_iteration.h"
#include "../matrix/lu_factorization.h"
#include <cmath>
#include <algorithm>
#include <iostream>
#include <vector>

int flattenIndex(int i, int j, int n) {
    return (i - 1) * n + (j - 1);
}

void evaluateFunction(double* F, double** flux, int n, double h, double rho0, double beta) {
    double h2 = h * h;

    for (int i = 1; i <= n; i++) {
        for (int j = 1; j <= n; j++) {
            int idx = flattenIndex(i, j, n);

            // Calculate the diffusion term
            double diffusionTerm = -(flux[i+1][j] - 2.0 * flux[i][j] + flux[i-1][j]) / h2
                                  -(flux[i][j+1] - 2.0 * flux[i][j] + flux[i][j-1]) / h2;

            // Calculate the removal term
            double removalTerm = (rho0 + beta / sqrt(std::max(flux[i][j], 1e-10))) * flux[i][j];

            // Evaluate the function: diffusion + removal - source = 0
            F[idx] = diffusionTerm + removalTerm - 1.0;
        }
    }
}

void constructJacobian(double** J, double** flux, int n, double h, double rho0, double beta) {
    double h2 = h * h;
    int size = n * n;

    // Initialize Jacobian to zero
    for (int i = 0; i < size; i++) {
        for (int j = 0; j < size; j++) {
            J[i][j] = 0.0;
        }
    }

    // Fill the Jacobian matrix
    for (int i = 1; i <= n; i++) {
        for (int j = 1; j <= n; j++) {
            int row = flattenIndex(i, j, n);

            // Diagonal elements (partial derivative with respect to φ_i,j)
            double phi = std::max(flux[i][j], 1e-10);
            double sqrt_phi = sqrt(phi);

            // Derivative of the nonlinear removal term
            double dRemoval_dPhi = rho0 + beta / sqrt_phi - beta * phi / (2.0 * sqrt_phi * sqrt_phi * sqrt_phi);

            // Diagonal elements include diffusion and removal terms
            J[row][row] = 4.0 / h2 + dRemoval_dPhi;

            // Off-diagonal elements for neighboring nodes
            // North neighbor (i-1, j)
            if (i > 1) {
                int col = flattenIndex(i-1, j, n);
                J[row][col] = -1.0 / h2;
            }

            // South neighbor (i+1, j)
            if (i < n) {
                int col = flattenIndex(i+1, j, n);
                J[row][col] = -1.0 / h2;
            }

            // West neighbor (i, j-1)
            if (j > 1) {
                int col = flattenIndex(i, j-1, n);
                J[row][col] = -1.0 / h2;
            }

            // East neighbor (i, j+1)
            if (j < n) {
                int col = flattenIndex(i, j+1, n);
                J[row][col] = -1.0 / h2;
            }
        }
    }
}

void newtonIteration(double** flux, int n, double h, double rho0, double beta,
                     double epsilon, int maxIter, int& iterations, double& error) {
    int size = n * n;

    // Allocate memory for vectors and matrices
    double* F = new double[size];         // Function values
    double* deltaFlux = new double[size]; // Flux update

    // Allocate Jacobian matrix
    double** J = new double*[size];
    for (int i = 0; i < size; i++) {
        J[i] = new double[size];
    }

    // Main Newton iteration loop
    for (iterations = 0; iterations < maxIter; iterations++) {
        // Evaluate nonlinear function
        evaluateFunction(F, flux, n, h, rho0, beta);

        // Construct Jacobian matrix
        constructJacobian(J, flux, n, h, rho0, beta);

        // Solve the linear system J * deltaFlux = -F using LU factorization
        // First, negate F
        for (int i = 0; i < size; i++) {
            F[i] = -F[i];
        }

        // Solve the system
        luSolve(J, F, deltaFlux, size);

        // Update flux values and calculate max relative error
        error = 0.0;
        for (int i = 1; i <= n; i++) {
            for (int j = 1; j <= n; j++) {
                int idx = flattenIndex(i, j, n);
                double oldFlux = flux[i][j];
                double newFlux = oldFlux + deltaFlux[idx];

                // Ensure flux doesn't go negative
                newFlux = std::max(newFlux, 1e-10);

                // Calculate relative error
                double relError = std::abs(newFlux / oldFlux - 1.0);
                error = std::max(error, relError);

                // Update flux
                flux[i][j] = newFlux;
            }
        }

        // Check for convergence
        if (error < epsilon) {
            break;
        }
    }

    // Free allocated memory
    delete[] F;
    delete[] deltaFlux;
    for (int i = 0; i < size; i++) {
        delete[] J[i];
    }
    delete[] J;
}