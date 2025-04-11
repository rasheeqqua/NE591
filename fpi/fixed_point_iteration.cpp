
// fixed_point_iteration.cpp
// Implementation of Fixed-Point Iteration solver for the nonlinear neutron diffusion equation
// April 11, 2025

#include "fixed_point_iteration.h"
#include <cmath>
#include <algorithm>

void fixedPointIteration(double** flux, int n, double h, double rho0, double beta,
                         double epsilon, int maxIter, int& iterations, double& error) {

    // Create a temporary array to store the new iteration
    double** newFlux = new double*[n+2];
    for (int i = 0; i < n+2; i++) {
        newFlux[i] = new double[n+2];
        // Initialize with zeros (for boundary conditions)
        for (int j = 0; j < n+2; j++) {
            newFlux[i][j] = 0.0;
        }
    }

    // Precompute h squared for efficiency
    double h2 = h * h;
    double h2Inv = 1.0 / h2;

    // Main iteration loop
    for (iterations = 0; iterations < maxIter; iterations++) {
        // Reset maximum error for this iteration
        error = 0.0;

        // Update interior nodes using the fixed-point iteration formula
        for (int i = 1; i <= n; i++) {
            for (int j = 1; j <= n; j++) {
                // Sum of neighbor fluxes divided by h²
                double neighborSum = (flux[i+1][j] + flux[i-1][j] + flux[i][j+1] + flux[i][j-1]) * h2Inv;

                // Calculate nonlinear removal term
                double nonlinearRho = rho0 + beta / sqrt(std::max(flux[i][j], 1e-10)); // Avoid division by zero

                // Apply fixed-point iteration formula:
                // φ_i,j = [(φ_i+1,j + φ_i−1,j + φ_i,j+1 + φ_i,j−1)/h² + 1] / [4/h² + (ρ₀ + β/√φ_i,j)]
                newFlux[i][j] = (neighborSum + 1.0) / (4.0 * h2Inv + nonlinearRho);

                // Calculate relative error for this node
                if (flux[i][j] != 0.0) {
                    double relError = std::abs(newFlux[i][j] / flux[i][j] - 1.0);
                    error = std::max(error, relError);
                }
            }
        }

        // Check for convergence
        if (error < epsilon) {
            break;
        }

        // Copy new values to the flux array for next iteration
        for (int i = 1; i <= n; i++) {
            for (int j = 1; j <= n; j++) {
                flux[i][j] = newFlux[i][j];
            }
        }
    }

    // If we exited due to max iterations, make sure final values are copied
    if (iterations == maxIter) {
        for (int i = 1; i <= n; i++) {
            for (int j = 1; j <= n; j++) {
                flux[i][j] = newFlux[i][j];
            }
        }
    }

    // Free temporary array
    for (int i = 0; i < n+2; i++) {
        delete[] newFlux[i];
    }
    delete[] newFlux;
}