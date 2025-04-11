// fixed_point_iteration.h
// Fixed-Point Iteration solver for the nonlinear neutron diffusion equation
// April 11, 2025

#ifndef FIXED_POINT_ITERATION_H
#define FIXED_POINT_ITERATION_H

/**
 * Solves the nonlinear neutron diffusion equation using Fixed-Point Iterations
 *
 * The equation solved is:
 * −(φ_i+1,j − 2φ_i,j + φ_i−1,j)/h² − (φ_i,j+1 − 2φ_i,j + φ_i,j−1)/h² + (ρ₀ + β/√φ_i,j)φ_i,j = 1
 *
 * @param flux      2D array for flux values (including boundary nodes)
 * @param n         Number of interior nodes per dimension
 * @param h         Mesh size
 * @param rho0      Linear component of removal term
 * @param beta      Nonlinear component of removal term
 * @param epsilon   Stopping criterion for convergence
 * @param maxIter   Maximum number of iterations allowed
 * @param iterations Reference to store the actual number of iterations performed
 * @param error     Reference to store the final error achieved
 */
void fixedPointIteration(double** flux, int n, double h, double rho0, double beta,
                         double epsilon, int maxIter, int& iterations, double& error);

#endif // FIXED_POINT_ITERATION_H