// newton_iteration.h
// Newton's Iteration solver for the nonlinear neutron diffusion equation
// April 13, 2025

#ifndef NEWTON_ITERATION_H
#define NEWTON_ITERATION_H

/**
 * Solves the nonlinear neutron diffusion equation using Newton's Iterations
 *
 * The equation solved is:
 * −(φ_i+1,j − 2φ_i,j + φ_i−1,j)/h² − (φ_i,j+1 − 2φ_i,j + φ_i,j−1)/h² + (ρ₀ + β/√φ_i,j)φ_i,j = 1
 *
 * This implementation constructs the Jacobian matrix, evaluates the nonlinear function vector,
 * and solves the resulting system using LU factorization to compute the change in flux.
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
void newtonIteration(double** flux, int n, double h, double rho0, double beta,
                     double epsilon, int maxIter, int& iterations, double& error);

/**
 * Maps 2D grid coordinates to 1D array index
 *
 * @param i Row index (1-based)
 * @param j Column index (1-based)
 * @param n Number of interior nodes per dimension
 * @return 1D array index (0-based)
 */
int flattenIndex(int i, int j, int n);

/**
 * Evaluates the nonlinear function F for a given flux
 *
 * @param F Output vector to store function values
 * @param flux 2D array of flux values
 * @param n Number of interior nodes per dimension
 * @param h Mesh size
 * @param rho0 Linear component of removal term
 * @param beta Nonlinear component of removal term
 */
void evaluateFunction(double* F, double** flux, int n, double h, double rho0, double beta);

/**
 * Constructs the Jacobian matrix for the Newton iteration
 *
 * @param J Output Jacobian matrix (stored as 1D array in row-major order)
 * @param flux 2D array of flux values
 * @param n Number of interior nodes per dimension
 * @param h Mesh size
 * @param rho0 Linear component of removal term
 * @param beta Nonlinear component of removal term
 */
void constructJacobian(double** J, double** flux, int n, double h, double rho0, double beta);

#endif // NEWTON_ITERATION_H