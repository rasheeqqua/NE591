# Outlab 12: Linear System Solvers with Eigenvalue Methods

## Compilation
```bash
make
```

## Execution
```bash
cd scripts
bsub < submit.sh
```

## Status
Operational

## Source Code Files
### Core Modules
- `matrix_modules/matrix_operations.h` - Matrix operation declarations
- `matrix_modules/matrix_operations.cpp` - Matrix operation implementations
- `matrix_modules/verify_positive_definite_matrix.h` - Matrix verification declarations
- `matrix_modules/verify_positive_definite_matrix.cpp` - Matrix verification implementations

### Solver Implementations
- `LUP/LUP_solver.h` - LUP decomposition method declarations
- `LUP/LUP_main.cpp` - LUP solver implementation
- `LUP/substitution.cpp` - Forward and back substitution implementations
- `PI/power_iterations.h` - Power Iterations method declarations
- `PI/power_iterations.cpp` - Power Iterations method implementation
- `II/inverse_iterations.h` - Inverse Iteration method declarations
- `II/inverse_iterations.cpp` - Inverse Iteration method implementation

### Main Programs
- `main.cpp` - Main program for solving individual systems

## Problem Description
This program implements and compares several methods:
1. LUP (Lower-Upper decomposition with Pivoting) - direct method for linear systems
2. SOR (Successive Over-Relaxation) - iterative method for linear systems
3. CG (Conjugate Gradient) - iterative method for linear systems
4. PCG (Preconditioned Conjugate Gradient) - iterative method with Jacobi preconditioner for linear systems
5. Power Iterations - iterative method for computing the dominant eigenvector and eigenvalue
6. Inverse Iteration - iterative method for computing an eigenvector close to a specified eigenvalue

### Input
- Input for main solver is read from `input.txt`
- Format for Inverse Iteration (flag 6):
  - Line 1: Flag (6 for Inverse Iteration) and approximate eigenvalue λ
  - Line 2: Stopping criterion epsilon, maximum iterations
  - Line 3: Matrix order n
  - Next n lines: Elements of the n×n matrix A
  - Last line: n elements of the initial guess vector

### Output
- Main solver output is written to `output.txt`
- Includes:
  - Header with problem description
  - Echo of input data
  - Solution information (iterations, residuals/errors)
  - For Power Iterations and Inverse Iteration:
    - Computed eigenvector
    - Computed eigenvalue
    - Residual vector and its infinity norm
  - Execution time

### Inverse Iteration Implementation
- Computes an eigenvector and eigenvalue near a specified shift λ
- Utilizes the LUP factorization from Lab 5 to factorize (A - λI)
- Reuses the same factorization for each iteration to solve linear systems efficiently
- Algorithm steps:
  1. Create shifted matrix (A - λI)
  2. Perform LUP factorization once at the beginning
  3. Initialize with user-provided initial guess vector
  4. Normalize using the L-infinity norm
  5. Iterate: solve (A - λI)y = x_{k-1} using the LUP factorization
  6. Normalize the result and compute eigenvalue using Rayleigh Quotient
  7. Check convergence on both eigenvector and eigenvalue
  8. Compute residual vector ρ = A·x - λ·x and its infinity norm
- Reports the eigenvalue of the original matrix A, not of the shifted matrix