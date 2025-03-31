# Outlab 11: Linear System Solvers with Power Iterations

## Compilation
```bash
make           # Build both main solver and performance test
```

## Execution
```bash
cd scripts
bsub < submit.sh    # Queues both `linear_solvers` and `performance_test`.
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
- `PI/power_iterations.h` - Power Iterations method declarations
- `PI/power_iterations.cpp` - Power Iterations method implementation

### Main Programs
- `main.cpp` - Main program for solving individual systems

## Problem Description
This program implements and compares five methods:
1. LUP (Lower-Upper decomposition with Pivoting) - direct method for linear systems
2. SOR (Successive Over-Relaxation) - iterative method for linear systems
3. CG (Conjugate Gradient) - iterative method for linear systems
4. PCG (Preconditioned Conjugate Gradient) - iterative method with Jacobi preconditioner for linear systems
5. Power Iterations - iterative method for computing the dominant eigenvector

### Input
- Input for main solver is read from `input.txt`
- Format for Power Iterations (flag 4):
  - Line 1: Flag (4 for Power Iterations) and unused parameter
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
  - Execution time

### Power Iterations Implementation
- Computes the dominant eigenvector of a matrix using the power method
- Uses the L-infinity norm (maximum absolute value) for vector normalization
- Convergence is determined by the L-infinity norm of the difference between consecutive normalized vectors
- Algorithm steps:
  1. Initialize with user-provided initial guess vector
  2. Normalize using the L-infinity norm
  3. Iterate: multiply by matrix A, normalize, check convergence
  4. Return normalized eigenvector and error metrics
- Input format (flag 4) uses the same structure as other solvers, but the right-hand side vector is interpreted as the initial guess
- Outputs the final eigenvector, iteration count, and convergence error

### Limitations
- Power Iterations computes only the dominant eigenvector and does not explicitly return the eigenvalue