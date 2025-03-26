# Outlab 10: Linear System Solvers Comparison

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
- `LUP/LUP_solver.h` - LUP decomposition method declarations
- `LUP/LUP_main.cpp` - LUP decomposition implementation
- `LUP/LUPFactorization.cpp` - LUP factorization routines
- `LUP/ApplyPermutationMatrix.cpp` - Permutation application routines
- `LUP/substitution.cpp` - Forward and back substitution routines
- `SOR/SOR_solver.h` - SOR method declarations
- `SOR/SOR_main.cpp` - SOR method implementation
- `CG/CG_solver.h` - CG method declarations
- `CG/CG_main.cpp` - CG method implementation

### Main Programs
- `main.cpp` - Main program for solving individual systems
- `performance_test.cpp` - Program for performance comparison of methods

## Problem Description
This program implements and compares three methods for solving linear systems Ax = b where A is symmetric positive definite:
1. LUP (Lower-Upper decomposition with Pivoting) - direct method
2. SOR (Successive Over-Relaxation) - iterative method
3. CG (Conjugate Gradient) - iterative method

### Input
- Input for main solver is read from `input.txt`
- Format:
  - Line 1: Flag (0 for LUP, 1 for SOR, 2 for CG) and SOR weight parameter
  - Line 2: Stopping criterion epsilon, maximum iterations
  - Line 3: Matrix order n
  - Next n lines: Elements of the n×n matrix A
  - Last line: n elements of the right-hand-side vector b

### Output
- Main solver output is written to `output.txt`
- Performance test results are written to `performance_results.txt`
- Input matrices for performance test are written to `input_n*.txt` files
- Includes:
  - Header with problem description
  - Verification of matrix properties
  - Echo of input data
  - Solution information (iterations, residuals)
  - Execution time

### Performance Testing
Performance test implements:
- Testing with matrix sizes n = 32, 64, 128, 512, 1024
- Generation of SPD matrices
- Comparison of iteration counts for iterative methods
- Comparison of execution times across all methods

### Limitations
- Matrix must be symmetric positive definite
- SOR convergence depends on weight parameter