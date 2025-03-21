# Inlab 10: Conjugate Gradient Method - Matrix Operations

## Compilation
```bash
make
```

## Execution
```bash
./cg_solver
```

## Status
Operational

## Source Code Files
- `matrix_operations.h` - Header file with matrix operation declarations
- `matrix_operations.cpp` - Implementation of matrix operations
- `main.cpp` - Main program file

## Problem Description
This program implements matrix operations required for the Conjugate Gradient method to solve systems of linear
equations Ax = b where A is a symmetric positive definite matrix.

### Input
- Input is read from `input.txt`
- Format:
    - Line 1: Stopping criterion epsilon, maximum iterations
    - Line 2: Matrix order n
    - Next n lines: Elements of the n×n matrix A
    - Last line: n elements of the right-hand-side vector b

### Output
- Output is written to `output.txt`
- Includes:
    - Header with problem description
    - Verification of matrix symmetry
    - Echo of input data
    - Execution time

### Operations Implemented
1. Scalar-vector multiplication
2. Vector addition
3. Scalar product of vectors
4. A-inner product
5. Matrix-vector product

### Limitations
- Matrix must be symmetric positive definite
- Program only checks for symmetry, not positive definiteness
- Full CG algorithm implementation is planned for Outlab 10