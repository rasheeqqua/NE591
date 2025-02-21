# Outlab 6: Solve Matrix Equation with Relaxation Methods

P.S.: To run the test, open the CMakeLists.txt file by typing in:
```bash
nano CMakeLists.txt
```
and then change line 6 from this:

```add_executable(outlab6 main.cpp)```

to this:

```add_executable(outlab6 test.cpp)```

and then press:
Ctrl + O, and then Enter to save
Ctrl + X to exit

#### a) Compilation Instructions:
```bash
mkdir build
cd build
cmake -S .. -B .
make
```

#### b) Execution Instructions:
```bash
./outlab6
cat output.txt
```

To modify input:
```bash
cd ..
nano input.txt
```

Then change the input values in the nano text editor. After you are done changing the values, press:
Ctrl + O, and then Enter to save
Ctrl + X to exit

Then type in terminal:
```bash
cd build
./outlab6
cat output.txt  # View updated results
```

#### c) Status: `operational`

#### d) Source Code Structure:
`outlab6/main.cpp` -> Contains program initialization, I/O handling, and method selection  
`outlab6/PointJacobi.cpp` -> Implements Point-Jacobi iteration method
`outlab6/GaussSeidel.cpp` -> Implements Gauss-Seidel method  
`outlab6/SOR.cpp` -> Implements for SOR method
`outlab6/CalculateResidual.cpp` -> Calculates residual norms and error metrics  
`outlab6/CheckDiagonalDominance.cpp` -> Verifies matrix diagonal dominance
`outlab6/input.txt` -> Contains method selection, convergence criteria, matrix A, and vector b  
`outlab6/build/output.txt` -> Contains iteration results, solution vector, and performance metrics
`outlab6/build/matrix_results.txt` -> Contains the results from the numerical experiments
`outlab6/Examples` -> Sample input/output files for testing
`outlab6/LUPFactorization` -> Contains the implementation of LUP factorization from outlab 5.

#### e) Problem Description
This program implements relaxation iterative methods for solving systems of linear equations Ax = b. Currently, it supports the Point-Jacobi method, with framework in place for Gauss-Seidel and SOR methods. The solution process involves iterative refinement until convergence criteria are met.

#### Variable Declarations
- `method`: Integer flag for method selection (0=Point-Jacobi, 1=Gauss-Seidel, 2=SOR)
- `omega`: Relaxation parameter for SOR method
- `epsilon`: Convergence criterion for relative error
- `maxIter`: Maximum allowed iterations
- `n`: Matrix order
- `A`: Coefficient matrix (n×n)
- `b`: Right-hand side vector
- `x`: Solution vector
- `actualIter`: Actual iterations performed
- `finalError`: Final achieved error

#### Function Declarations
- `bool isDiagonallyDominant(A)`: Checks matrix diagonal dominance
- `double calculateError(x_new, x_old)`: Computes infinity norm of relative change
- `double calculateResidual(A, x, b)`: Computes maximum absolute residual
- `bool pointJacobi(A, b, x, epsilon, maxIter, actualIter, finalError)`: Implements Point-Jacobi iteration
- - `bool gaussSeidel(A, b, x, epsilon, maxIter, actualIter, finalError)`: Implements Gauss-Seidel iteration
- - `bool sor(A, b, x, omega, epsilon, maxIter, actualIter, finalError)`: Implements SOR iteration

#### Implementation Steps

1. Input Processing and Validation:
   - Read method selection and parameters from `input.txt` (method, omega, epsilon, maxIter, n)
   - Validate input parameters (matrix size, maximum iterations, stopping criterion)
   - Read coefficient matrix A and vector b
   - Print input information to `output.txt`

2. Method Execution:
   - Initialize solution vector x with zeros
   - Based on method selection, perform the appropriate iterative method:
      - **Point-Jacobi**: Implement Point-Jacobi iteration
      - **Gauss-Seidel**: Implement Gauss-Seidel iteration
      - **SOR**: Implement Successive Over-Relaxation (SOR) iteration (validate omega ∈ (0,2))
   - For the selected method:
      - Perform iterative updates
      - Check convergence criteria (relative error and maximum iterations)
      - Track actual iterations and final error
      - Measure execution time for the iterative method

3. Direct Solution via LUP Factorization:
   - Perform LUP factorization on matrix A to obtain matrices L, U, and permutation matrix P
   - Solve the system using LUP factorization:
      - Apply permutation matrix to vector b to obtain Pb
      - Perform forward substitution to solve Ly = Pb
      - Perform back substitution to solve Ux = y
   - Obtain solution vector x_lup using LUP method
   - Measure execution time for LUP method

4. Result Comparison and Residual Calculation:
   - Compare the solution vectors obtained from the iterative method and LUP method
   - Calculate residuals for both solutions (maximum absolute residual)
   - Output solution vectors and residuals to `output.txt`
   - Report convergence status and performance metrics (iterations, final error, execution times)

5. Output Generation:
   - Write iteration results and convergence information to `output.txt`
   - Print solution vectors from both methods
   - Include residuals and execution times in the output