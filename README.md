# Outlab 6: Solve Matrix Equation with Relaxation Methods

#### a) To compile the source code, open a terminal in the directory containing the source files and type the following commands:
```bash
mkdir build
cd build
cmake -S .. -B .
make
```

#### b) To run the program, use the command:
```bash
./outlab6
cat output.txt
```
If you want to update the input file then use the following commands:
```bash
cd ..
nano input.txt
```
Then change the input values in the `nano` Text Editor's window.
After you are done changing the values, then press `ctrl + o` to save the input file.
This will prompt you to ensure the filename. Just press `enter`.
After this press `ctrl + x` to exit the nano editor.

Now type in the following commands:
```bash
cd build
./outlab6
```
The code should now acknowledge the changed inputs and if you now type in:
```bash
cat output.txt
```
You should see the updated output file.

#### c) The code is: `operational`.

# Outlab 6: Solve Matrix Equation with Relaxation Methods

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
# Make changes in nano
# Ctrl + O, Enter to save
# Ctrl + X to exit
cd build
./outlab6
cat output.txt  # View updated results
```

#### c) Status: `operational`

#### d) Source Code Structure:
`outlab6/main.cpp` -> Contains program initialization, I/O handling, and method selection  
`outlab6/PointJacobi.cpp` -> Implements Point-Jacobi iteration method  
`outlab6/CalculateResidual.cpp` -> Calculates residual norms and error metrics  
`outlab6/CheckDiagonalDominance.cpp` -> Verifies matrix diagonal dominance  
`outlab6/GaussSeidel.cpp` -> Placeholder for Gauss-Seidel method (future implementation)  
`outlab6/SOR.cpp` -> Placeholder for SOR method (future implementation)  
`outlab6/input.txt` -> Contains method selection, convergence criteria, matrix A, and vector b  
`outlab6/build/output.txt` -> Contains iteration results, solution vector, and performance metrics  
`outlab6/Examples` -> Sample input/output files for testing

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
- Placeholder functions for Gauss-Seidel and SOR methods

#### Implementation Steps
1. Input Processing:
   - Read method selection and parameters
   - Read matrix A and vector b
   - Validate input data

2. Method Selection:
   - Execute appropriate iterative method
   - Track convergence and iterations

3. Solution Computation:
   - Perform iterative updates
   - Check convergence criteria
   - Calculate residuals

4. Output Generation:
   - Write iteration results
   - Report solution vector
   - Include performance metrics