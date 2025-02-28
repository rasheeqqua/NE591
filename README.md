# Project Milestone 2: Steady State One-Speed Diffusion Equation Solver

#### a) To compile the source code, open a terminal in the directory containing the source files and type the following commands:
```bash
chmod +x compile.sh
./compile.sh
```

#### b) To queue the program, use the command:
```bash
bsub < submit.sh
```

THIS WILL RUN ALL 3 PROGRAMS: i) Diffusion Solver program, ii) Verification program (Rotational and Reflective
verification), and iii) Performance program

To read only the output obtained from the diffusion solver program, type in:
```bash
cat output.txt
```

To read only the output from the verification program, type in:
```bash
cat reflective_*.txt
cat rotational_*.txt
```

To read only the output from the performance program, type in:
```bash
cat performance_*.txt
```

The input files for the verification test and performance test are generated automatically. If you want to
edit the input file for the Diffusion Solver program, then type in:

```bash
nano input.txt
```

Then change the input values in the `nano` Text Editor's window.
After you are done changing the values, then press `ctrl + o` to save the input file.
This will prompt you to ensure the filename. Just press `enter`.
After this press `ctrl + x` to exit the nano editor.
And then type in:

```bash
bsub < submit.sh
```

to queue the jobs again.

If you want to remove the output files, then type in:
```bash
chmod +x remove_output.sh
./remove_output.sh
```


#### c) The code is: `operational`.

#### d) Pathname of the Source Codes:
`DiffusionSolver.cpp` -> Contains the DiffusionSolver class implementation with all solver methods
`main.cpp` -> Contains the main program module
`PerformanceAnalyzer.cpp` -> Contains code for evaluating performance of different solution methods
`VerificationTest.cpp` -> Contains code for verification tests with different boundary conditions
`LUP/LUPFactorization.cpp` -> Contains the main algorithm for LUP decomposition
`LUP/substitution.cpp` -> Implements forward and back substitution algorithms
`LUP/ApplyPermutationMatrix.cpp` -> Handles permutation matrix operations
`LUP/MatrixVectorProduct.cpp` -> Implementation of matrix vector product
`Examples/General Input Output/` -> Contains the input and output files for the diffusion solver
`Examples/Performance Results/` -> Contains performance results for the different solver methods
`Examples/Verification Results/` -> Contains verification test results
`compile.sh` -> Script to compile all executables
`submit.sh` -> Script to run all executables on the Hazel cluster
`remove_output.sh` -> Script to clean up output files

#### e) Brief Description of the Solved Problem
This program solves the steady-state, one-speed diffusion equation in a 2D rectangular region with vacuum boundary conditions. The program implements and compares four solution methods:
- Direct solution using LUP Decomposition
- Point Jacobi iterative method
- Gauss-Seidel iterative method
- Successive Over-Relaxation (SOR) iterative method

The diffusion equation solved is:
- -D∇²φ(x,y) + Σₐφ(x,y) = q(x,y)

where φ is the neutron flux, D is the diffusion coefficient, Σₐ is the macroscopic absorption cross-section, and q is the source term.

#### Variable Declarations

- `flag`: Integer flag indicating which solution method to use (0: LUP, 1: Jacobi, 2: Gauss-Seidel, 3: SOR)
- `a, b`: Rectangle dimensions in cm
- `m, n`: Grid dimensions (number of interior points)
- `D`: Diffusion coefficient (cm)
- `sigma_a`: Macroscopic removal cross section (cm⁻¹)
- `q`: 2D vector storing source term values
- `delta, gamma`: Grid spacing in x and y directions
- `maxIterations`: Maximum number of iterations allowed for iterative methods
- `tolerance`: Convergence tolerance for iterative methods
- `omega`: Relaxation parameter for SOR method
- `phi`: 2D vector storing the scalar neutron flux solution
- `iterations`: Number of iterations performed
- `finalError`: Final error after iterations
- `converged`: Boolean indicating whether the solution converged
- `executionTime`: Time taken to execute the solution

#### Function Declarations and Purposes

- `bool readInput(string)`: Reads input parameters from file
- `void setParameters(...)`: Manually sets parameters (for testing)
- `vector<vector<double>> solve(int&, double&, bool&)`: Main solve function that dispatches to appropriate method
- `bool solvePointJacobi(vector<vector<double>>&, int&, double&)`: Solves using Point Jacobi method
- `bool solveGaussSeidel(vector<vector<double>>&, int&, double&)`: Solves using Gauss-Seidel method
- `bool solveSOR(vector<vector<double>>&, int&, double&)`: Solves using SOR method
- `vector<vector<double>> solveLUP()`: Solves using LUP direct method
- `void writeOutput(...)`: Writes solution and statistics to output file
- `double calculateError(const vector<vector<double>>&, const vector<vector<double>>&)`: Calculates error between iterations
- `double calculateMaxResidual(const vector<vector<double>>&)`: Calculates maximum residual
- `double compareSolutions(const vector<vector<double>>&, const vector<vector<double>>&)`: Compares two solutions
- `void applyDiffusionOperator(const vector<vector<double>>&, vector<vector<double>>&)`: Matrix-free operator application

#### Step-by-Step Explanation of the Code

1. Input Processing:
    - The program reads an input file containing the solution method flag, iteration parameters, problem geometry, grid dimensions, physical parameters, and source term values.
    - The input validation ensures all parameters are within valid ranges.

2. Solution Method Selection:
    - Based on the flag value, the program dispatches to the appropriate solution method:
        - Flag 0: Direct solution using LUP decomposition
        - Flag 1: Point Jacobi iterative method
        - Flag 2: Gauss-Seidel iterative method
        - Flag 3: Successive Over-Relaxation (SOR) method

3. Iterative Solution (for flags 1-3):
    - The solution starts with an initial guess of zero (vacuum boundary conditions).
    - For each iteration:
        - The program updates flux values based on the specific method's formula.
        - Error is calculated between consecutive iterations.
        - If the error falls below the tolerance, the solution is considered converged.
        - If max iterations are reached without convergence, the solution terminates.

4. Direct Solution (for flag 0):
    - The system is solved using LUP decomposition from Milestone 1.
    - The coefficient matrix and RHS vector are constructed from the diffusion equation.
    - LUP factorization is performed, followed by forward and back substitution.

5. Performance Analysis:
    - The program measures execution time.
    - For iterative methods, it tracks the number of iterations and final error.
    - It calculates the maximum residual to verify solution accuracy.

6. Output Generation:
    - Results are written to an output file, including:
        - Solution method and parameters
        - Problem specifications
        - Iteration information (for iterative methods)
        - Performance metrics
        - The scalar flux solution

7. Verification Testing (optional):
    - For iterative methods with small grids, the program can compare the solution with the LUP direct method for verification.
    - Maximum relative difference between solutions is calculated.