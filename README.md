# Inlab 6: Solve Matrix Equation with Point Jacobi Method

#### a) To compile the source code, open a terminal in the directory containing the source files and type the following commands:
```bash
mkdir build
cd build
cmake -S .. -B .
make
```

#### b) To run the program, use the command:
```bash
./inlab6
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
./inlab6
```
The code should now acknowledge the changed inputs and if you now type in:
```bash
cat output.txt
```
You should see the updated output file.

#### c) The code is: `operational`.

Here's a refactored README based on the code we've developed:

#### d) Pathname of the Source Codes:
`inlab6/main.cpp` -> Contains the main program module
`inlab6/LUPFactorization.cpp` -> Contains the main algorithm for LUP decomposition.
`inlab6/substitution.cpp` -> Implements forward and back substitution algorithms
`inlab6/ApplyPermutationMatrix.cpp` -> Handles permutation matrix operations
`inlab6/CalculateResiduals.cpp` -> Calculates the residual of Ax-b.
`inlab6/LUFactorization.cpp` -> Implements algorithm for LU decomposition in case pivoting is not used.
`inlab6/MatrixVectorProduct.cpp` -> Implementation of matrix vector product brought from Out-Lab 1.
`inlab6/input.txt` -> Contains n, usePivoting, A matrix and b vector.
`inlab6/build/output.txt` -> Contains A, L, U, P matrices and b vector, and solution vector x and the residuals of Ax-b.
`inlab6/Examples` -> Folder with sample input and output files

#### e) Brief Description of the Solved Problem
This program solves a system of linear equations Ax = b using an algorithm based on LUP Factorization.
Assuming the original matrix (A), the permutation matrix (P), and the right-hand side vector (b) are all given
in an input text file, the program performs LU factorization with pivoting and solves for the vector x in two stages:
- Forward Substitution to solve Ly = Pb
- Back Substitution to solve Ux = y

#### Variable Declarations

- `n`: The order of the square matrices (A, L, U, and P), i.e., the number of equations in the system.
- `A`: The original coefficient matrix A of size n x n.
- `L`: The unit lower triangular matrix L of size n x n.
- `U`: The upper triangular matrix U of size n x n.
- `P`: The permutation matrix P of size n x n.
- `b`: The right-hand side vector b of length n.
- `y`: The intermediate vector y obtained after forward substitution.
- `x`: The solution vector x of the system Ax = b.
- `usePivoting`: Boolean flag indicating whether to use pivoting in the factorization.

#### Function Declarations and Purposes

- `bool lupFactorize(A, L, U, P, usePivoting)`: Performs LUP factorization on matrix A, returns success status.
- `bool verifyLUPFactorization(A, L, U, P)`: Verifies if PA = LU and checks matrix properties.
- `void forwardSubstitution(L, Pb, y)`: Solves the system Ly = Pb using forward substitution.
- `void backSubstitution(U, y, x)`: Solves the system Ux = y using back substitution.
- `vector<double> calculateResidual(A, x, b)`: Calculates the residual vector Ax - b.

#### Step-by-Step Explanation of the Code

1. Task 1:
    - Writes a header to the output file `output.txt`, including the program title, author, affiliation, date, and a separator line.

2. Task 2:
    - Input Reading:
        - Opens `input.txt` for reading.
        - Reads the matrix order `n` and pivoting flag.
        - Reads the elements of matrix A.
        - Reads the elements of permutation matrix P.
        - Reads the right-hand side vector b.

3. Task 3:
    - LUP Factorization:
        - Performs LUP factorization using the lupFactorize function.
        - Verifies the factorization using verifyLUPFactorization.
        - If factorization fails, writes an error message and terminates.
        - If successful, proceeds with the solution.

4. Task 4:
    - Forward Substitution:
        - Applies permutation matrix P to vector b to get Pb.
        - Calls `forwardSubstitution(L, Pb, y)` to solve Ly = Pb.
    - Back Substitution:
        - Calls `backSubstitution(U, y, x)` to solve Ux = y.

5. Task 5:
    - Solution Verification:
        - Calculates residual vector using calculateResidual.
        - Computes maximum absolute residual.
    - Output Results:
        - Writes matrices A, L, U, and P to output file.
        - Writes solution vector x to output file.
        - Writes residual vector and maximum residual to output file.
