# In-Lab Assignment 5: LUP Factorization Solution Program

#### Notes:
- Customizing the Input:
    - The program reads input from a file named `input.txt`.
    - The input file must be formatted as specified in the `input.txt` file.

- Viewing the results:
    - All the results are printed to `output.txt` file inside the `build` folder. Open the `output.txt` file to see the messages and results.

#### a) To compile the source code, open a terminal in the directory containing the source files and type the following commands:
```bash
mkdir build
cd build
cmake -S .. -B .
make

```

#### b) To run the program, use the command:
```bash
./inlab4
```

#### c) The code is: `operational`.

#### d) Pathname of the source codes:
`inlab4/main.cpp` -> contains the main module excluding the forward and back substitution algorithms.
`inlab4/substitution.cpp` -> Contains the algorithms for forward and back substitution algorithms.
`inlab4/input.txt` -> Contains the order of the matrices, the elements of L, U, and b matrices.
`inlab4/InputOutput` -> This folder contains sample input and output files.

#### e) Brief Description of the Solved Problem
This program solves a system of linear equations Ax = b using an algorithm based on LU Factorization.
Assuming the lower triangular matrix (L), the upper triangular matrix (U) and the right-hand side vector (b) are all given
in an input text file, the program solves for the vector x in two stages:
- Forward Substitution to solve Ly = b
- Back Substitution to solve Ux = y

#### Variable Declarations

- `n`: The order of the square matrices (L and U), i.e., the number of equations in the system.
- `L`: The unit lower triangular matrix L of size  n x n.
- `U`: The upper triangular matrix U of size n x n.
- `b`: The right-hand side vector b of length n.
- `y`: The intermediate vector y obtained after forward substitution.
- `x`: The solution vector x of the system  Ax = b.
- `valid`: Flag indicating whether the input data is valid (e.g., no zero diagonal elements in U).

#### Function Declarations and Purposes

- `void forwardSubstitution(L, b, y)`: Solves the system Ly = b using forward substitution.
- `void backSubstitution(U, y, x)`:Solves the system Ux = y using back substitution.
- `bool checkInputData(n, U, outputFile)`: Checks the correctness of the input data, specifically that the diagonal elements of  U  are non-zero.

#### Step-by-Step Explanation of the Code

1. Task 1:
    - Writes a header to the output file `output.txt`, including the program title, author, affiliation, date, and a separator line.

2. Task 2:
    - Input Reading:
        - Opens `input.txt` for reading.
        - Reads the matrix order `n`.
        - Initializes the matrices `L` and `U` as  nxn  with all the elements set to zero.
        - Sets the diagonal elements of `L` to 1.
        - Reads the non-zero elements of `L` below the diagonal, row by row.
        - Reads the non-zero elements of `U` on and above the diagonal, row by row.
        - Reads the `n` elements of the right-hand side vector `b`.
    - Input Format:
        - For Matrix  L :
            - For row `i` from 1 to `n - 1`, read `i` elements for columns `0` to `i - 1`.
        - For Matrix  U :
            - For row `i` from 0 to `n - 1`, read `n - i` elements for columns `i` to `n - 1`.
        - For Vector `b`:
            - Read `n` elements.

3. Task 3:
    - Input Validation:
        - Checks that `n` is positive.
        - Uses `checkInputData` function to ensure diagonal elements of `U` are non-zero.
        - If invalid data is detected, writes an informative error message to the output file and terminates.
        - If all data is correct, congratulates the user and proceeds.
    - Echoes Input Data:
        - Writes the matrix order, and the matrices `L`, `U`, and vector `b` to the output file in a professional format.

4. Task 4:
    - Forward Substitution:
        - Calls `forwardSubstitution(L, b, y)` to solve  Ly = b  and obtains the intermediate vector `y`.
    - Back Substitution:
        - Calls `backSubstitution(U, y, x)` to solve  Ux = y  and obtains the solution vector `x`.

5. Task 5:
    - Output Solution:
        - Writes the solution vector `x` to the output file `output.txt` with appropriate formatting.