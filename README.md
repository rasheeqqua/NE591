# NE591
Repository for in-Lab and out-Lab assignments of NE:591 - Mathematical and Computational methods in NE.

#### a) To compile the code:
        1. Go inside the decompressed outlab1 folder and open the terminal in that directory.
        2. Type in the following commands in the terminal:
           ```mkdir build```
           ```cd build```
           ```cmake -S .. -B .```
           ```make```

#### b) To run the code:
        1. Run the following command in terminal:
           ```./outlab1```
        2. Then enter the values as you wish!

#### c) The code is `Operational`.

#### d) The code can be explained in the following steps:
        Step 1: The code asks for the constant value k. Matrix A is multiplied by k later on, and the resultant matrix
                is stored in matrix D. Then the code asks for 3 values: M, N and J. The dimension of matrices A and B
                are (M x N), meaning there are M number of rows and N number of columns. The dimension of the F matrix
                is (N x J).
        Step 2: If the user provided dimensions of the matrices are non-positive, then the code throws an error saying
                that the dimensions have to be bigger than 0.
        Step 3: According to the user provided dimensions, the code then resizes the already variables (A, B, F) that
                are used for storing the original (A, B, F) and resultant (C, D, E) matrices.
        Step 4: The original matrices (A, B, F) are then populated using a nested loop that has 2 loops in it. The outer
                loop increments the value `m/n` until the maximum row number is reached and the inner loop increments the
                value `n/j` until the maximum column number is reached.
                Matrix A -> The diagonal values have been assigned to 1.0 and the rest of the values are equal to 0.5
                Matrix B -> If the row number is bigger than the column number then the value is 0.75, otherwise 0.25
                Matrix F -> All the values are equal to { 1 / (row no. + column no.) } except the first one, since the
                first values is indexed at (0, 0), which will result in { 1 / (0 + 0) } or, 1 / 0 = infinity. So instead
                1 is added to the summation of n + j.
        Step 5: Matrix addition is performed simply by looping through the rows(m) of matrices A and B in the outer loop
                and through the columns (n) in the inner loop. The operation follows this simple formula:
                C(m, n) = A(m, n) + B(m, n)
        Step 6: Matrix A is multiplied by the constant value k using the formula:
                D(m, n) = A(m, n) * k;
        Step 7: Suppose the dimension of matrix A is (2 x 3) and the dimension of F is (3 x 4)
                In this case, the dimension of the resultant matrix, E = A x F, will be (2 x 4)
                and, the total number of multiplication operations will be 2 x 3 x 4 = 24
                So, we have a nested loop with 3 loops in it. The outer loop changes according to the row number of A.
                The second loop changes according to the column number of A or row number of F. Since, according to our
                example, the resultant matrix will have 8 elements in it due to its 2 x 4 dimension. So that means, if
                we loop through these first 2 loops, we will make 2 x 4 = 8 operations. So whatever summation we get in
                the innermost loop, we assign that value to the E matrix, exactly 8 times.
                The innermost loop performs the actual matrix multiplication using the following formula:
                E(m, j) = sum of all A(m, n) x F(n, j) where m, n, j = 1 to M, N, J respectively.
                According to our example, the outer loop runs 2 times, the second loop runs 3 times, and the inner loop
                runs 4 times, doing exactly (2 x 3 x 4) = 24 multiplication operations.
        Step 8: We loop through the rows and columns of each resultant matrices (C, D and E) and print the results in a
                presentable way.