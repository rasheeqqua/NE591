# Out-Lab Assignment 2

#### a) To compile the source code, go into the decompressed `outlab2` folder and open a terminal there. Then type in the following terminal commands:
```mkdir build```
```cd build```
```cmake -S .. -B .```
```make```

#### b) To run the code type this command in the terminal:
```./outlab2```

#### c) The code is: `operational`.

#### d) Brief Description of the Solved Problem

#### Lagrange Interpolation Program
This project implements the Lagrange Interpolation Polynomial to approximate the function f(x) = e^x over the interval [-1, 1].

#### Variable Declaration
- `n`: Number of interpolation points provided by the user.
- `m`: Number of evaluation points where the interpolation polynomial is evaluated.
- `userX`: Vector containing the user-provided interpolation points \( x_i \).
- `x`: Vector storing interpolation points used for computation (may come from file or `userX`).
- `y`: Vector storing corresponding function values \( y_i \).
- `choice`: User's choice for data input method (1: from file, 2: evaluate \( f(x) \)).
- `a`, `b`: Start and end points of the interpolation interval.
- `xEval`: Vector containing evaluation points generated between `a` and `b`.
- `EPSILON`: Small constant used for floating-point comparison to check equality of values.
- `wIndex`, `wX`, `wY`, `wInterp`, `wError`: Variables used to set column widths for formatted output.

#### Function Declarations and Purposes
- `double lagrangeInterpolate(double xEval, const std::vector<double>& x, const std::vector<double>& y)`:  
  Computes the interpolated value at `xEval` using the Lagrange interpolation formula with given data points `x` and `y`.

- `std::vector<double> evaluate(const std::vector<double>& input)`:  
  Evaluates the function \( f(x) = e^x \) at each point in the `input` vector and returns the resulting values.

- `std::vector<double> createEvaluationPoints(int m, double a, double b)`:  
  Generates `m` equally spaced evaluation points between `a` and `b`.

#### Step by step explanation of the code

1. Program Introduction:  
   Displays program information, including the assignment title, author, date, and purpose.

2. User Input for `n` and `m`:  
   Prompts the user to input the number of interpolation points (`n`) and evaluation points (`m`). Validates that both are positive integers.

3. Input Interpolation Points `x_i`:
    - Prompts the user to input `n` interpolation points in ascending order.
    - Stores these points in `userX` and sorts them to ensure correct ordering.

4. Choose Data Input Method:
    - Option 1: Read `(x, y)` pairs from a text file.
    - Option 2: Use the provided `x_i` values and compute `y_i` by evaluating f(x) = e^x.
    - Stores the user's choice in `choice`.

5. Data Preparation:
    - If Option 1:
        - Reads data from a specified file.
        - Parses and stores the `x` and `y` values in their respective vectors.
        - Determines the interpolation interval `[a, b]` from the data.
    - If Option 2:
        - Sets `x` to `userX`.
        - Computes `y` by evaluating f(x) for each `x_i`.
        - Determines the interpolation interval `[a, b]` from `x`.

6. Confirmation and Proceeding:  
   Displays a summary of the inputs and confirms proceeding with the interpolation process.

7. Generate Evaluation Points:  
   Calls `createEvaluationPoints` to generate `m` equally spaced points between `a` and `b`, storing them in `xEval`.

8. Interpolation and Error Calculation:
    - For each evaluation point in `xEval`:
        - Computes the interpolated value using `lagrangeInterpolate`.
        - If Option 2 (Function Known):
            - Computes the actual function value at `xEval[i]`.
            - Calculates the error for all evaluation points as the difference between the actual and interpolated values.
        - If Option 1 (Data from File):
            - Checks if `xEval[i]` matches any `x_i` in the data.
            - If it matches, calculates the error; otherwise, indicates that the error is not computed.

9. Display Results:  
   Formats and displays the results in a table, including the index, `x(i)`, actual function or file value, interpolated value, and error.

10. Completion Message:  
    Displays a message indicating that the program has finished executing.
