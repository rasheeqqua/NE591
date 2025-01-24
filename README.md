# In-Lab Assignment 3: Numerical Integration Program

#### a) To compile the source code, open a terminal in the directory containing the source files and type the following commands:
```bash
mkdir build
cd build
cmake -S .. -B .
make
```

#### b) To run the program, use the command:
```bash
./inlab3
```

#### c) The code is: `operational`.

#### d) Brief Description of the Solved Problem
This program implements numerical integration using composite quadrature rules to approximate the definite integral of a user-defined function over a specified interval.
The available methods are the Composite Trapezoidal Rule and the Composite Simpson's Rule.

#### Variable Declarations

- `a`: The lower limit of integration.
- `b`: The upper limit of integration.
- `m`: The number of equal intervals into which the integration interval [a, b] is divided.
- `selector`: The user's choice for the quadrature rule:
    - `1`: Composite Trapezoidal Rule
    - `2`: Composite Simpson's Rule
    - `3`: Gaussian Quadrature (not implemented)
- `result`: The computed value of the approximate integral.
- `inputValid`: Flag indicating whether the input data is valid.
- `h`: The width of each sub-interval { (b - a) / m }.

#### Function Declarations and Purposes

- `double integrand(double x)`:
Defined in `integrand.cpp` module, this function represents the integrand f(x) to be integrated over the interval [a, b].

- `double compositeTrapezoidal(double a, double b, int m)`:  
Computes the approximate integral using the Composite Trapezoidal Rule. It divides the interval [a, b] into `m` sub-intervals and
sums the areas of the trapezoids under the curve.

- `double compositeSimpson(double a, double b, int m)`:  
Computes the approximate integral using the Composite Simpson's Rule. It divides the interval [a, b] into `m` sub-intervals (requires `m` to be even) and
sums the areas using parabolic approximations.

#### Step-by-Step Explanation of the Code

1. Task 1: 
Prints program information to the terminal, including the program title, author, affiliation, date, and a brief description of its purpose.

2. Task 2:
Prompts the user to enter the lower and upper limits of integration (`a` and `b`), the number of intervals (`m`), and
the choice of quadrature rule (`selector`). Validates that `a < b`, `m > 0`, and `selector` is within the valid range.

3. Task 3: 
Checks the correctness of the input data. If incorrect data is detected, prints informative error messages to the terminal and terminates the program.
If all data is correct, congratulates the user and echoes the inputs in a professional style.

4. Task 4: 
The function f(x) to be integrated. This is defined in `integrand.cpp` module.
This allows users to modify the integrand without changing the main program.

5. Task 5: 
Based on the user's choice (`selector`), the program computes the approximate integral using the selected quadrature rule.
If Gaussian Quadrature is selected, prints an informative message and terminates.

- Composite Trapezoidal Rule:
      ```
      double compositeTrapezoidal(double a, double b, int m)
      // calculates: I = h/2 * [f(a) + f(a+h) + f(a+2h)... + f(b)]
      ```

- Composite Simpson's Rule:
      ```
      double compositeSimpson(double a, double b, int m)
      // calculates: I = h/3 * [f(a) + 4f(a+h) + 2f(a+2h) + 4f(a+3h) + 2f(a+4h)... + f(b)]
      ```

6. Task 6:  
Prints a summary of the results to the terminal, including the integration interval, the number of intervals,
the selected quadrature rule, and the computed value of the integral.

NOTE:
- Customizing the Integrand:  
To change the function being integrated, modify the `integrand(double x)` function in `integrand.cpp`.
For example, to integrate f(x) = x^2, you can change the function to:
  ```
  double integrand(double x) {
      return x * x;
  }
  ```