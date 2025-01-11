Step 1: Getting the input from the user
    The code first prompts the user to provide 3 values:
      a) the angle in radians (whose sine value needs to be calculated)
      b) stopping criterion
      c) maximum series index

Step 2: Throwing errors for invalid inputs
    The program will crash and throw an error if:
      a) the absolute value of the angle is bigger than or equal to 1
      b) the stopping criterion is less than or equal to 0
      c) the index is less than or equal to 0

Step 3: Calculating the truncated Taylor series
    There are 2 nested loops for carrying out the calculation.
    i) The outer loop only increases the value of series index, N, 
    until the value of N equals to the provided maximum series index.
    ii) The first inner loop sums up the Taylor series up until N and stores the summation in current_value.
    iii) The second inner loop sums up the Taylor series up until N-1 and stores the summation in previous_value.

Step 4: Checking if the stopping criterion is met
    If the absolute difference between the current_value and the previous_value is less than or equal to the 
    stopping criterion value then we can say that the calculation has converged. The previous value will then be
    stored as the truncated_value since there is practically no noticeable difference between the previous_value
    and current_value. We are essentially saying that the iteration reached converged in the last index.
    We change the `converged` flag to `true`.

Step 5: Show message in case stopping criterion is not met
    If the iteration doesn't converge within the given series index, we let the user know that the calculation is
    not going to converge.

Step 6: Show the results in case of convergence
    In case the convergence is reached, we show the value of:
      a) the intrinsic function, sin(x), calculated using the `cmath` library.
      b) the truncated_value
      c) the absolute differences between the values of sin(x) and the truncated sin(x)