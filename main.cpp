#include <cmath>
#include <iomanip>
#include <iostream>
#include <stdexcept>

/*--------------------Variable Declaration--------------------*/
/*
 * angle - Stores the radian angle provided by the user.
 * tolerance - Stores the stopping criterion provided by the user.
 * index - Stores the maximum series index provided by the user.
 * converged - Boolean flag. Indicates the solution has converged, if set to TRUE.
 * truncated_value - Stores the truncated Taylor series value, if the solution converges.
 * intrinsic_value - Stores the sin(x) value calculated using C++ library.
 * current_value - Stores the summation of terms calculated up to N.
 * previous_value - Stores the summation of terms calculated up to N-1.
 * int factorial(const int n) - Function for calculating the factorial of a real number, n.
 */
/*------------------End Variable Declaration------------------*/

int factorial(const int n) {
  int value = 1;
  for(int i = 1; i <= 2*n+1; i++) {
    value *= i;
  }
  return value;
}

int main() {
  double angle;
  double tolerance;
  int index;

  bool converged = false;
  double truncated_value = 0;

  // Step 1: Getting the input from the user
  std::cout << "Enter angle (in radians) whose SINE you wish to evaluate: ";
  std::cin >> angle;
  const double intrinsic_value = sin(angle);

  std:: cout << "Enter the stopping criterion: ";
  std::cin >> tolerance;

  std::cout << "Enter the maximum value of the series index: ";
  std::cin >> index;

  // Step 2: Throwing errors for invalid inputs
  if(std::abs(angle) >= 1) {
    throw std::invalid_argument("The absolute value of the angle has to be less than 1");
  } else if(tolerance <= 0) {
    throw std::invalid_argument("The stopping criterion has to be bigger than 0");
  } else if(index <= 0) {
    throw std::invalid_argument("The index has to be bigger than 0");
  } else {
    std::cout << "All the provided inputs are correct! Calculating the Taylor series now, give me a few moments!" << std::endl;
  }

  // Step 3: Calculating the truncated Taylor series
  for(int N = 1; N <= index; N++) {
    double current_value = 0;
    double previous_value = 0;

    for (int k = 0; k <= N; k++) {
      current_value += pow(-1, k) * (pow(angle, 2*k+1) / factorial(k));
    }

    for (int k = 0; k < N; k++)
    {
      previous_value += pow(-1, k) * (pow(angle, 2*k+1) / factorial(k));
    }

    // Step 4: Checking if the stopping criterion is met
    if(std::abs(current_value - previous_value) <= tolerance) {
      converged = true;
      std::cout << "Taylor series converged successfully at series index: " << N << std::endl;
      truncated_value = current_value;
      break;
    }
  }

  // Step 5: Show message in case stopping criterion is not met
  if(converged == false) {
    std::cout << "Taylor series did not converge successfully at series index: " << index << std::endl;
  } else {
    // Step 6: Show the results in case of convergence
    std::cout << "Results: " << std::endl;
    std::cout << "sin(x) computed with intrinsic function: " << intrinsic_value << std::endl;
    std::cout << "sin(x) computed with truncated Taylor series: " << std::setprecision (8) << truncated_value << std::endl;
    std::cout << "The difference between these two values: " << std::scientific << std::abs(intrinsic_value - truncated_value) << std::endl;
  }

  return 0;
}