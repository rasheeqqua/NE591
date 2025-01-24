/*--------------------In-Lab Assignment 3--------------------*/
/*
 * Created by Hasibul Hossain Rasheeq on Sunday, 1/24/25.
*/

#include <iostream>
#include <ctime>
#include "CompositeTrapezoidal.cpp"
#include "CompositeSimpson.cpp"

int main() {
    std::cout << "Numerical Integration Program\n";
    std::cout << "Author: Hasibul H. Rasheeq, NC State University\n";

    // Task 1
    time_t now = time(0);
    char* dt = ctime(&now);
    std::cout << "Date: " << dt << "\n";

    std::cout << "This program computes the approximate integral using numerical quadrature rules\n";
    std::cout << "-----------------------------------------------------------------------------------\n\n";

    // Task 2
    double a, b;
    int m, selector;

    std::cout << "Enter the lower limit of integration (a): ";
    std::cin >> a;

    std::cout << "Enter the upper limit of integration (b): ";
    std::cin >> b;

    std::cout << "Enter the number of intervals (m): ";
    std::cin >> m;

    std::cout << "Select the Quadrature Rule:\n";
    std::cout << "1. Composite Trapezoidal Rule\n";
    std::cout << "2. Composite Simpson's Rule\n";
    std::cout << "3. Gaussian Quadrature\n";
    std::cout << "Enter your choice (1-3): ";
    std::cin >> selector;

    // Task 3
    bool inputValid = true;

    if (a >= b) {
        std::cout << "Error: The lower limit 'a' must be less than the upper limit 'b'.\n";
        inputValid = false;
    }

    if (m <= 0) {
        std::cout << "Error: The number of intervals 'm' must be a positive integer.\n";
        inputValid = false;
    }

    if (selector < 1 || selector > 3) {
        std::cout << "Error: Invalid selection for the Quadrature Rule.\n";
        inputValid = false;
    }

    if (!inputValid) {
        std::cout << "Please correct the input data and try again.\n";
        return 1; // Terminate the program due to invalid input
    }

    std::cout << "All input data is correct. Proceeding with the computation.\n\n";
    std::cout << "Integration Interval: [" << a << ", " << b << "]\n";
    std::cout << "Number of Intervals (m): " << m << "\n";
    std::cout << "Selected Quadrature Rule: ";
    if (selector == 1)
        std::cout << "Composite Trapezoidal Rule\n";
    else if (selector == 2)
        std::cout << "Composite Simpson's Rule\n";
    else
        std::cout << "Gaussian Quadrature\n";
    std::cout << "----------------------------------------------------------\n\n";

    // Task 4 implementation is in `integrand.cpp` file.

    // Task 5
    double result = 0.0;

    switch (selector) {
    case 1:
        result = compositeTrapezoidal(a, b, m);
        break;
    case 2:
        if (m % 2 != 0) {
            std::cout << "Error: Simpson's Rule requires an even number of intervals.\n";
            return 1;
        }
        result = compositeSimpson(a, b, m);
        break;
    case 3:
        std::cout << "Gauss-Legendre Quadrature not available yet.\n";
        return 0; // Terminate execution as per the assignment
    default:
        std::cout << "Invalid selection.\n";
        return 1;
    }

    // Task 6
    std::cout << "Integration Result\n";
    std::cout << "------------------\n";
    std::cout << "Integration Interval: [" << a << ", " << b << "]\n";
    std::cout << "Number of Intervals (m): " << m << "\n";
    std::cout << "Quadrature Rule Used: ";
    if (selector == 1)
        std::cout << "Composite Trapezoidal Rule\n";
    else if (selector == 2)
        std::cout << "Composite Simpson's Rule\n";
    else
        std::cout << "Gaussian Quadrature\n";
    std::cout << "Actual value of integration of exp(x) over [0,4] is: 53.59819.\n";
    std::cout << "Computed Value of the Integral: " << result << "\n";

    return 0;
}
