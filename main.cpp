/*--------------------Out-Lab Assignment 3--------------------*/
/*
 * Created by Hasibul H. Rasheeq on Thursday, 1/30/25.
*/
#include <iostream>
#include <ctime>
#include <fstream> // Included for file operations
#include "CompositeSimpson.cpp"
#include "CompositeTrapezoidal.cpp"
#include "ExactIntegral.cpp"
#include "GaussianQuadrature.cpp"

int main() {
    // Open the output file
    std::ofstream outfile("output.txt");
    if (!outfile) {
        std::cerr << "Error: Cannot open output.txt for writing." << std::endl;
        return 1;
    }

    // Write initial information to the output file
    outfile << "Numerical Integration Program\n";
    outfile << "Author: Hasibul H. Rasheeq, NC State University\n";

    // Task 1
    time_t now = time(0);
    char* dt = ctime(&now);
    outfile << "Date: " << dt << "\n";

    outfile << "This program computes the approximate integral using numerical quadrature rules\n";
    outfile << "-----------------------------------------------------------------------------------\n\n";

    // Task 2
    double a, b;
    int m, selector;

    // User inputs remain on the console
    std::cout << "Enter the lower limit of integration (a): ";
    std::cin >> a;

    std::cout << "Enter the upper limit of integration (b): ";
    std::cin >> b;

    if (a >= b) {
        std::cout << "Error: The lower limit 'a' must be less than the upper limit 'b'.\n";
        return 1;
    }

    std::cout << "Enter the number of intervals (m): ";
    std::cin >> m;

    if (m <= 0) {
        std::cout << "Error: The number of intervals 'm' must be a positive integer.\n";
        return 1;
    }

    std::cout << "Select the Quadrature Rule:\n";
    std::cout << "1. Composite Trapezoidal Rule\n";
    std::cout << "2. Composite Simpson's Rule\n";
    std::cout << "3. Gaussian Quadrature\n";
    std::cout << "Enter your choice (1-3): ";
    std::cin >> selector;
    if (selector < 1 || selector > 3) {
        std::cout << "Error: Invalid selection for the Quadrature Rule.\n";
        return 1;
    }

    // Write confirmation and details first to the screen and then to the output file
    std::cout << "All input data is correct. The output will be saved in output.txt file inside the build folder" <<
        std::endl;
    outfile << "All input data is correct. Proceeding with the computation.\n\n";
    outfile << "Integration Interval: [" << a << ", " << b << "]\n";
    outfile << "Number of Intervals (m): " << m << "\n";
    outfile << "Selected Quadrature Rule: ";
    switch (selector) {
        case 1:
            outfile << "Composite Trapezoidal Rule\n";
            break;
        case 2:
            outfile << "Composite Simpson's Rule\n";
            break;
        case 3:
            outfile << "Gauss-Legendre Quadrature\n";
            break;
    }
    outfile << "----------------------------------------------------------\n\n";

    // Task 4 implementation is in `integrand.cpp` file.

    // Task 5
    double result = 0.0;

    switch (selector) {
    case 1:
        result = compositeTrapezoidal(a, b, m);
        break;
    case 2:
        if (m % 2 != 0) {
            outfile << "Error: Simpson's Rule requires an even number of intervals.\n";
            return 1;
        }
        result = compositeSimpson(a, b, m);
        break;
    case 3:
        result = gaussLegendreQuadrature(a, b, m);
        break;
    }

    double exactValue = exactIntegral(a, b);

    // Task 6
    outfile << "Integration Result\n";
    outfile << "------------------\n";
    outfile << "Integration Interval: [" << a << ", " << b << "]\n";
    outfile << "Number of Intervals (m): " << m << "\n";
    outfile << "Quadrature Rule Used: ";
    switch (selector) {
    case 1:
        outfile << "Composite Trapezoidal Rule\n";
        break;
    case 2:
        outfile << "Composite Simpson's Rule\n";
        break;
    case 3:
        outfile << "Gauss-Legendre Quadrature\n";
        break;
    }
    outfile << "Actual value of the integral: " << exactValue << "\n";
    outfile << "Computed value of the integral: " << result << "\n";
    outfile << "Error = Calculated value - Actual value = " << result - exactValue << "\n";

    // Close the output file
    outfile.close();

    return 0;
}
