//
// Created by Hasibul H. Rasheeq on 1/30/25.
//
/*--------------------Numerical Experiment--------------------*/
#include <iostream>
#include <fstream>
#include <vector>
#include <iomanip>
#include <cmath>
#include "CompositeSimpson.cpp"
#include "CompositeTrapezoidal.cpp"
#include "ExactIntegral.cpp"
#include "GaussianQuadrature.cpp"

int main() {
    // Open the output file
    std::ofstream outfile("NumericalOutput.txt");
    if (!outfile) {
        std::cerr << "Error: Cannot open NumericalOutput.txt for writing." << std::endl;
        return 1;
    }

    double a = -1.0;
    double b = 1.0;

    double exactValue = exactIntegral(a, b);

    // Vectors to store m or n values
    std::vector<int> m_values = {2, 4, 8, 16, 32};
    std::vector<int> n_values = {2, 3, 4, 5};

    // Headers
    outfile << "Results for Composite Trapezoidal Rule\n";
    outfile << "m\tApproximate Value\tIntegration Error\n";

    // Trapezoidal Rule
    for (int m : m_values) {
        double approx = compositeTrapezoidal(a, b, m);
        double error = approx - exactValue;
        outfile << m << "\t" << std::setprecision(10) << approx << "\t" << error << "\n";
    }

    outfile << "\nResults for Composite Simpson's Rule\n";
    outfile << "m\tApproximate Value\tIntegration Error\n";

    // Simpson's Rule
    for (int m : m_values) {
        if (m % 2 != 0) continue; // Simpson's Rule requires even m
        double approx = compositeSimpson(a, b, m);
        double error = approx - exactValue;
        outfile << m << "\t" << std::setprecision(10) << approx << "\t" << error << "\n";
    }

    outfile << "\nResults for Gauss-Legendre Quadrature\n";
    outfile << "n\tApproximate Value\tIntegration Error\n";

    // Gauss-Legendre Quadrature
    for (int n : n_values) {
        double approx = gaussLegendreQuadrature(a, b, n);
        double error = approx - exactValue;
        outfile << n << "\t" << std::setprecision(10) << approx << "\t" << error << "\n";
    }

    outfile.close();

    std::cout << "Numerical experiments completed. Results saved in NumericalOutput.txt" << std::endl;

    return 0;
}
