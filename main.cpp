#include <algorithm>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include <cmath>
#include <iomanip>
#include "polynomial.cpp"
#include "function.cpp"

// Create m evaluation points between a and b
std::vector<double> createEvaluationPoints(int m, double a, double b) {
    std::vector<double> evalPoints(m);
    if (m == 1) {
        evalPoints[0] = a;
        return evalPoints;
    }
    double step = (b - a) / (m - 1);
    for (int i = 0; i < m; ++i) {
        evalPoints[i] = a + i * step;
    }
    return evalPoints;
}

int main() {
    // Step #1: Program introduction
    std::cout << "*********************************************************\n";
    std::cout << "    NE-591 IN-LAB ASSIGNMENT 2\n";
    std::cout << "    Author  : Hasibul H. Rasheeq, NCSU\n";
    std::cout << "    Date    : 1/17/25\n";
    std::cout << "    Purpose : Lagrange Interpolation\n";
    std::cout << "*********************************************************\n\n";

    // Step #2: Get n and m from the user
    int n, m;
    std::cout << "Enter number of interpolation points (n): ";
    std::cin >> n;
    std::cout << "Enter number of evaluation points (m): ";
    std::cin >> m;

    // Step #3: Validate n and m
    if (n < 1 || m < 1) {
        std::cerr << "ERROR: You cannot provide non-positive values for n or m.\n";
        return 1;
    }

    // Step #4: Get the x_i values from the user
    std::vector<double> userX(n);
    std::cout << "Please enter your " << n << " interpolation points in ascending order:\n";
    for (int i = 0; i < n; ++i) {
        std::cin >> userX[i];
    }
    std::sort(userX.begin(), userX.end());

    // Step #5: Choose data input method
    int choice;
    std::cout << "\nOption 1: Read ALL (x, y) pairs from a text file.\n";
    std::cout << "Option 2: Use the user-provided x-values and evaluate f(x).\n";
    std::cout << "Enter your choice (1 or 2): ";
    std::cin >> choice;

    std::vector<double> x; // x_i values
    std::vector<double> y; // y_i values
    double a, b; // interval [a, b]

    // Step #6: Read y_i values
    if (choice == 1) {
        std::vector<std::pair<double, double>> dataPairs;
        std::string filename;
        std::cout << "Enter the name of the file: ";
        std::cin >> filename;
        std::ifstream file(filename);
        if (!file.is_open()) {
            std::cerr << "ERROR: Cannot open the file: " << filename << "\n";
            return 1;
        }
        std::string line;
        while (std::getline(file, line)) {
            std::size_t pos = line.find(',');
            if (pos == std::string::npos) {
                continue;
            }
            std::string left = line.substr(0, pos);
            std::string right = line.substr(pos + 1);
            double xx = std::stod(left);
            double yy = std::stod(right);
            dataPairs.emplace_back(xx, yy);
        }
        file.close();

        if (dataPairs.empty()) {
            std::cerr << "ERROR: The file appears to have no valid (x,y) data!\n";
            return 1;
        }

        // Sort data pairs by x values
        std::sort(dataPairs.begin(), dataPairs.end());

        // Separate x and y values
        for (const auto& pair : dataPairs) {
            x.push_back(pair.first);
            y.push_back(pair.second);
        }

        a = x.front();
        b = x.back();
    } else if (choice == 2) {
        x = userX; // x_i values provided by the user
        y = evaluate(userX); // y_i values computed from the function

        std::cout << "\nValues computed from the function:\n";
        for (size_t i = 0; i < x.size(); ++i) {
            std::cout << "x = " << x[i] << ", y = " << y[i] << "\n";
        }

        a = x.front();
        b = x.back();
    } else {
        std::cerr << "ERROR: Invalid choice.\n";
        return 1;
    }

    // Step #7: Confirm data and proceed
    std::cout << "\nAll input data looks correct. Proceeding...\n";
    std::cout << "Your n = " << n << ", m = " << m << "\n";

    if (choice == 1) {
        std::cout << "You chose Option 1: we have " << x.size() << " points read from file.\n";
    } else {
        std::cout << "You chose Option 2: we have " << x.size() << " points from userX/evaluate.\n";
    }

    std::cout << "Interpolation interval is [" << a << ", " << b << "]\n\n";
    std::cout << "Lagrange Interpolation will be implemented here.\n\n";

    // Step #8: Create evaluation points
    std::vector<double> xEval = createEvaluationPoints(m, a, b);

    // Step #9: Compute and display interpolated values and errors
    const int wIndex = 5;
    const int wX = 12;
    const int wY = 12;
    const int wInterp = 12;
    const int wError = 15; // Increased width to accommodate scientific notation
    const double EPSILON = 1e-9; // For floating-point comparison

    std::cout << "------------------------------------------------------------------\n";
    std::cout << std::setw(wIndex) << "i"
              << std::setw(wX) << "x(i)"
              << std::setw(wY) << (choice == 1 ? "y(file)" : "f(x)")
              << std::setw(wInterp) << "Interp"
              << std::setw(wError) << "Error"
              << "\n";
    std::cout << "------------------------------------------------------------------\n";

    for (int i = 0; i < m; ++i) {
        double xx = xEval[i];
        std::cout << std::setw(wIndex) << (i + 1);
        std::cout << std::setw(wX) << xx;

        double Lx = lagrangeInterpolate(xx, x, y); // Interpolated value

        if (choice == 2) {
            // Option 2: Compute f(xx) and error for all points
            double yActual = exp(xx);
            double error = yActual - Lx;
            std::cout << std::setw(wY) << yActual;
            std::cout << std::setw(wInterp) << Lx;
            std::cout << std::setw(wError) << std::scientific << std::setprecision(6) << error << "\n";
            std::cout << std::fixed;
        } else {
            // Option 1: Only compute error if xx matches x_i
            bool isDataPoint = false;
            double yActual = 0.0;
            for (size_t k = 0; k < x.size(); ++k) {
                if (std::abs(xx - x[k]) < EPSILON) {
                    isDataPoint = true;
                    yActual = y[k];
                    break;
                }
            }

            if (isDataPoint) {
                double error = yActual - Lx;
                std::cout << std::setw(wY) << yActual;
                std::cout << std::setw(wInterp) << Lx;
                std::cout << std::setw(wError) << std::scientific << std::setprecision(6) << error << "\n";
                std::cout << std::fixed;
            } else {
                std::cout << std::setw(wY) << " ";
                std::cout << std::setw(wInterp) << Lx;
                std::cout << std::setw(wError) << "-" << "\n";
            }
        }
    }

    // Step #10: Finish
    std::cout << "\nProgram finished.\n";
    return 0;
}
