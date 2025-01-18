#include <algorithm>
#include <fstream>
#include <iostream>
#include <map>
#include <string>
#include <vector>
#include <cmath>
#include <iomanip>
#include "function.cpp"
#include "polynomial.cpp"

// Step #8
std::vector<double> createEvaluationPoints(int m, double a, double b) {
    std::vector<double> evalPoints(m);
    if (m == 1) {
        evalPoints[0] = a;
        return evalPoints;
    }
    double step = (b - a) / (m - 1);
    for (int i = 0; i < m; ++i) {
        evalPoints[i] = std::round((a + i * step) * 10) / 10.0;
    }
    return evalPoints;
}

int main() {
    // Step #1
    std::cout << "*********************************************************\n";
    std::cout << "    NE-591 IN-LAB ASSIGNMENT 2\n";
    std::cout << "    Author  : Hasibul H. Rasheeq, NCSU\n";
    std::cout << "    Date    : 1/17/25\n";
    std::cout << "    Purpose : Lagrange Interpolation\n";
    std::cout << "*********************************************************\n\n";

    // Step #2
    int n, m;
    std::cout << "Enter number of interpolation points (n): ";
    std::cin >> n;
    std::cout << "Enter number of evaluation points (m): ";
    std::cin >> m;

    // Step #3
    if (n < 1 || m < 1) {
        std::cerr << "ERROR: You cannot provide non-positive values for n or m.\n";
        return 1;
    }

    // Step #4
    std::vector<double> userX(n);
    std::cout << "Please enter your " << n << " interpolation points in ascending order:\n";
    for (int i = 0; i < n; ++i) {
        std::cin >> userX[i];
    }
    std::sort(userX.begin(), userX.end());

    // Step #5
    int choice;
    std::cout << "\nOption 1: Read ALL (x, y) pairs from a text file.\n";
    std::cout << "Option 2: Use the user-provided x-values and evaluate f(x).\n";
    std::cout << "Enter your choice (1 or 2): ";
    std::cin >> choice;

    std::map<double, double> xy;
    double a, b;

    // Step #6
    if (choice == 1) {
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
            xy[xx] = yy;
        }
        file.close();
        if (xy.empty()) {
            std::cerr << "ERROR: The file appears to have no valid (x,y) data!\n";
            return 1;
        }

        a = xy.begin()->first;
        b = xy.rbegin()->first;
    } else if (choice == 2) {
        std::vector<double> computedY = evaluate(userX);
        for (size_t i = 0; i < userX.size(); ++i) {
            xy[userX[i]] = computedY[i];
        }
        std::cout << "\nValues computed from the function:\n";
        for (auto &kv : xy) {
            std::cout << "x = " << kv.first << ", y = " << kv.second << "\n";
        }
        a = userX.front();
        b = userX.back();
    } else {
        std::cerr << "ERROR: Invalid choice.\n";
        return 1;
    }

    // Step #7
    std::cout << "\nAll input data looks correct. Proceeding...\n";
    std::cout << "Your n = " << n << ", m = " << m << "\n";

    if (choice == 1) {
        std::cout << "You chose Option 1: we have " << xy.size() << " points read from file.\n";
    } else {
        std::cout << "You chose Option 2: we have " << xy.size() << " points from userX/evaluate.\n";
    }

    std::cout << "Interpolation interval is [" << a << ", " << b << "]\n\n";
    std::cout << "Lagrange Interpolation will be implemented here (placeholder)\n\n";

    // Step #8
    std::vector<double> xEval = createEvaluationPoints(m, a, b);

    // Step #9
    const int wIndex = 5;
    const int wX = 12;
    const int wY = 12;
    const int wInterp = 12;
    const int wError = 12;
    if (choice == 1) {
        std::cout << "--------------------------------------------------------------\n";
        std::cout << std::setw(wIndex) << "i"
                  << std::setw(wX) << "x(i)"
                  << std::setw(wY) << "y(file)"
                  << std::setw(wInterp) << "Interp"
                  << std::setw(wError) << "Error"
                  << "\n";
        std::cout << "--------------------------------------------------------------\n";
        for (int i = 0; i < m; ++i) {
            double xx = xEval[i];
            std::cout << std::setw(wIndex) << (i + 1);
            std::cout << std::setw(wX) << xx;
            auto it = xy.find(xx);
            if (it != xy.end()) {
                std::cout << std::setw(wY) << it->second;
            } else {
                std::cout << std::setw(wY) << " ";
            }
            double Lx = lagrangeInterpolate(xx);
            std::cout << std::setw(wInterp) << Lx;
            std::cout << std::setw(wError) << "-";
            std::cout << "\n";
        }
    } else {
        std::cout << "---------------------------------------------\n";
        std::cout << std::setw(wIndex) << "i"
                  << std::setw(wX) << "x(i)"
                  << std::setw(wY) << "f(x)"
                  << std::setw(wInterp) << "Interp"
                  << "\n";
        std::cout << "---------------------------------------------\n";
        for (int i = 0; i < m; ++i) {
            double xx = xEval[i];
            std::cout << std::setw(wIndex) << (i + 1);
            std::cout << std::setw(wX) << xx;
            double fx = std::exp(xx);
            std::cout << std::setw(wY) << fx;
            double Lx = lagrangeInterpolate(xx);
            std::cout << std::setw(wInterp) << Lx;
            std::cout << "\n";
        }
    }

    // Step #10
    std::cout << "\nProgram finished.\n";
    return 0;
}
