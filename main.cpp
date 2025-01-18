#include <algorithm>
#include <fstream>
#include <iostream>
#include <map>
#include <sstream>
#include <string>
#include <vector>
#include "function.cpp" // for evaluate()

int main() {
    std::cout << "          NE-591 IN-LAB ASSIGNMENT 2          " << std::endl;
    std::cout << "          Credits: Hasibul H. Rasheeq, NCSU          " << std::endl;
    std::cout << "          Jan 17, 2025          " << std::endl;
    std::cout << "\n\n          LAGRANGE INTERPOLATION PROGRAM          " << std::endl;
    std::cout << "************************************************ \n\n" << std::endl;

    int n, m, choice;
    std::vector<double> x;
    std::vector<double> y;
    std::map<double, double> xy;

    std::cout << "Number of interpolation points, n: ";
    std::cin >> n;
    std::cout << "Number of evaluation points, m: ";
    std::cin >> m;
    if (n < 1 || m < 1) {
        std::cout << "You cannot provide non-positive values for m and n!" << std::endl;
        return 1;
    }

    std::cout << "Please enter your " << n << " interpolation points, x: " << std::endl;
    for (int i = 0; i < n; i++) {
        double temp;
        std::cin >> temp;
        x.push_back(temp);
    }
    std::sort(x.begin(), x.end());

    std::cout << "\nOption 1: Read ALL (x, y) pairs from a text file.\n";
    std::cout << "Option 2: Use the user-provided x-values and compute y = f(x) (currently exp(x)).\n";
    std::cout << "P.S.: You can change f(x) by editing function.cpp.\n";
    std::cout << "Enter your choice: ";
    std::cin >> choice;

    if (choice == 1) {
        std::string filename;
        std::cout << "Enter the name of the file: ";
        std::cin >> filename;

        std::ifstream file(filename);
        if (!file.is_open()) {
            std::cout << "Cannot open the file. Check name or file integrity." << std::endl;
            return 1;
        }

        // We assume the file has lines like: -1.0, 0.37 (x, y)
        std::string line;
        while (std::getline(file, line)) {
            std::stringstream ss(line);
            double x_in_file, y_in_file;
            char comma;

            if (ss >> x_in_file) {
                if (ss.peek() == ',') {
                    ss >> comma;
                }
                if (ss >> y_in_file) {
                    xy[x_in_file] = y_in_file;
                }
            }
        }
        file.close();

    } else if (choice == 2) {
        y = evaluate(x);
        for (size_t i = 0; i < x.size(); i++) {
            xy[x[i]] = y[i];
        }
    } else {
        std::cout << "Invalid choice." << std::endl;
        return 1;
    }

    return 0;
}
