/*--------------------In-Lab-2--------------------*/
/*
 * Created by Hasibul H. Rasheeq on Friday, 1/17/25.
*/

#include <algorithm>
#include <functional>
#include <bits/stdc++.h>
#include "function.cpp"
#include <iostream>
#include <string>
#include <vector>
#include <bits/std_function.h>

int main() {
    std::cout << "          NE-591 IN-LAB ASSIGNMENT 2          " << std::endl;
    std::cout << "          Credits: Hasibul H. Rasheeq, NCSU          " << std::endl;
    std::cout << "          Jan 17, 2025          " << std::endl;
    std::cout << "\n\n          LAGRANGE INTERPOLATION PROGRAM          " <<std::endl;
    std::cout << "************************************************ \n\n" << std::endl;

    int n, m, choice;
    std::vector<double> x, y, step_values;

    std::cout << "Number of interpolation points, n: ";
    std::cin >> n;
    std::cout << "Number of evaluation points, m: ";
    std::cin >> m;
    if (n < 1 || m < 1) {
        std::cout << "You cannot provide non-positive values for m and n!" << std::endl;
        return 1;
    }

    std::cout << "Interpolation points, x: " << std::endl;
    for (int i = 0; i < n; i++) {
        double temp;
        std::cin >> temp;
        x.push_back(temp);
    }
    std::sort(x.begin(), x.end());

    std::cout << "Option 1: Enter the values from a text file.";
    std::cout << "Option 2: Calculate the values from f(x)." << std::endl;
    std::cout << "P.S.: You can change the function by editing line 11 of function.cpp file" << std::endl;
    std::cout << "Enter your choice: ";
    std::cin >> choice;
    if (choice == 1) {
        std::string filename;
        std::cout << "Enter the name of the file: " << std::endl;
        std::cin >> filename;

        std::ifstream file(filename);
        if (!file.is_open()) {
            std::cout << "Sorry cannot open the file. Either the name of the file is wrong or it is corrupted." << std::endl;
            return 1;
        }

        std::string s;
        while (getline(file, s, '\n')) {
            double evaluated_value = std::stod(s);
            y.push_back(evaluated_value);
        }
    } else if (choice == 2) {
        y = evaluate(x);
    } else {
        std::cout << "Invalid choice." << std::endl;
        return 1;
    }

    std::cout << "Looks like all the provided values are correct so far! Give me a few moments to return the interpolated values!" << std::endl;

    return 0;
}
