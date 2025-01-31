/*--------------------Inlab4--------------------*/
/*
 * Created by Hasibul Hossain Rasheeq on Friday, 1/31/25.
*/

#include <fstream>
#include <vector>
#include <iomanip>
#include "substitution.cpp"

bool checkInputData(int n, const std::vector<std::vector<double>>& U, std::ofstream& outputFile) {
    bool valid = true;
    for (int i = 0; i < n; ++i) {
        if (U[i][i] == 0.0) {
            outputFile << "Error: Zero detected on the diagonal of U at position (" << i + 1 << "," << i + 1 << ").\n";
            valid = false;
        }
    }
    return valid;
}

int main() {
    std::ifstream inputFile("../input.txt");
    std::ofstream outputFile("output.txt");

    // Task 1: Write the header to the output file
    outputFile << "LU Factorization Solution Program\n";
    outputFile << "Author: Hasibul H. Rasheeq\n";
    outputFile << "Affiliation: NC State University\n";
    outputFile << "Date: " << __DATE__ << "\n";
    outputFile << "----------------------------------------\n\n";

    if (!inputFile) {
        outputFile << "Error: Cannot open input file.\n";
        return 1;
    }

    // Task 2: Read the input text file
    int n;
    inputFile >> n;

    // Task 3: Input validation
    if (n <= 0) {
        outputFile << "Error: Matrix order must be a positive integer.\n";
        return 1;
    }

    // Task 2: Initialize L and U matrices
    std::vector<std::vector<double>> L(n, std::vector<double> (n, 0.0));
    std::vector<std::vector<double>> U(n, std::vector<double> (n,0.0));
    std::vector<double> b(n);

    // Task 2: Set diagonal elements of L to 1
    for (int i = 0; i < n; ++i) {
        L[i][i] = 1.0;
    }

    // Task 2: Read non-zero elements of L (below the diagonal)
    for (int i = 1; i < n; ++i) {
        for (int j = 0; j < i; ++j) {
            inputFile >> L[i][j];
        }
    }

    // Task 2: Read non-zero elements of U (on and above the diagonal)
    for (int i = 0; i < n; ++i) {
        for (int j = i; j < n; ++j) {
            inputFile >> U[i][j];
        }
    }

    // Task 2: Read Right-Hand-Side (RHS) vector b
    for (int i = 0; i < n; ++i) {
        inputFile >> b[i];
    }

    // Task 3: Check the correctness of the input data
    if (!checkInputData(n, U, outputFile)) {
        outputFile << "Error: Invalid input data detected.\n";
        return 1;
    } else {
        outputFile << "Congratulations! All input data are correct.\n\n";
    }

    // Task 3: Echo the input data to the output file
    outputFile << "Matrix order: n = " << n << "\n\n";

    outputFile << "Lower Triangular Matrix L:\n";
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            outputFile << std::setw(10) << L[i][j] << " ";
        }
        outputFile << "\n";
    }
    outputFile << "\n";

    outputFile << "Upper Triangular Matrix U:\n";
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            outputFile << std::setw(10) << U[i][j] << " ";
        }
        outputFile << "\n";
    }
    outputFile << "\n";

    outputFile << "Right-Hand Side Vector b:\n";
    for (int i = 0; i < n; ++i) {
        outputFile << std::setw(10) << b[i] << "\n";
    }
    outputFile << "\n";

    // Task 4: Solve Ly = b using forward substitution
    std::vector<double> y(n);
    forwardSubstitution(L, b, y);

    // Task 4: Solve Ux = y using back substitution
    std::vector<double> x(n);
    backSubstitution(U, y, x);

    // Task 5: Write the solution vector to the output file
    outputFile << "Solution Vector x:\n";
    for (int i = 0; i < n; ++i) {
        outputFile << "x[" << i + 1 << "] = " << std::setw(10) << x[i] << "\n";
    }

    // Close files
    inputFile.close();
    outputFile.close();

    return 0;
}