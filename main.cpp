/*--------------------Outlab5--------------------*/
/*
 * Created by Hasibul Hossain Rasheeq on 02/08/25.
*/

#include <algorithm>
#include <fstream>
#include <vector>
#include <iomanip>
#include <iostream>
#include <sstream>
#include "ApplyPermutationMatrix.cpp"
#include "CalculateResiduals.cpp"
#include "LUFactorization.cpp"
#include "LUPFactorization.cpp"
#include "substitution.cpp"

void printMatrix(std::ofstream& outputFile, const std::vector<std::vector<double>>& matrix, const std::string& name) {
    outputFile << name << ":\n";
    for (const auto& row : matrix) {
        for (double val : row) {
            outputFile << std::setw(10) << val << " ";
        }
        outputFile << "\n";
    }
    outputFile << "\n";
}

int main() {
    std::ifstream inputFile("../input.txt");
    std::ofstream outputFile("output.txt");

    outputFile << "LUP Factorization Solution Program\n";
    outputFile << "Author: Hasibul H. Rasheeq\n";
    outputFile << "Date: " << __DATE__ << "\n";
    outputFile << "----------------------------------------\n\n";

    if (!inputFile) {
        outputFile << "Error: Cannot open input file.\n";
        return 1;
    }

    std::string line;
    // Read matrix order and pivoting flag
    std::getline(inputFile, line);
    std::stringstream ss(line);
    int n;
    bool usePivoting;
    ss >> n >> usePivoting;

    if (n <= 0) {
        outputFile << "Error: Matrix order must be a positive integer.\n";
        return 1;
    }

    // Read matrix A
    std::vector<std::vector<double>> A(n, std::vector<double>(n));
    for (int i = 0; i < n; ++i) {
        std::getline(inputFile, line);
        std::stringstream rowStream(line);
        for (int j = 0; j < n; ++j) {
            rowStream >> A[i][j];
        }
    }

    // Read vector b
    std::vector<double> b(n);
    for (int i = 0; i < n; ++i) {
        std::getline(inputFile, line);
        b[i] = std::stod(line);
    }

    // Print input matrix A
    printMatrix(outputFile, A, "Input Matrix A");

    // Perform either LU factorization or else LUP factorization
    std::vector<std::vector<double>> L(n, std::vector<double>(n));
    std::vector<std::vector<double>> U(n, std::vector<double>(n));
    std::vector<std::vector<double>> P(n, std::vector<double>(n));
    std::vector<double> y(n), x(n), Pb(n);
    if (!usePivoting) {
        if (!luFactorize(A, L, U)) {
            outputFile << "Error: LU factorization failed.\n";
            return 1;
        }

        // Print L and U matrices
        printMatrix(outputFile, L, "Lower Triangular Matrix L");
        printMatrix(outputFile, U, "Upper Triangular Matrix U");

        forwardSubstitution(L, b, y);
        backSubstitution(U, y, x);
    } else if (usePivoting) {
        if (!lupFactorize(A, L, U, P)) {
            outputFile << "Error: LUP Factorization failed.\n";
            return 1;
        }

        // Print L and U matrices
        printMatrix(outputFile, L, "Lower Triangular Matrix L");
        printMatrix(outputFile, U, "Upper Triangular Matrix U");
        printMatrix(outputFile, P, "Upper Triangular Matrix U");

        applyPermutationMatrix(P, b, Pb);

        forwardSubstitution(L, Pb, y);
        backSubstitution(U, y, x);
    }

    // Print solution vector
    outputFile << "Solution Vector x:\n";
    for (int i = 0; i < n; ++i) {
        outputFile << "x[" << i + 1 << "] = " << std::setw(10) << x[i] << "\n";
    }
    outputFile << "\n";

    // Calculate and print residuals
    std::vector<double> residuals = calculateResiduals(A, x, b);
    outputFile << "Residual errors (Ax-b):\n";
    for (double residual : residuals) {
        outputFile << std::fixed << std::setprecision(1) << residual << " ";
    }
    outputFile << "\n\n";

    // Print maximum residual
    double maxResidual = *std::max_element(residuals.begin(), residuals.end());
    outputFile << "Maximum Absolute Residual: " << std::scientific << maxResidual << "\n";

    inputFile.close();
    outputFile.close();

    return 0;
}