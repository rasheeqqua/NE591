/*--------------------Inlab4--------------------*/
/*
 * Created by Hasibul Hossain Rasheeq on Thursday, 02/06/25.
*/

#include <fstream>
#include <vector>
#include <iomanip>
#include <sstream>
#include <cmath>
#include "LUFactorization.cpp"
#include "substitution.cpp"

// Matrix-vector multiplication
std::vector<double> matrixVectorProduct(const std::vector<std::vector<double>>& A,
                                      const std::vector<double>& x) {
    int n = A.size();
    std::vector<double> result(n, 0.0);

    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            result[i] += A[i][j] * x[j];
        }
    }
    return result;
}

// Calculate residual
double calculateResidual(const std::vector<std::vector<double>>& A,
                        const std::vector<double>& x,
                        const std::vector<double>& b) {
    std::vector<double> Ax = matrixVectorProduct(A, x);
    double maxResidual = 0.0;

    for (int i = 0; i < b.size(); ++i) {
        double residual = std::abs(Ax[i] - b[i]);
        maxResidual = std::max(maxResidual, residual);
    }
    return maxResidual;
}

int main() {
    std::ifstream inputFile("../input.txt");
    std::ofstream outputFile("output.txt");

    outputFile << "LU Factorization Solution Program\n";
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

    // Perform LU factorization
    std::vector<std::vector<double>> L(n, std::vector<double>(n));
    std::vector<std::vector<double>> U(n, std::vector<double>(n));
    if (!luFactorize(A, L, U, usePivoting)) {
        outputFile << "Error: LU factorization failed.\n";
        return 1;
    }

    // Solve the system
    std::vector<double> y(n), x(n);
    forwardSubstitution(L, b, y);
    backSubstitution(U, y, x);

    // Calculate residual
    double maxResidual = calculateResidual(A, x, b);

    // Output results
    outputFile << "Solution Vector x:\n";
    for (int i = 0; i < n; ++i) {
        outputFile << "x[" << i + 1 << "] = " << std::setw(10) << x[i] << "\n";
    }
    outputFile << "\nMaximum Absolute Residual: " << maxResidual << "\n";

    inputFile.close();
    outputFile.close();

    return 0;
}