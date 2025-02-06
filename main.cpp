/*--------------------Inlab4--------------------*/
/*
 * Created by Hasibul Hossain Rasheeq on Thursday, 02/06/25.
*/

#include <algorithm>
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

std::vector<double> calculateResiduals(const std::vector<std::vector<double>>& A,
                                     const std::vector<double>& x,
                                     const std::vector<double>& b) {
    std::vector<double> Ax = matrixVectorProduct(A, x);
    std::vector<double> residuals(b.size());
    for (int i = 0; i < b.size(); ++i) {
        residuals[i] = std::abs(Ax[i] - b[i]);
    }
    return residuals;
}

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

    // Print input matrix A
    printMatrix(outputFile, A, "Input Matrix A");

    // Perform LU factorization
    std::vector<std::vector<double>> L(n, std::vector<double>(n));
    std::vector<std::vector<double>> U(n, std::vector<double>(n));
    if (!luFactorize(A, L, U, usePivoting)) {
        outputFile << "Error: LU factorization failed.\n";
        return 1;
    }

    // Print L and U matrices
    printMatrix(outputFile, L, "Lower Triangular Matrix L");
    printMatrix(outputFile, U, "Upper Triangular Matrix U");

    // Solve the system
    std::vector<double> y(n), x(n);
    forwardSubstitution(L, b, y);
    backSubstitution(U, y, x);

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