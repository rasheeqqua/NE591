/*--------------------Inlab5--------------------*/
/*
 * Created by Hasibul Hossain Rasheeq on Friday, 02/07/25.
*/

#include <fstream>
#include <vector>
#include <iomanip>
#include <string>
#include "substitution.cpp"
#include "RowSwap.cpp"

bool checkInputData(int n,
                    const std::vector<std::vector<double>>& L,
                    const std::vector<std::vector<double>>& U,
                    const std::vector<std::vector<double>>& P,
                    std::ofstream& outputFile) {
    bool valid = true;

    // Check L matrix diagonal
    for (int i = 0; i < n; ++i) {
        if (L[i][i] != 1.0) {
            outputFile << "Error: L matrix diagonal must be 1.0\n";
            valid = false;
        }
    }

    // Check U matrix diagonal
    for (int i = 0; i < n; ++i) {
        if (U[i][i] == 0.0) {
            outputFile << "Error: Zero detected on the diagonal of U\n";
            valid = false;
        }
    }

    // Check Permutation Matrix
    for (int i = 0; i < n; ++i) {
        double rowSum = 0;
        double colSum = 0;
        for (int j = 0; j < n; ++j) {
            rowSum += P[i][j];
            colSum += P[j][i];
        }
        if (rowSum != 1.0 || colSum != 1.0) {
            outputFile << "Error: Invalid Permutation Matrix\n";
            valid = false;
            break;
        }
    }

    return valid;
}

void printMatrix(std::ofstream& outputFile,
                 const std::vector<std::vector<double>>& matrix,
                 const std::string& matrixName) {
    outputFile << matrixName << ":\n";
    for (const auto& row : matrix) {
        for (const auto& elem : row) {
            outputFile << std::setw(10) << elem << " ";
        }
        outputFile << "\n";
    }
    outputFile << "\n";
}

void printVector(std::ofstream& outputFile,
                 const std::vector<double>& vec,
                 const std::string& vectorName) {
    outputFile << vectorName << ":\n";
    for (const auto& elem : vec) {
        outputFile << std::setw(10) << elem << "\n";
    }
    outputFile << "\n";
}

int main() {
    std::ifstream inputFile("../input.txt");
    std::ofstream outputFile("output.txt");

    // Task 1: Write header
    outputFile << "LUP Factorization Solution Program\n";
    outputFile << "Author: Hasibul H. Rasheeq\n";
    outputFile << "Affiliation: NC State University\n";
    outputFile << "Date: " << __DATE__ << "\n";
    outputFile << "----------------------------------------\n\n";

    if (!inputFile) {
        outputFile << "Error: Cannot open input file.\n";
        return 1;
    }

    std::string line;
    std::getline(inputFile, line);
    std::stringstream ss(line);
    int n;
    ss >> n;

    // Input validation
    if (n <= 0) {
        outputFile << "Error: Matrix order must be positive.\n";
        return 1;
    }

    // Read L matrix entirely
    std::vector<std::vector<double>> L(n, std::vector<double>(n));
    for (int i = 0; i < n; ++i) {
        std::getline(inputFile, line);
        std::stringstream rowStream(line);
        for (int j = 0; j < n; ++j) {
            rowStream >> L[i][j];
        }
    }

    // Read U matrix entirely
    std::vector<std::vector<double>> U(n, std::vector<double>(n));
    for (int i = 0; i < n; ++i) {
        std::getline(inputFile, line);
        std::stringstream rowStream(line);
        for (int j = 0; j < n; ++j) {
            rowStream >> U[i][j];
        }
    }

    // Read Permutation Matrix
    std::vector<std::vector<double>> P(n, std::vector<double>(n));
    for (int i = 0; i < n; ++i) {
        std::getline(inputFile, line);
        std::stringstream rowStream(line);
        for (int j = 0; j < n; ++j) {
            rowStream >> P[i][j];
        }
    }

    // Read RHS vector b
    std::vector<double> b(n);
    for (int i = 0; i < n; ++i) {
        std::getline(inputFile, line);
        b[i] = std::stod(line);
    }

    // Task 3: Check input data
    if (!checkInputData(n, L, U, P, outputFile)) {
        outputFile << "Error: Invalid input data.\n";
        return 1;
    }

    // Task 4: Apply Permutation Matrix to b
    std::vector<double> Pb(n);
    applyPermutationMatrix(P, b, Pb);

    // Task 5: Solve system
    std::vector<double> y(n);
    std::vector<double> x(n);

    forwardSubstitution(L, Pb, y);
    backSubstitution(U, y, x);

    printMatrix(outputFile, L, "Lower Triangular Matrix L");
    printMatrix(outputFile, U, "Upper Triangular Matrix U");
    printMatrix(outputFile, P, "Permutation Matrix P");
    printVector(outputFile, b, "Right-Hand Side Vector b");

    // Task 6: Write solution vector
    outputFile << "Solution Vector x:\n";
    for (int i = 0; i < n; ++i) {
        outputFile << "x[" << i + 1 << "] = " << std::setw(10) << x[i] << "\n";
    }

    inputFile.close();
    outputFile.close();

    return 0;
}