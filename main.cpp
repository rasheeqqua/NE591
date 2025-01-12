/*--------------------Out-Lab Assignemnt 1--------------------*/
/*
 * Created by Hasibul Hossain Rasheeq on Sunday, 1/12/25.
*/

/*----------------------------VARIABLE DECLARATION----------------------------*/
/*
 * k - Stores the constant that will be multiplied with Matrix A.
 * M, N - Stores the number of rows and columns respectively, of matrix A and B.
 * N, J - Stores the number of rows and columns of matrix F.
 * A, B, F - Stores the generated matrices A, B, and F.
 * C - Stores the resultant matrix of operation A + B
 * D - Stores the resultant matrix of operation k.A
 * E - Stores the resultant matrix of operation A x F
 */
/*-------------------------VARIABLE DECLARATION ENDS--------------------------*/

#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <vector>

int main() {
    double k;
    int M, N, J;

    std::vector<std::vector<double>> A;
    std::vector<std::vector<double>> B;
    std::vector<std::vector<double>> F;

    std::vector<std::vector<double>> C;
    std::vector<std::vector<double>> D;
    std::vector<std::vector<double>> E;

    // Step 1: Prompting the user for providing the values of k, M, N and J
    std::cout << "Please enter a real value (k): " << std::endl;
    std::cin >> k;
    std::cout << "Please enter the dimensions (M, N, J) for the matrices: " << std::endl;
    std::cin >> M >> N >> J;

    // Step 2: Throwing an error in case one of the user provided dimensions is non-positive
    if (M <= 0 || N <= 0 || J <= 0) {
        throw std::invalid_argument("Invalid dimensions! Dimensions have to be bigger than 0");
    } else {
        std:: cout << "There are no errors in the input! Give me a few moments to perform the Matrix operations!" << std::endl;
    }

    // Step 3: Resize matrices according to the user provided dimensions
    A.resize(M, std::vector<double>(N));
    B.resize(M, std::vector<double>(N));
    C.resize(M, std::vector<double>(N));
    D.resize(M, std::vector<double>(N));
    E.resize(M, std::vector<double>(J));
    F.resize(N, std::vector<double>(J));

    // Step 4: Populate the matrices
    for (int m = 0; m < M; m++) {
        for (int n = 0; n < N; n++) {
            if (m == n) {
                A[m][n] = 1.0;
            } else {
                A[m][n] = 0.5;
            }
        }
    }

    for (int m = 0; m < M; m++) {
        for (int n = 0; n < N; n++) {
            if (m > n) {
                B[m][n] = 0.75;
            } else {
                B[m][n] = 0.25;
            }
        }
    }

    for (int n = 0; n < N; n++) {
        for (int j = 0; j < J; j++) {
            if (n + j == 0) {
                F[n][j] = 1.0/(n+j+1);
            } else {
                F[n][j] = 1.0/(n+j);
            }

        }
    }

    // Step 5: Performing matrix addition operation
    for (int m = 0; m < M; m++) {
        for (int n = 0; n < N; n++) {
            C[m][n] = A[m][n] + B[m][n];
        }
    }

    // Step 6: Performing scalar multiplication operation
    for (int m = 0; m < M; m++) {
        for (int n = 0; n < N; n++) {
            D[m][n] = k * A[m][n];
        }
    }

    // Step 7: Performing matrix multiplication operation
    for (int m = 0; m < M; m++) {
        for (int j = 0; j < J; j++) {
            double sum = 0;
            for (int n = 0; n < N; n++) {
                sum += A[m][n] * F[n][j];
            }
            E[m][j] = sum;
        }
    }

    // Step 8: Displaying the results in a professional way
    std::cout << "------------------------------RESULTS-------------------------" << std::endl;
    std::cout << "Matrix C (A + B):" << std::endl;
    for (int m = 0; m < M; m++) {
        for (int n = 0; n < N; n++) {
            std::cout << C[m][n] << " ";
        }
        std::cout << std::endl;
    }

    std::cout << "\nMatrix D (k * A):" << std::endl;
    for (int m = 0; m < M; m++) {
        for (int n = 0; n < N; n++) {
            std::cout << std::setprecision(2) << D[m][n] << " ";
        }
        std::cout << std::endl;
    }

    std::cout << "\nMatrix E (A * F):" << std::endl;
    for (int m = 0; m < M; m++) {
        for (int j = 0; j < J; j++) {
            std::cout << std::setprecision(2) << E[m][j] << " ";
        }
        std::cout << std::endl;
    }

    std::cout << "------------------------------RESULTS END-------------------------" << std::endl;
    return 0;
}
