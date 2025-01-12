/*--------------------Out-Lab Assignemnt 1--------------------*/
/*
 * Created by Hasibul Hossain Rasheeq on Sunday, 1/12/25.
*/

#include <iostream>
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

    std::cout << "Please enter a real value (k): " << std::endl;
    std::cin >> k;
    std::cout << "Please enter the dimensions (M, N, J) for the matrices: " << std::endl;
    std::cin >> M >> N >> J;

    if (M <= 0 || N <= 0 || J <= 0) {
        throw std::invalid_argument("Invalid dimensions! Dimensions have to be bigger than 0");
    } else {
        std:: cout << "There are no errors in the input! Give me a few moments to perform the Matrix operations!" << std::endl;
    }

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
            if (n > m) {
                B[m][n] = 0.25;
            } else {
                B[m][n] = 0.75;
            }
        }
    }

    for (int n = 0; n < N; n++) {
        for (int j = 0; j < J; j++) {
            F[n][j] = 1.0/(n+j);
        }
    }

    for (int m = 0; m < M; m++) {
        for (int n = 0; n < N; n++) {
            C[m][n] = A[m][n] + B[m][n];
        }
    }

    for (int m = 0; m < M; m++) {
        for (int n = 0; n < N; n++) {
            D[m][n] = k * A[m][n];
        }
    }

    for (int m = 0; m < M; m++) {
        for (int j = 0; j < J; j++) {
            double sum = 0;
            for (int n = 0; n < N; n++) {
                sum += A[m][n] * F[n][j];
            }
            E[m][j] = sum;
        }
    }

    std::cout << "------------------------------RESULTS-------------------------" << std::endl;

    for (int m = 0; m < M; m++) {
        for (int n = 0; n < N; n++) {
            std::cout << C[m][n] << " ";
        }
        std::cout << std::endl;
    }

    for (int m = 0; m < M; m++) {
        for (int n = 0; n < N; n++) {
            std::cout << D[m][n] << " ";
        }
        std::cout << std::endl;
    }

    for (int m = 0; m < M; m++) {
        for (int j = 0; j < J; j++) {
            std::cout << E[m][j] << " ";
        }
        std::cout << std::endl;
    }

    std::cout << "------------------------------RESULTS END-------------------------" << std::endl;
    return 0;
}
