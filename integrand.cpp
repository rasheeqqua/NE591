//
// Created by Hasibul H. Rasheeq on 1/30/25.
//
#include <cmath>

// Task 4
double integrand(double x) {
    if (x < 0.25) {
        return exp(-x + 0.25); // For x in [-1, 0.25)
    } else {
        return exp(x - 0.25);  // For x in [0.25, 1]
    }
}