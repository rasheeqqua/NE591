//
// Created by Hasibul H. Rasheeq on 1/30/25.
//
#include "integrand.cpp"
double compositeSimpson(double a, double b, int m) {
    double h = (b - a) / m;
    double sum = integrand(a) + integrand(b);
    for (int i = 1; i < m; ++i) {
        double x = a + i * h;
        if (i % 2 == 0) {
            sum += 2 * integrand(x);
        } else {
            sum += 4 * integrand(x);
        }
    }
    return (h / 3) * sum;
}