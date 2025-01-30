//
// Created by Hasibul H. Rasheeq on 1/30/25.
//
double compositeTrapezoidal(double a, double b, int m) {
    double h = (b - a) / m;
    double sum = 0.5 * (integrand(a) + integrand(b));
    for (int i = 1; i < m; ++i) {
        double x = a + i * h;
        sum += integrand(x);
    }
    return h * sum;
}