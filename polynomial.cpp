#include <vector>

double lagrangeInterpolate(double xEval, const std::vector<double>& x, const std::vector<double>& y) {
    int n = x.size();
    double Lx = 0.0;

    for (int i = 0; i < n; ++i) {
        double term = y[i];
        for (int j = 0; j < n; ++j) {
            if (j != i) {
                term *= (xEval - x[j]) / (x[i] - x[j]);
            }
        }
        Lx += term;
    }

    return Lx;
}