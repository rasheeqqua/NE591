//
// Created by Hasibul H. Rasheeq on 1/30/25.
//

// Adaptive Simpson's Method
double adaptiveSimpson(double a, double b, double epsilon, double S, double fa, double fb, double fc, int depth) {
    double c = (a + b) / 2.0;
    double h = b - a;
    double d = (a + c) / 2.0;
    double e = (c + b) / 2.0;
    double fd = integrand(d);
    double fe = integrand(e);
    double Sleft = (h / 12.0) * (fa + 4 * fd + fc);
    double Sright = (h / 12.0) * (fc + 4 * fe + fb);
    double S2 = Sleft + Sright;
    if (depth <= 0 || fabs(S2 - S) <= 15 * epsilon) {
        return S2 + (S2 - S) / 15.0;
    }
    return adaptiveSimpson(a, c, epsilon / 2.0, Sleft, fa, fc, fd, depth - 1) +
           adaptiveSimpson(c, b, epsilon / 2.0, Sright, fc, fb, fe, depth - 1);
}

double exactIntegral(double a, double b, double epsilon = 1E-9, int maxRecursionDepth = 20) {
    double fa = integrand(a);
    double fb = integrand(b);
    double c = (a + b) / 2.0;
    double fc = integrand(c);
    double S = (b - a) * (fa + 4 * fc + fb) / 6.0;
    return adaptiveSimpson(a, b, epsilon, S, fa, fb, fc, maxRecursionDepth);
}
