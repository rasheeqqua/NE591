//
// Created by Hasibul H. Rasheeq on 1/30/25.
//
#include <vector>

void gaussLegendreNodesAndWeights(int n, std::vector<double>& nodes, std::vector<double>& weights) {
    nodes.clear();
    weights.clear();

    if (n == 2) {
        nodes = { -0.5773502692, 0.5773502692 };
        weights = { 1.0, 1.0 };
    } else if (n == 3) {
        nodes = { -0.7745966692, 0.0, 0.7745966692 };
        weights = { 0.5555555556, 0.8888888889, 0.5555555556 };
    } else if (n == 4) {
        nodes = { -0.8611363116, -0.3399810436, 0.3399810436, 0.8611363116 };
        weights = { 0.3478548451, 0.6521451549, 0.6521451549, 0.3478548451 };
    } else if (n == 5) {
        nodes = { -0.9061798459, -0.5384693101, 0.0, 0.5384693101, 0.9061798459 };
        weights = { 0.2369268851, 0.4786286705, 0.5688888889, 0.4786286705, 0.2369268851 };
    } else {
        std::cout << "Gauss-Legendre Quadrature implemented for n = 2 to 5 only.\n";
        std::cout << "Using n = 2 by default.\n";
        nodes = { -0.5773502692, 0.5773502692 };
        weights = { 1.0, 1.0 };
    }
}

double gaussLegendreQuadrature(double a, double b, int n) {
    std::vector<double> nodes, weights;
    gaussLegendreNodesAndWeights(n, nodes, weights);

    // Print nodes and weights to the terminal
    std::cout << "Quadrature Points (Nodes) and Weights for n = " << n << ":\n";
    for (int i = 0; i < nodes.size(); ++i) {
        std::cout << "Node " << i+1 << ": " << nodes[i] << ", Weight: " << weights[i] << "\n";
    }

    // Apply the Gauss-Legendre Quadrature formula
    double sum = 0.0;
    for (int i = 0; i < nodes.size(); ++i) {
        // Change variable from [-1,1] to [a,b]
        double xi = ((b - a) / 2.0) * nodes[i] + (a + b) / 2.0;
        sum += weights[i] * integrand(xi);
    }
    return (b - a) / 2.0 * sum;
}