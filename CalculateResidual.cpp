//
// Created by Hasibul H. Rasheeq on 02/20/25.
//

// Updated calculateError function with denominator check
double calculateError(const std::vector<double>& x_new, const std::vector<double>& x_old) {
    double maxDiff = 0.0;
    const double TOLERANCE = 1e-15;  // Small number to check for near-zero values

    for (size_t i = 0; i < x_new.size(); i++) {
        // Check if denominator is not too close to zero
        if (std::abs(x_new[i]) > TOLERANCE) {
            double diff = std::abs(x_new[i] - x_old[i]) / std::abs(x_new[i]);
            maxDiff = std::max(maxDiff, diff);
        } else {
            // If denominator is too small, use absolute difference instead
            double diff = std::abs(x_new[i] - x_old[i]);
            maxDiff = std::max(maxDiff, diff);
        }
    }
    return maxDiff;
}

// Calculate maximum residual
double calculateResidual(const std::vector<std::vector<double>>& A,
                        const std::vector<double>& x,
                        const std::vector<double>& b) {
    double maxResidual = 0.0;
    int n = A.size();
    for (int i = 0; i < n; i++) {
        double sum = 0.0;
        for (int j = 0; j < n; j++) {
            sum += A[i][j] * x[j];
        }
        maxResidual = std::max(maxResidual, std::abs(sum - b[i]));
    }
    return maxResidual;
}
