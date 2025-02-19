//
// Created by Hasibul H. Rasheeq on 02/20/25.
//

bool gaussSeidel(const std::vector<std::vector<double>>& A,
                const std::vector<double>& b,
                std::vector<double>& x,
                double epsilon,
                int maxIter,
                int& actualIter,
                double& finalError) {
    int n = A.size();
    std::vector<double> x_old(n);

    // Main iteration loop
    for (actualIter = 0; actualIter < maxIter; actualIter++) {
        // Store previous iteration values for error calculation
        x_old = x;

        // Update each component of x
        for (int i = 0; i < n; i++) {
            double sum1 = 0.0;  // Sum for already updated components
            double sum2 = 0.0;  // Sum for components not yet updated

            // Use already updated values from current iteration
            for (int j = 0; j < i; j++) {
                sum1 += A[i][j] * x[j];
            }

            // Use values from previous iteration for remaining components
            for (int j = i + 1; j < n; j++) {
                sum2 += A[i][j] * x_old[j];
            }

            // Update current component
            x[i] = (b[i] - sum1 - sum2) / A[i][i];
        }

        // Calculate error for convergence check
        finalError = calculateError(x, x_old);
        if (finalError < epsilon) {
            return true;  // Converged
        }
    }

    return false;  // Failed to converge within maxIter
}