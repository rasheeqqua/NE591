//
// Created by Hasibul H. Rasheeq on 02/20/25.
//

bool sor(const std::vector<std::vector<double>>& A,
         const std::vector<double>& b,
         std::vector<double>& x,
         double omega,
         double epsilon,
         int maxIter,
         int& actualIter,
         double& finalError) {
    int n = A.size();
    std::vector<double> x_old(n);

    // Validate relaxation parameter
    if (omega <= 0.0 || omega >= 2.0) {
        return false;  // Invalid relaxation parameter
    }

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

            // Calculate new value with relaxation
            double x_new = (b[i] - sum1 - sum2) / A[i][i];
            // Apply SOR formula: x_new = (1-ω)x_old + ω*x_new
            x[i] = (1.0 - omega) * x_old[i] + omega * x_new;
        }

        // Calculate error for convergence check
        finalError = calculateError(x, x_old);
        if (finalError < epsilon) {
            return true;  // Converged
        }
    }

    return false;  // Failed to converge within maxIter
}