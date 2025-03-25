#include <vector>

void forwardSubstitution(const std::vector<std::vector<double>>& L, const std::vector<double>& b, std::vector<double>& y) {
    int n = L.size();
    for (int i = 0; i < n; ++i) {
        y[i] = b[i];
        for (int j = 0; j < i; ++j) {
            y[i] -= L[i][j] * y[j];
        }
        // Assuming L[i][i] == 1, so division is not necessary
    }
}

void backSubstitution(const std::vector<std::vector<double>>& U, const std::vector<double>& y, std::vector<double>& x) {
    int n = U.size();
    for (int i = n - 1; i >= 0; --i) {
        x[i] = y[i];
        for (int j = i + 1; j < n; ++j) {
            x[i] -= U[i][j] * x[j];
        }
        x[i] /= U[i][i];
    }
}
