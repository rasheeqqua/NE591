//
// Created by Hasibul H. Rasheeq on 02/08/25.
//

void applyPermutationMatrix(const std::vector<std::vector<double>>& P,
                             const std::vector<double>& b,
                             std::vector<double>& Pb) {
    int n = P.size();

    for (int i = 0; i < n; ++i) {
        // Find which row this row should be swapped with
        for (int j = 0; j < n; ++j) {
            if (P[i][j] == 1) {
                Pb[i] = b[j];
                break;
            }
        }
    }
}
