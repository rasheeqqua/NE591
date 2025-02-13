//
// Created by hrashee on 2/13/25.
//

#include "LUPFactorization.cpp"
#include "substitution.cpp"
#include "ApplyPermutationMatrix.cpp"

class DiffusionSolver {
public:
    int flag;           // Solution method flag
    double a, b;        // Rectangle dimensions
    int m, n;          // Grid dimensions
    double D;          // Diffusion coefficient
    double sigma_a;    // Macroscopic removal cross section
    std::vector<std::vector<double>> q;  // Source term
    double delta, gamma;  // Grid spacing

    // Helper function to convert 2D indices to 1D index
    int idx(int i, int j) const {
        return i * n + j;
    }

    // Function to build the coefficient matrix
    std::vector<std::vector<double>> buildMatrix() const {
        int size = m * n;
        std::vector<std::vector<double>> A(size, std::vector<double>(size, 0.0));

        for (int i = 1; i <= m; ++i) {
            for (int j = 1; j <= n; ++j) {
                int row = idx(i-1, j-1);

                // Diagonal term
                double diag = 2.0*D/(delta*delta) + 2.0*D/(gamma*gamma) + sigma_a;
                A[row][row] = diag;

                // Off-diagonal terms
                if (i > 1) A[row][idx(i-2, j-1)] = -D/(delta*delta);
                if (i < m) A[row][idx(i, j-1)] = -D/(delta*delta);
                if (j > 1) A[row][idx(i-1, j-2)] = -D/(gamma*gamma);
                if (j < n) A[row][idx(i-1, j)] = -D/(gamma*gamma);
            }
        }
        return A;
    }

    // Function to build the right-hand side vector
    std::vector<double> buildRHS() const {
        std::vector<double> b;
        for (int i = 0; i < m; ++i) {
            for (int j = 0; j < n; ++j) {
                b.push_back(q[i][j]);
            }
        }
        return b;
    }

public:
    // Function to read input from file
    bool readInput(const std::string& filename) {
        std::ifstream input(filename);
        if (!input.is_open()) {
            std::cerr << "Error: Cannot open input file " << filename << std::endl;
            return false;
        }

        // Read input parameters
        input >> flag;
        if (flag != 0) {
            std::cerr << "Error: Invalid flag value. Expected 0 for direct solution." << std::endl;
            return false;
        }

        input >> a >> b;
        if (a <= 0 || b <= 0) {
            std::cerr << "Error: Invalid rectangle dimensions." << std::endl;
            return false;
        }

        input >> m >> n;
        if (m <= 0 || n <= 0) {
            std::cerr << "Error: Invalid grid dimensions." << std::endl;
            return false;
        }

        input >> D >> sigma_a;
        if (D <= 0 || sigma_a < 0) {
            std::cerr << "Error: Invalid physical parameters." << std::endl;
            return false;
        }

        // Initialize source term matrix
        q.resize(m, std::vector<double>(n));
        for (int i = 0; i < m; ++i) {
            for (int j = 0; j < n; ++j) {
                input >> q[i][j];
            }
        }

        // Calculate grid spacing
        delta = a / (m + 1);
        gamma = b / (n + 1);

        return true;
    }

    // Function to solve the diffusion equation
    std::vector<std::vector<double>> solve() {
        // Build the coefficient matrix and RHS vector
        auto A = buildMatrix();
        auto b = buildRHS();

        // Initialize matrices for LUP factorization
        int size = m * n;
        std::vector<std::vector<double>> L(size, std::vector<double>(size));
        std::vector<std::vector<double>> U(size, std::vector<double>(size));
        std::vector<std::vector<double>> P(size, std::vector<double>(size));

        // Perform LUP factorization
        if (!lupFactorize(A, L, U, P)) {
            throw std::runtime_error("LUP factorization failed");
        }

        // Solve the system
        std::vector<double> Pb(size);
        applyPermutationMatrix(P, b, Pb);

        std::vector<double> y(size), x(size);
        forwardSubstitution(L, Pb, y);
        backSubstitution(U, y, x);

        // Convert solution back to 2D grid
        std::vector<std::vector<double>> phi(m + 2, std::vector<double>(n + 2, 0.0));
        for (int i = 1; i <= m; ++i) {
            for (int j = 1; j <= n; ++j) {
                phi[i][j] = x[idx(i-1, j-1)];
            }
        }

        return phi;
    }

    // Function to write output to file
    void writeOutput(const std::vector<std::vector<double>>& phi, const std::string& filename) {
        std::ofstream output(filename);
        output << std::scientific << std::setprecision(6);

        output << "Steady State One-Speed Diffusion Equation Solution\n";
        output << "Grid dimensions: " << m << " x " << n << "\n";
        output << "Rectangle dimensions: " << a << " cm x " << b << " cm\n";
        output << "Physical parameters: D = " << D << " cm, sigma_a = " << sigma_a << " cm-1\n\n";

        output << "Scalar Flux Solution:\n";
        for (int i = 0; i <= m + 1; ++i) {
            for (int j = 0; j <= n + 1; ++j) {
                output << std::setw(15) << phi[i][j];
            }
            output << "\n";
        }
    }
};
