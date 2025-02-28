//
// Steady State One-Speed Diffusion Equation Solver - Milestone 2
// Author: Hasibul H. Rasheeq
// Date: February 27, 2025
// Version: 2.0
//
// This program solves the steady-state, one-speed diffusion equation
// in a 2D rectangular region with vacuum boundary conditions using
// either direct (LUP) or iterative methods (Jacobi, Gauss-Seidel, SOR).
//

#include <iostream>
#include <fstream>
#include <iomanip>
#include <vector>
#include <string>
#include <cmath>
#include <chrono>
#include <algorithm>

// LUP solver from Milestone 1 (only used for comparison and verification)
#include "LUP/LUPFactorization.cpp"
#include "LUP/substitution.cpp"
#include "LUP/ApplyPermutationMatrix.cpp"

class DiffusionSolver {
private:
    // Problem parameters
    int flag;               // Solution method flag
    double a, b;            // Rectangle dimensions
    int m, n;               // Grid dimensions
    double D;               // Diffusion coefficient
    double sigma_a;         // Macroscopic removal cross section
    std::vector<std::vector<double>> q;  // Source term
    double delta, gamma;    // Grid spacing

    // Iteration parameters
    int maxIterations;      // Maximum number of iterations
    double tolerance;       // Convergence tolerance
    double omega;           // Relaxation parameter for SOR

    // Helper function to convert 2D indices to 1D index
    int idx(int i, int j) const {
        return (i-1) * n + (j-1);
    }

    // Matrix-free implementation of diffusion equation
    // This computes Ax for the diffusion operator A and flux vector x
    void applyDiffusionOperator(const std::vector<std::vector<double>>& phi,
                               std::vector<std::vector<double>>& result) const {
        for (int i = 1; i <= m; ++i) {
            for (int j = 1; j <= n; ++j) {
                // Apply the diffusion operator at point (i,j)
                double center = phi[i][j];
                double left = (i > 1) ? phi[i-1][j] : 0.0;     // Account for boundary
                double right = (i < m) ? phi[i+1][j] : 0.0;    // Account for boundary
                double bottom = (j > 1) ? phi[i][j-1] : 0.0;   // Account for boundary
                double top = (j < n) ? phi[i][j+1] : 0.0;      // Account for boundary

                // Calculate the result of diffusion operator
                result[i][j] = -D * ((left - 2*center + right)/(delta*delta) +
                                    (bottom - 2*center + top)/(gamma*gamma)) +
                                    sigma_a * center;
            }
        }
    }

    // Calculate error between two consecutive iterations
    double calculateError(const std::vector<std::vector<double>>& current,
                        const std::vector<std::vector<double>>& previous) const {
        double maxError = 0.0;
        double eps = 1e-15; // To avoid division by zero

        for (int i = 1; i <= m; ++i) {
            for (int j = 1; j <= n; ++j) {
                double absVal = std::abs(current[i][j]);
                if (absVal > eps) {
                    double relError = std::abs(current[i][j] - previous[i][j]) / absVal;
                    maxError = std::max(maxError, relError);
                } else {
                    double absError = std::abs(current[i][j] - previous[i][j]);
                    maxError = std::max(maxError, absError);
                }
            }
        }

        return maxError;
    }

    // Calculate the maximum residual: ||Ax - b||_∞
    double calculateMaxResidual(const std::vector<std::vector<double>>& phi) const {
        std::vector<std::vector<double>> Ax(m+2, std::vector<double>(n+2, 0.0));
        applyDiffusionOperator(phi, Ax);

        double maxResidual = 0.0;
        for (int i = 1; i <= m; ++i) {
            for (int j = 1; j <= n; ++j) {
                double residual = std::abs(Ax[i][j] - q[i-1][j-1]);
                maxResidual = std::max(maxResidual, residual);
            }
        }

        return maxResidual;
    }

    // Matrix-free implementation of Point Jacobi method
    bool solvePointJacobi(std::vector<std::vector<double>>& phi,
                         int& iterations, double& finalError) {
        // Create arrays for old and new solutions
        std::vector<std::vector<double>> phi_old(m+2, std::vector<double>(n+2, 0.0));

        // Diagonal term coefficient (constant for homogeneous media)
        double diagCoef = 2.0*D/(delta*delta) + 2.0*D/(gamma*gamma) + sigma_a;

        // Iterate until convergence or max iterations
        for (iterations = 0; iterations < maxIterations; ++iterations) {
            // Save current solution for error calculation
            phi_old = phi;

            // Update each interior point
            for (int i = 1; i <= m; ++i) {
                for (int j = 1; j <= n; ++j) {
                    // Gather contributions from neighboring points
                    double left = (i > 1) ? phi_old[i-1][j] : 0.0;
                    double right = (i < m) ? phi_old[i+1][j] : 0.0;
                    double bottom = (j > 1) ? phi_old[i][j-1] : 0.0;
                    double top = (j < n) ? phi_old[i][j+1] : 0.0;

                    // Apply Jacobi update formula
                    double numerator = q[i-1][j-1] + D*(left + right)/(delta*delta) +
                                      D*(bottom + top)/(gamma*gamma);
                    phi[i][j] = numerator / diagCoef;
                }
            }

            // Check for convergence
            finalError = calculateError(phi, phi_old);
            if (finalError < tolerance) {
                return true;
            }
        }

        return false; // Did not converge within max iterations
    }

    // Matrix-free implementation of Gauss-Seidel method
    bool solveGaussSeidel(std::vector<std::vector<double>>& phi,
                         int& iterations, double& finalError) {
        // Create array for old solution
        std::vector<std::vector<double>> phi_old(m+2, std::vector<double>(n+2, 0.0));

        // Diagonal term coefficient (constant for homogeneous media)
        double diagCoef = 2.0*D/(delta*delta) + 2.0*D/(gamma*gamma) + sigma_a;

        // Iterate until convergence or max iterations
        for (iterations = 0; iterations < maxIterations; ++iterations) {
            // Save current solution for error calculation
            phi_old = phi;

            // Update each interior point
            for (int i = 1; i <= m; ++i) {
                for (int j = 1; j <= n; ++j) {
                    // Gather contributions from neighboring points
                    // Use updated values where available
                    double left = (i > 1) ? phi[i-1][j] : 0.0;        // Already updated
                    double right = (i < m) ? phi_old[i+1][j] : 0.0;   // Not yet updated

                    // For this row, points to the left are updated
                    double bottom = (j > 1) ? phi[i][j-1] : 0.0;      // Already updated
                    double top = (j < n) ? phi_old[i][j+1] : 0.0;     // Not yet updated

                    // Apply Gauss-Seidel update formula
                    double numerator = q[i-1][j-1] + D*(left + right)/(delta*delta) +
                                      D*(bottom + top)/(gamma*gamma);
                    phi[i][j] = numerator / diagCoef;
                }
            }

            // Check for convergence
            finalError = calculateError(phi, phi_old);
            if (finalError < tolerance) {
                return true;
            }
        }

        return false; // Did not converge within max iterations
    }

    // Matrix-free implementation of SOR method
    bool solveSOR(std::vector<std::vector<double>>& phi,
                 int& iterations, double& finalError) {
        // Create array for old solution
        std::vector<std::vector<double>> phi_old(m+2, std::vector<double>(n+2, 0.0));

        // Diagonal term coefficient (constant for homogeneous media)
        double diagCoef = 2.0*D/(delta*delta) + 2.0*D/(gamma*gamma) + sigma_a;

        // Iterate until convergence or max iterations
        for (iterations = 0; iterations < maxIterations; ++iterations) {
            // Save current solution for error calculation
            phi_old = phi;

            // Update each interior point
            for (int i = 1; i <= m; ++i) {
                for (int j = 1; j <= n; ++j) {
                    // Gather contributions from neighboring points
                    // Use updated values where available
                    double left = (i > 1) ? phi[i-1][j] : 0.0;        // Already updated
                    double right = (i < m) ? phi_old[i+1][j] : 0.0;   // Not yet updated

                    // For this row, points to the left are updated
                    double bottom = (j > 1) ? phi[i][j-1] : 0.0;      // Already updated
                    double top = (j < n) ? phi_old[i][j+1] : 0.0;     // Not yet updated

                    // Calculate Gauss-Seidel value
                    double numerator = q[i-1][j-1] + D*(left + right)/(delta*delta) +
                                      D*(bottom + top)/(gamma*gamma);
                    double phi_gs = numerator / diagCoef;

                    // Apply SOR formula: phi_new = (1-ω)*phi_old + ω*phi_gs
                    phi[i][j] = (1.0 - omega) * phi_old[i][j] + omega * phi_gs;
                }
            }

            // Check for convergence
            finalError = calculateError(phi, phi_old);
            if (finalError < tolerance) {
                return true;
            }
        }

        return false; // Did not converge within max iterations
    }

    // Direct solver using LUP decomposition from Milestone 1
    // Only used for comparison/verification
    std::vector<std::vector<double>> solveLUP() const {
        // Build coefficient matrix and RHS vector
        int size = m * n;
        std::vector<std::vector<double>> A(size, std::vector<double>(size, 0.0));
        std::vector<double> b(size, 0.0);

        // Fill matrix and vector
        for (int i = 1; i <= m; ++i) {
            for (int j = 1; j <= n; ++j) {
                int row = idx(i, j);

                // Diagonal term
                double diag = 2.0*D/(delta*delta) + 2.0*D/(gamma*gamma) + sigma_a;
                A[row][row] = diag;

                // Off-diagonal terms
                if (i > 1) A[row][idx(i-1, j)] = -D/(delta*delta);
                if (i < m) A[row][idx(i+1, j)] = -D/(delta*delta);
                if (j > 1) A[row][idx(i, j-1)] = -D/(gamma*gamma);
                if (j < n) A[row][idx(i, j+1)] = -D/(gamma*gamma);

                // RHS vector
                b[row] = q[i-1][j-1];
            }
        }

        // Initialize matrices for LUP factorization
        std::vector<std::vector<double>> L(size, std::vector<double>(size, 0.0));
        std::vector<std::vector<double>> U(size, std::vector<double>(size, 0.0));
        std::vector<std::vector<double>> P(size, std::vector<double>(size, 0.0));

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
        std::vector<std::vector<double>> phi(m+2, std::vector<double>(n+2, 0.0));
        for (int i = 1; i <= m; ++i) {
            for (int j = 1; j <= n; ++j) {
                phi[i][j] = x[idx(i, j)];
            }
        }

        return phi;
    }

public:
    // Function to read input file
    bool readInput(const std::string& filename) {
        std::ifstream input(filename);
        if (!input.is_open()) {
            std::cerr << "Error: Cannot open input file " << filename << std::endl;
            return false;
        }

        // Read solution method flag
        input >> flag;

        // Read method-specific parameters
        switch (flag) {
            case 0: // LUP (direct)
                input >> maxIterations >> tolerance >> omega;
                // LUP doesn't need additional parameters
                maxIterations = 0;
                tolerance = 0.0;
                omega = 0.0;
                break;

            case 1: // Point Jacobi
            case 2: // Gauss-Seidel
                // Read max iterations and tolerance
                input >> maxIterations >> tolerance >> omega;
                omega = 0.0; // Not used
                break;

            case 3: // SOR
                // Read max iterations, tolerance, and relaxation parameter
                input >> maxIterations >> tolerance >> omega;
                // Validate omega
                if (omega <= 0.0 || omega >= 2.0) {
                    std::cerr << "Error: SOR relaxation parameter must be in range (0,2)" << std::endl;
                    return false;
                }
                break;

            default:
                std::cerr << "Error: Invalid solution method flag. Must be 0-3." << std::endl;
                return false;
        }

        // Read problem geometry
        input >> a >> b;
        if (a <= 0 || b <= 0) {
            std::cout << flag << maxIterations << tolerance << omega << std::endl;
            std::cout << a << b << std::endl;
            std::cerr << "Error: Invalid rectangle dimensions. Both a and b must be positive." << std::endl;
            return false;
        }

        // Read grid dimensions
        input >> m >> n;
        if (m <= 0 || n <= 0) {
            std::cerr << "Error: Invalid grid dimensions. Both m and n must be positive integers." << std::endl;
            return false;
        }

        // Read physical parameters
        input >> D >> sigma_a;
        if (D <= 0) {
            std::cerr << "Error: Invalid diffusion coefficient. D must be positive." << std::endl;
            return false;
        }
        if (sigma_a < 0) {
            std::cerr << "Error: Invalid removal cross section. Sigma_a must be non-negative." << std::endl;
            return false;
        }

        // Initialize source term matrix
        q.resize(m, std::vector<double>(n));
        for (int i = 0; i < m; ++i) {
            for (int j = 0; j < n; ++j) {
                input >> q[i][j];
                // Check if source term is non-negative
                if (q[i][j] < 0) {
                    std::cerr << "Error: Invalid source term at position (" << i+1 << "," << j+1
                              << "). Source values must be non-negative." << std::endl;
                    return false;
                }
            }
        }

        // Calculate grid spacing
        delta = a / (m + 1);
        gamma = b / (n + 1);

        std::cout << "Successfully read input file with grid size " << m << "x" << n << std::endl;
        return true;
    }

    // Add this method to your DiffusionSolver class:

    // Method to set parameters directly (for testing purposes)
    void setParameters(int methodFlag,
                      int maxIter,
                      double tol,
                      double w,
                      double rectA, double rectB,
                      int gridM, int gridN,
                      double diffCoef, double removalXS,
                      const std::vector<std::vector<double>>& sourceTerms) {
        // Set method flag
        flag = methodFlag;

        // Set iteration parameters
        maxIterations = maxIter;
        tolerance = tol;
        omega = w;

        // Set geometry parameters
        a = rectA;
        b = rectB;
        m = gridM;
        n = gridN;

        // Set physical parameters
        D = diffCoef;
        sigma_a = removalXS;

        // Set source terms
        q = sourceTerms;

        // Calculate grid spacing
        delta = a / (m + 1);
        gamma = b / (n + 1);
    }

    // Main solve function that dispatches to appropriate method
    std::vector<std::vector<double>> solve(int& iterations, double& finalError, bool& converged) {
        // Initialize solution with zeros (vacuum boundary conditions)
        std::vector<std::vector<double>> phi(m+2, std::vector<double>(n+2, 0.0));

        // Select solution method based on flag
        switch (flag) {
            case 0: // LUP (direct)
                try {
                    phi = solveLUP();
                    iterations = 0;
                    finalError = 0.0;
                    converged = true;
                } catch (const std::exception& e) {
                    std::cerr << "Error in LUP solver: " << e.what() << std::endl;
                    converged = false;
                }
                break;

            case 1: // Point Jacobi
                converged = solvePointJacobi(phi, iterations, finalError);
                break;

            case 2: // Gauss-Seidel
                converged = solveGaussSeidel(phi, iterations, finalError);
                break;

            case 3: // SOR
                converged = solveSOR(phi, iterations, finalError);
                break;

            default:
                std::cerr << "Error: Invalid solution method" << std::endl;
                converged = false;
        }

        return phi;
    }

    // Function to write output to file
    void writeOutput(const std::vector<std::vector<double>>& phi,
                    const std::string& filename,
                    int iterations,
                    double finalError,
                    bool converged,
                    double executionTime) {
        std::ofstream output(filename);
        if (!output.is_open()) {
            std::cerr << "Error: Unable to open output file: " << filename << std::endl;
            return;
        }

        // Get current date/time
        std::time_t now = std::time(nullptr);
        char timeStr[100];
        std::strftime(timeStr, sizeof(timeStr), "%b %d %Y %H:%M:%S", std::localtime(&now));

        output << "Steady State One-Speed Diffusion Equation Solver - Milestone 2" << std::endl;
        output << "Version: 2.0" << std::endl;
        output << "Author: Hasibul H. Rasheeq" << std::endl;
        output << "Date: " << timeStr << std::endl;
        output << "----------------------------------------" << std::endl << std::endl;

        // Method used
        output << "Solution Method: ";
        switch (flag) {
            case 0: output << "LUP Decomposition (Direct)"; break;
            case 1: output << "Point Jacobi Iteration"; break;
            case 2: output << "Gauss-Seidel Iteration"; break;
            case 3: output << "Successive Over-Relaxation (SOR) with ω = " << omega; break;
        }
        output << std::endl << std::endl;

        // Problem parameters
        output << "Problem Parameters:" << std::endl;
        output << "  Grid dimensions: " << m << " x " << n << " (" << m*n << " unknowns)" << std::endl;
        output << "  Rectangle dimensions: " << a << " cm x " << b << " cm" << std::endl;
        output << "  Grid spacing: δ = " << delta << " cm, γ = " << gamma << " cm" << std::endl;
        output << "  Physical parameters: D = " << D << " cm, Σa = " << sigma_a << " cm-1" << std::endl << std::endl;

        // Iteration information (for iterative methods)
        if (flag > 0) {
            output << "Iteration Information:" << std::endl;
            output << "  Maximum iterations allowed: " << maxIterations << std::endl;
            output << "  Convergence tolerance: " << std::scientific << std::setprecision(6) << tolerance << std::endl;
            if (converged) {
                output << "  Converged after " << iterations << " iterations" << std::endl;
                output << "  Final error: " << std::scientific << std::setprecision(6) << finalError << std::endl;
            } else {
                output << "  Failed to converge within " << maxIterations << " iterations" << std::endl;
                output << "  Final error: " << std::scientific << std::setprecision(6) << finalError << std::endl;
            }
            output << std::endl;
        }

        // Performance metrics
        output << "Performance Metrics:" << std::endl;
        output << "  Execution time: " << std::fixed << std::setprecision(8) << executionTime << " seconds" << std::endl;

        // Calculate and print maximum residual
        double maxResidual = calculateMaxResidual(phi);
        output << "  Maximum absolute residual: " << std::scientific << std::setprecision(6)
               << maxResidual << std::endl << std::endl;

        // Solution
        output << "Scalar Flux Solution:" << std::endl;
        output << std::scientific << std::setprecision(6);
        for (int i = 0; i <= m + 1; ++i) {
            for (int j = 0; j <= n + 1; ++j) {
                output << std::setw(15) << phi[i][j];
            }
            output << "\n";
        }

        output.close();
    }

    // Compare two solutions (e.g., LUP vs iterative)
    double compareSolutions(const std::vector<std::vector<double>>& phi1,
                           const std::vector<std::vector<double>>& phi2) {
        double maxDiff = 0.0;
        double eps = 1e-15;

        for (int i = 1; i <= m; ++i) {
            for (int j = 1; j <= n; ++j) {
                double absVal = std::max(std::abs(phi1[i][j]), std::abs(phi2[i][j]));
                if (absVal > eps) {
                    double relDiff = std::abs(phi1[i][j] - phi2[i][j]) / absVal;
                    maxDiff = std::max(maxDiff, relDiff);
                } else {
                    double absDiff = std::abs(phi1[i][j] - phi2[i][j]);
                    maxDiff = std::max(maxDiff, absDiff);
                }
            }
        }

        return maxDiff;
    }

    // Getter methods for problem parameters
    int getFlag() const { return flag; }
    int getMaxIterations() const { return maxIterations; }
    double getTolerance() const { return tolerance; }
    double getOmega() const { return omega; }
    std::pair<int, int> getGridDimensions() const { return {m, n}; }
};
