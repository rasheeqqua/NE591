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

#include "DiffusionSolver.cpp"

int main() {
    // Initialize solver
    DiffusionSolver solver;

    // Define the input and output files
    std::string inputFile = "input.txt";
    std::string outputFile = "output.txt";

    // Read input file
    if (!solver.readInput(inputFile)) {
        std::cerr << "Error reading input file. Exiting." << std::endl;
        return 1;
    }

    // Variables to store results
    int iterations;
    double finalError;
    bool converged;

    // Measure execution time
    auto startTime = std::chrono::high_resolution_clock::now();

    // Solve the system using the selected method
    auto solution = solver.solve(iterations, finalError, converged);

    // End timing
    auto endTime = std::chrono::high_resolution_clock::now();
    double executionTime = std::chrono::duration<double>(endTime - startTime).count();

    // Write output
    solver.writeOutput(solution, outputFile, iterations, finalError, converged, executionTime);

    // Optional: If using an iterative method, compare with LUP solution for verification
    if (solver.getFlag() > 0 && converged && solver.getGridDimensions().first * solver.getGridDimensions().second <= 400) {
        std::cout << "Comparing iterative solution with LUP solution for verification..." << std::endl;

        // Save the original flag
        int originalFlag = solver.getFlag();

        // Create a new solver instance for LUP
        DiffusionSolver lupSolver;
        if (lupSolver.readInput(inputFile)) {
            int lupIterations;
            double lupError;
            bool lupConverged;

            // Force LUP method (flag = 0)
            auto lupSolution = lupSolver.solve(lupIterations, lupError, lupConverged);

            if (lupConverged) {
                // Compare solutions
                double maxDifference = solver.compareSolutions(solution, lupSolution);
                std::cout << "Maximum relative difference between iterative and LUP solutions: "
                          << std::scientific << std::setprecision(6) << maxDifference << std::endl;

                // Write LUP solution to a comparison file
                lupSolver.writeOutput(lupSolution, "lup_solution.txt", 0, 0.0, true, 0.0);
            }
        }
    }

    std::cout << "Solution completed. Results written to " << outputFile << std::endl;

    return 0;
}