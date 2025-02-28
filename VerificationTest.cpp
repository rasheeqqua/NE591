// Diffusion Equation Solver Verification Tests
// Version: 1.0
// Author: Hasibul H. Rasheeq
// Date: February 27, 2025
//
// This program implements verification tests for the diffusion equation solver
// to validate its correctness through symmetry properties and known solutions.

#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>
#include <vector>
#include <cmath>
#include <algorithm>
#include "DiffusionClass.cpp"

// Function to check if a solution has reflective symmetry across a specified axis
bool verifyReflectiveSymmetry(const std::vector<std::vector<double>>& solution,
                             int m, int n, bool checkX = true, bool checkY = true) {
    double epsilon = 1e-10; // Tolerance for floating-point comparisons

    // Check symmetry across the x-axis (horizontal)
    if (checkY) {
        int midY = (n + 1) / 2; // Middle row index
        for (int i = 0; i <= m + 1; ++i) {
            for (int j = 0; j <= midY; ++j) {
                int mirrorJ = n + 1 - j;
                if (std::abs(solution[i][j] - solution[i][mirrorJ]) > epsilon) {
                    std::cout << "Failed Y-symmetry at (" << i << "," << j << ") and ("
                             << i << "," << mirrorJ << "): "
                             << solution[i][j] << " != " << solution[i][mirrorJ] << std::endl;
                    return false;
                }
            }
        }
    }

    // Check symmetry across the y-axis (vertical)
    if (checkX) {
        int midX = (m + 1) / 2; // Middle column index
        for (int i = 0; i <= midX; ++i) {
            for (int j = 0; j <= n + 1; ++j) {
                int mirrorI = m + 1 - i;
                if (std::abs(solution[i][j] - solution[mirrorI][j]) > epsilon) {
                    std::cout << "Failed X-symmetry at (" << i << "," << j << ") and ("
                             << mirrorI << "," << j << "): "
                             << solution[i][j] << " != " << solution[mirrorI][j] << std::endl;
                    return false;
                }
            }
        }
    }

    return true;
}

// Function to check if a solution has rotational symmetry around the center
bool verifyRotationalSymmetry(const std::vector<std::vector<double>>& solution, int m, int n) {
    double epsilon = 1e-10; // Tolerance for floating-point comparisons

    // Check if grid is square (required for rotational symmetry)
    if (m != n) {
        std::cout << "Grid must be square to check rotational symmetry" << std::endl;
        return false;
    }

    // Center point of the grid
    int centerI = (m + 1) / 2;
    int centerJ = (n + 1) / 2;

    // Check 90-degree rotational symmetry
    for (int i = 0; i <= m + 1; ++i) {
        for (int j = 0; j <= n + 1; ++j) {
            // Calculate distance from center
            int di = i - centerI;
            int dj = j - centerJ;

            // 90-degree rotation: (x,y) -> (-y,x)
            int newI = centerI - dj;
            int newJ = centerJ + di;

            // Check if the rotated point is within bounds
            if (newI >= 0 && newI <= m + 1 && newJ >= 0 && newJ <= n + 1) {
                if (std::abs(solution[i][j] - solution[newI][newJ]) > epsilon) {
                    std::cout << "Failed 90° rotational symmetry at (" << i << "," << j << ") and ("
                             << newI << "," << newJ << "): "
                             << solution[i][j] << " != " << solution[newI][newJ] << std::endl;
                    return false;
                }
            }
        }
    }

    // Check 180-degree rotational symmetry
    for (int i = 0; i <= m + 1; ++i) {
        for (int j = 0; j <= n + 1; ++j) {
            // 180-degree rotation: (x,y) -> (-x,-y)
            int newI = 2 * centerI - i;
            int newJ = 2 * centerJ - j;

            // Check if the rotated point is within bounds
            if (newI >= 0 && newI <= m + 1 && newJ >= 0 && newJ <= n + 1) {
                if (std::abs(solution[i][j] - solution[newI][newJ]) > epsilon) {
                    std::cout << "Failed 180° rotational symmetry at (" << i << "," << j << ") and ("
                             << newI << "," << newJ << "): "
                             << solution[i][j] << " != " << solution[newI][newJ] << std::endl;
                    return false;
                }
            }
        }
    }

    return true;
}

// Function to create a test case with a symmetric source distribution
void createSymmetricTestCase(const std::string& filename, int m, int n,
                           double a, double b, double D, double sigma_a,
                           const std::string& symmetryType) {
    std::ofstream file(filename);

    file << "0\n";               // flag
    file << a << " " << b << "\n"; // rectangle dimensions
    file << m << " " << n << "\n"; // grid dimensions
    file << D << " " << sigma_a << "\n"; // physical parameters

    // Create source distribution based on desired symmetry type
    if (symmetryType == "reflective_x") {
        // Source symmetric about x-axis (horizontally symmetric)
        int midY = n / 2;
        for (int i = 0; i < m; ++i) {
            for (int j = 0; j < n; ++j) {
                if (j == midY) {
                    file << "10.0 "; // Source along middle row
                } else {
                    file << "0.0 ";
                }
            }
            file << "\n";
        }
    }
    else if (symmetryType == "reflective_y") {
        // Source symmetric about y-axis (vertically symmetric)
        int midX = m / 2;
        for (int i = 0; i < m; ++i) {
            for (int j = 0; j < n; ++j) {
                if (i == midX) {
                    file << "10.0 "; // Source along middle column
                } else {
                    file << "0.0 ";
                }
            }
            file << "\n";
        }
    }
    else if (symmetryType == "reflective_xy") {
        // Source symmetric about both x and y axes
        for (int i = 0; i < m; ++i) {
            for (int j = 0; j < n; ++j) {
                file << "10.0 "; // Uniform source
            }
            file << "\n";
        }
    }
    else if (symmetryType == "rotational") {
        // Source with rotational symmetry
        int midX = m / 2;
        int midY = n / 2;

        // For rotational symmetry, we need a square grid
        if (m != n) {
            std::cerr << "Error: Rotational symmetry requires square grid (m=n)" << std::endl;
            file.close();
            return;
        }

        for (int i = 0; i < m; ++i) {
            for (int j = 0; j < n; ++j) {
                // Place sources at positions symmetric under 90° rotation
                if ((i == midX-1 && j == midY-1) ||
                    (i == midX-1 && j == midY+1) ||
                    (i == midX+1 && j == midY-1) ||
                    (i == midX+1 && j == midY+1)) {
                    file << "10.0 ";
                } else {
                    file << "0.0 ";
                }
            }
            file << "\n";
        }
    }

    file.close();
}

// Function to document and analyze the verification test results
void documentTestResults(const std::string& testName,
                       const std::string& testDescription,
                       const std::string& hypothesis,
                       bool passed,
                       const std::string& outputFile) {
    std::ofstream doc(testName + "_verification.md");

    doc << "# Verification Test: " << testName << "\n\n";
    doc << "## Test Description\n";
    doc << testDescription << "\n\n";

    doc << "## Hypothesis\n";
    doc << hypothesis << "\n\n";

    doc << "## Test Results\n";
    if (passed) {
        doc << "**PASSED** ✅\n\n";
        doc << "The test successfully verified the hypothesis. The numerical solution correctly "
            << "exhibits the expected symmetry properties.\n\n";
    } else {
        doc << "**FAILED** ❌\n\n";
        doc << "The test failed to verify the hypothesis. The numerical solution does not "
            << "exhibit the expected symmetry properties. Further investigation is needed.\n\n";
    }

    doc << "## Input/Output Files\n";
    doc << "- Input file: `" << testName << ".txt`\n";
    doc << "- Output file: `" << outputFile << "`\n\n";

    doc << "## Conclusion\n";
    if (passed) {
        doc << "This verification test confirms that the diffusion equation solver correctly "
            << "preserves the symmetry of the physical problem. This increases confidence in "
            << "the correctness of the implementation.\n";
    } else {
        doc << "This verification test reveals potential issues in the diffusion equation solver "
            << "implementation that need to be addressed. The solver should preserve the symmetry "
            << "of the physical problem, but the numerical results do not exhibit the expected symmetry.\n";
    }

    doc.close();
}

int main() {
    std::cout << "Running Diffusion Equation Solver Verification Tests\n";
    std::cout << "===================================================\n\n";

    // Test parameters
    double D = 0.142;      // Diffusion coefficient
    double sigma_a = 0.0222; // Macroscopic removal cross section

    //--------------------------------------------------------------------------
    // Test 1: Reflective Symmetry (X and Y axes)
    //--------------------------------------------------------------------------
    std::cout << "Test 1: Reflective Symmetry (X and Y axes)\n";
    std::cout << "----------------------------------------\n";

    // Create test case with uniform source (symmetric in X and Y)
    int m1 = 5, n1 = 5;
    double a1 = 60.0, b1 = 60.0;
    std::string test1 = "reflective_xy_symmetry";
    std::string output1 = test1 + "_output.txt";

    createSymmetricTestCase(test1 + ".txt", m1, n1, a1, b1, D, sigma_a, "reflective_xy");

    // Solve the problem
    DiffusionSolver solver1;
    solver1.readInput(test1 + ".txt");
    auto solution1 = solver1.solve();
    solver1.writeOutputWithDetails(solution1, output1, 0.0, 0, {}, 0.0);

    // Verify reflective symmetry
    bool test1Passed = verifyReflectiveSymmetry(solution1, m1, n1, true, true);
    std::cout << "Test 1 " << (test1Passed ? "PASSED" : "FAILED") << "\n\n";

    // Document the test
    std::string desc1 = "This test verifies that the solver correctly preserves reflective symmetry "
                      "in both X and Y directions when the source distribution is uniform across "
                      "the domain. A uniform source with value 10.0 is used throughout the domain. "
                      "The physical parameters are D = " + std::to_string(D) + " and Σₐ = " +
                      std::to_string(sigma_a) + ". The rectangular domain is " +
                      std::to_string(a1) + " cm × " + std::to_string(b1) + " cm with a " +
                      std::to_string(m1) + "×" + std::to_string(n1) + " grid.";

    std::string hyp1 = "Since the problem is symmetric about both X and Y axes (the source, geometry, "
                     "and boundary conditions are all symmetric), the solution should exhibit "
                     "reflective symmetry about both axes. This means that φ(i,j) = φ(m+1-i,j) "
                     "and φ(i,j) = φ(i,n+1-j) for all valid i,j.";

    documentTestResults(test1, desc1, hyp1, test1Passed, output1);

    //--------------------------------------------------------------------------
    // Test 2: Rotational Symmetry
    //--------------------------------------------------------------------------
    std::cout << "Test 2: Rotational Symmetry\n";
    std::cout << "-------------------------\n";

    // Create test case with rotationally symmetric source
    int m2 = 8, n2 = 8;  // Must be square grid
    double a2 = 100.0, b2 = 100.0;
    std::string test2 = "rotational_symmetry";
    std::string output2 = test2 + "_output.txt";

    createSymmetricTestCase(test2 + ".txt", m2, n2, a2, b2, D, sigma_a, "rotational");

    // Solve the problem
    DiffusionSolver solver2;
    solver2.readInput(test2 + ".txt");
    auto solution2 = solver2.solve();
    solver2.writeOutputWithDetails(solution2, output2, 0.0, 0, {}, 0.0);

    // Verify rotational symmetry
    bool test2Passed = verifyRotationalSymmetry(solution2, m2, n2);
    std::cout << "Test 2 " << (test2Passed ? "PASSED" : "FAILED") << "\n\n";

    // Document the test
    std::string desc2 = "This test verifies that the solver correctly preserves rotational symmetry "
                      "when the source distribution has rotational symmetry about the center of "
                      "the domain. Four point sources with value 10.0 are placed at positions "
                      "that are symmetric under 90° rotation. The physical parameters are D = " +
                      std::to_string(D) + " and Σₐ = " + std::to_string(sigma_a) + ". The square "
                      "domain is " + std::to_string(a2) + " cm × " + std::to_string(b2) +
                      " cm with a " + std::to_string(m2) + "×" + std::to_string(n2) + " grid.";

    std::string hyp2 = "Since the problem has rotational symmetry (the source, geometry, and boundary "
                     "conditions are all symmetric under 90° rotation about the center), the solution "
                     "should exhibit the same rotational symmetry. This means that the flux values "
                     "at positions that are 90° rotations of each other should be equal.";

    documentTestResults(test2, desc2, hyp2, test2Passed, output2);

    //--------------------------------------------------------------------------
    // Summary
    //--------------------------------------------------------------------------
    std::cout << "Verification Test Summary\n";
    std::cout << "=======================\n";
    std::cout << "Test 1 (Reflective Symmetry): " << (test1Passed ? "PASSED" : "FAILED") << "\n";
    std::cout << "Test 2 (Rotational Symmetry): " << (test2Passed ? "PASSED" : "FAILED") << "\n";

    return 0;
}
