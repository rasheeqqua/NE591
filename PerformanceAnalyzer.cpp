#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <cstdlib>
#include <chrono>
// #include "DiffusionSolver.cpp"  // You might no longer directly link solver here if you do system calls

class PerformanceAnalyzer {
public:
    // Example function that runs parallel tests for n=128,256,512,1024 and procs=1,4,16,64
    void runParallelTests() {
        // These are the grid sizes for the test
        std::vector<int> gridSizes = {128, 256, 512, 1024};
        // Processor counts
        std::vector<int> procCounts = {1, 4, 16, 64};
            // Example diffusion parameters
        double D = 1.0;
        double SigmaA = 0.1;
        double a = 10.0;
        double b = 10.0;
        int maxIter = 10000;
        double tol = 1.0e-5;

        // We loop over procCounts, then over gridSizes
        for (auto p : procCounts) {
            for (auto n : gridSizes) {
                // 1) Create an input file dynamically:
                std::string inputName = "input_parallel_n" + std::to_string(n) + "_p" + std::to_string(p) + ".txt";
                {
                    std::ofstream finput(inputName);
                    // method=4 => Parallel Jacobi
                    finput << 4 << "\n";
                    // maxIter, tol, an unused omega
                    finput << maxIter << " " << tol << " 1.0\n";
                    // a, b
                    finput << a << " " << b << "\n";
                    // m=n, n=n
                    finput << n << " " << n << "\n";
                    // D, sigma_a
                    finput << D << " " << SigmaA << "\n";
                    // now fill the source terms, for example all 10.0 or something
                    for (int i = 0; i < n; i++) {
                        for (int j = 0; j < n; j++) {
                            // for demonstration, put a constant or zero
                            finput << "10.0 ";
                        }
                        finput << "\n";
                    }
                }

                // 2) Run the solver with mpirun -n p
                std::string cmd = "mpirun -n " + std::to_string(p) + " ./solver "
                                  + " > out_n" + std::to_string(n)
                                  + "_p" + std::to_string(p) + ".log";

                std::cout << "Running n=" << n << " with " << p << " processes...\n";
                auto t0 = std::chrono::high_resolution_clock::now();
                // system call
                int retVal = std::system(cmd.c_str());
                auto t1 = std::chrono::high_resolution_clock::now();
                double elapsed = std::chrono::duration<double>(t1 - t0).count();

                // 3) Save or print the performance data
                std::cout << "Completed in " << elapsed << " seconds. Return code=" << retVal << "\n";
                // You can also parse “out_nX_pY.log” or “output.txt” for iteration count, final error, etc.
            }
        }
    }
};

int main() {
    PerformanceAnalyzer analyzer;
    analyzer.runParallelTests();
    return 0;
}