#include <mpi.h>
#include <iostream>
#include <chrono>
#include "DiffusionSolver.cpp"

int main(int argc, char* argv[]) {
    MPI_Init(&argc, &argv);

    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    // We will create a DiffusionSolver instance on each rank,
    // but only rank=0 actually reads from file. Then we broadcast.
    DiffusionSolver solver;

    bool inputValid = true;
    int methodFlag = -1;  // to broadcast

    // Only rank 0 tries to read input
    if (rank == 0) {
        std::string inputFile = "input.txt"; // or pass via argv
        if (!solver.readInput(inputFile)) {
            inputValid = false;
        }
        methodFlag = solver.getFlag();
    }

    // Broadcast whether input is valid
    MPI_Bcast(&inputValid, 1, MPI_C_BOOL, 0, MPI_COMM_WORLD);
    if (!inputValid) {
        if (rank == 0) {
            std::cerr << "Invalid input encountered. Aborting.\n";
        }
        MPI_Finalize();
        return 1;
    }

    // Broadcast the method flag
    MPI_Bcast(&methodFlag, 1, MPI_INT, 0, MPI_COMM_WORLD);

    // If method !=4, we do a “serial” run on rank=0
    // but if multiple ranks are launched, ranks>0 basically do nothing.
    // If method=4, then all ranks participate in the parallel approach.

    // Next, we need to broadcast the rest of the solver’s parameters
    // from rank 0 to all ranks so that each solver object has the same data.
    // Easiest is to gather them up in local variables on rank 0, then broadcast.
    int maxIter=0, gridM=0, gridN=0;
    double tol=0.0, w=0.0, rectA=0.0, rectB=0.0, diffCoef=0.0, removalXS=0.0;
    if (rank == 0) {
        maxIter    = solver.getMaxIterations();
        tol        = solver.getTolerance();
        w          = solver.getOmega();
        auto dims  = solver.getGridDimensions();
        gridM      = dims.first;
        gridN      = dims.second;
        rectA      = solver.compareSolutions; // fix incorrect usage below
        // Actually, we’ll just parse a/b, D, sigma_a from the solver by adding small getters
        // or we can add them quickly here. Let’s do that.
    }

    // We need to add short getter functions for "a, b, D, sigma_a, q"
    // or we can replicate the logic.  For demonstration:
    // Let’s just add the needed getters in the solver for demonstration.
    // We'll modify the solver code slightly by adding:
    //   double getA() const {return a;}
    //   double getB() const {return b;}
    //   double getD() const {return D;}
    //   double getSigmaA() const {return sigma_a;}
    //   const std::vector<std::vector<double>>& getQ() const {return q;}

    // ... now rank 0 collects them:
    double tmpA=0.0, tmpB=0.0, tmpD=0.0, tmpSigmaA=0.0;
    std::vector<std::vector<double>> tmpQ;

    if (rank == 0) {
        tmpA      = solver.getA();
        tmpB      = solver.getB();
        tmpD      = solver.getD();
        tmpSigmaA = solver.getSigmaA();
        tmpQ      = solver.getQ();  // entire source
        auto dims = solver.getGridDimensions();
        gridM     = dims.first;
        gridN     = dims.second;
        maxIter   = solver.getMaxIterations();
        tol       = solver.getTolerance();
        w         = solver.getOmega();
    }

    // First broadcast all scalars
    MPI_Bcast(&maxIter,   1, MPI_INT,    0, MPI_COMM_WORLD);
    MPI_Bcast(&gridM,     1, MPI_INT,    0, MPI_COMM_WORLD);
    MPI_Bcast(&gridN,     1, MPI_INT,    0, MPI_COMM_WORLD);
    MPI_Bcast(&tol,       1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&w,         1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&tmpA,      1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&tmpB,      1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&tmpD,      1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&tmpSigmaA, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    // Next broadcast the source array. This is size M*N.
    // We'll flatten it for broadcast.
    if (methodFlag == 4) {
        // We do want all ranks to have q.
        tmpQ.resize(gridM, std::vector<double>(gridN, 0.0));
    }
    for (int i = 0; i < gridM; i++) {
        if (rank == 0) {
            MPI_Bcast(tmpQ[i].data(), gridN, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        } else {
            MPI_Bcast(tmpQ[i].data(), gridN, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        }
    }

    // Now each rank can do solver.setParameters() with the broadcasted data
    solver.setParameters(methodFlag, maxIter, tol, w,
                         tmpA, tmpB, gridM, gridN,
                         tmpD, tmpSigmaA, tmpQ);

    // measure time (only rank 0 or we do global)
    auto startTime = std::chrono::high_resolution_clock::now();

    int iterations = 0;
    double finalError = 0.0;
    bool converged = false;

    std::vector<std::vector<double>> solution;

    // If methodFlag != 4 => only rank=0 calls solve
    if (methodFlag != 4) {
        if (rank == 0) {
            solution = solver.solve(iterations, finalError, converged);
        }
    } else {
        // methodFlag=4 => parallel approach
        // All ranks call solver.solve, which in turn calls solveParallelPointJacobiMPI,
        // distributing the data.  Only rank=0 will receive the final global solution.
        solution = solver.solve(iterations, finalError, converged);
    }

    auto endTime = std::chrono::high_resolution_clock::now();
    double executionTime = std::chrono::duration<double>(endTime - startTime).count();

    // Rank 0 writes output
    if (rank == 0) {
        solver.writeOutput(solution, "output.txt", iterations, finalError, converged, executionTime);
        // Optionally compare with LUP for small grids
        if ((methodFlag != 0) && converged && (gridM*gridN <= 400)) {
            // compare code from your older snippet
            DiffusionSolver lupSolver;
            // Re-read input or replicate solver’s data
            lupSolver.setParameters(0,0,0.0,0.0, tmpA,tmpB, gridM,gridN, tmpD, tmpSigmaA, tmpQ);
            int lupIt=0; double lupErr=0.0; bool lupConv=false;
            auto lupSol = lupSolver.solve(lupIt, lupErr, lupConv);
            if (lupConv) {
                double maxDiff = solver.compareSolutions(solution, lupSol);
                std::cout << "Max rel difference vs LUP: " << maxDiff << "\n";
            }
        }
        std::cout << "Done. Results in output.txt\n";
    }

    MPI_Finalize();
    return 0;
}