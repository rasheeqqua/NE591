//
// Steady State One-Speed Diffusion Equation Solver - Milestone 4
// Author: Hasibul H. Rasheeq
// Date: March 22, 2025
// Version: 2.0
//
// This program solves the steady-state, one-speed diffusion equation
// in a 2D rectangular region with vacuum boundary conditions using
// either direct (LUP) or iterative methods (Point-Jacobi, Gauss-Seidel, SOR) or parallel iterative methods (parallel P-J).
//

#include <iostream>
#include <fstream>
#include <iomanip>
#include <vector>
#include <string>
#include <cmath>
#include <chrono>
#include <algorithm>
#include <ctime>

// Include MPI
#include <mpi.h>

// LUP solver components (unchanged from Milestone 3)
#include "../LUP/LUPFactorization.cpp"
#include "../LUP/substitution.cpp"
#include "../LUP/ApplyPermutationMatrix.cpp"

class DiffusionSolver {
private:
    // Same member data as before
    int flag;
    double a, b;
    int m, n;
    double D;
    double sigma_a;
    std::vector<std::vector<double>> q;
    double delta, gamma;

    int maxIterations;
    double tolerance;
    double omega;

    // Helper function to convert 2D indices to 1D (used by LUP)
    int idx(int i, int j) const {
        return (i - 1) * n + (j - 1);
    }

    // Applies diffusion operator (used locally by older methods)
    void applyDiffusionOperator(const std::vector<std::vector<double>>& phi,
                                std::vector<std::vector<double>>& result) const
    {
        for (int i = 1; i <= m; ++i) {
            for (int j = 1; j <= n; ++j) {
                double center = phi[i][j];
                double left   = (i > 1) ? phi[i-1][j] : 0.0;
                double right  = (i < m) ? phi[i+1][j] : 0.0;
                double bottom = (j > 1) ? phi[i][j-1] : 0.0;
                double top    = (j < n) ? phi[i][j+1] : 0.0;

                result[i][j] = -D * ((left - 2*center + right) / (delta*delta) +
                                     (bottom - 2*center + top) / (gamma*gamma))
                               + sigma_a * center;
            }
        }
    }

    // Calculate relative or absolute error (whichever is larger)
    double calculateError(const std::vector<std::vector<double>>& current,
                          const std::vector<std::vector<double>>& previous) const
    {
        double maxError = 0.0;
        double eps = 1e-15; // protect from division by zero

        for (int i = 1; i <= m; ++i) {
            for (int j = 1; j <= n; ++j) {
                double denom = std::abs(current[i][j]);
                if (denom > eps) {
                    double relError = std::abs(current[i][j] - previous[i][j]) / denom;
                    maxError = std::max(maxError, relError);
                } else {
                    // fallback to absolute difference if near zero
                    double absError = std::abs(current[i][j] - previous[i][j]);
                    maxError = std::max(maxError, absError);
                }
            }
        }
        return maxError;
    }

    // Returns max residual ||A*phi - q||∞ for a given phi
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

    //------------------------------------------------------------------
    //  Serial Iterative Methods (Jacobi, Gauss-Seidel, SOR)
    //------------------------------------------------------------------
    bool solvePointJacobi(std::vector<std::vector<double>>& phi,
                          int& iterations, double& finalError)
    {
        std::vector<std::vector<double>> phi_old(m+2, std::vector<double>(n+2, 0.0));
        double diagCoef = 2.0*D/(delta*delta) + 2.0*D/(gamma*gamma) + sigma_a;

        for (iterations = 0; iterations < maxIterations; ++iterations) {
            phi_old = phi;
            // Jacobi update
            for (int i = 1; i <= m; ++i) {
                for (int j = 1; j <= n; ++j) {
                    double left   = (i > 1) ? phi_old[i-1][j] : 0.0;
                    double right  = (i < m) ? phi_old[i+1][j] : 0.0;
                    double bottom = (j > 1) ? phi_old[i][j-1] : 0.0;
                    double top    = (j < n) ? phi_old[i][j+1] : 0.0;

                    double numerator = q[i-1][j-1]
                                       + D*(left + right)/(delta*delta)
                                       + D*(bottom + top)/(gamma*gamma);
                    phi[i][j] = numerator / diagCoef;
                }
            }

            finalError = calculateError(phi, phi_old);
            if (finalError < tolerance) {
                return true;
            }
        }
        return false;
    }

    bool solveGaussSeidel(std::vector<std::vector<double>>& phi,
                          int& iterations, double& finalError)
    {
        std::vector<std::vector<double>> phi_old(m+2, std::vector<double>(n+2, 0.0));
        double diagCoef = 2.0*D/(delta*delta) + 2.0*D/(gamma*gamma) + sigma_a;

        for (iterations = 0; iterations < maxIterations; ++iterations) {
            phi_old = phi;
            for (int i = 1; i <= m; ++i) {
                for (int j = 1; j <= n; ++j) {
                    double left   = (i > 1) ? phi[i-1][j] : 0.0;
                    double right  = (i < m) ? phi_old[i+1][j] : 0.0;
                    double bottom = (j > 1) ? phi[i][j-1] : 0.0;
                    double top    = (j < n) ? phi_old[i][j+1] : 0.0;

                    double numerator = q[i-1][j-1]
                                       + D*(left + right)/(delta*delta)
                                       + D*(bottom + top)/(gamma*gamma);

                    phi[i][j] = numerator / diagCoef;
                }
            }

            finalError = calculateError(phi, phi_old);
            if (finalError < tolerance) {
                return true;
            }
        }
        return false;
    }

    bool solveSOR(std::vector<std::vector<double>>& phi,
                  int& iterations, double& finalError)
    {
        std::vector<std::vector<double>> phi_old(m+2, std::vector<double>(n+2, 0.0));
        double diagCoef = 2.0*D/(delta*delta) + 2.0*D/(gamma*gamma) + sigma_a;

        for (iterations = 0; iterations < maxIterations; ++iterations) {
            phi_old = phi;
            for (int i = 1; i <= m; ++i) {
                for (int j = 1; j <= n; ++j) {
                    double left   = (i > 1) ? phi[i-1][j] : 0.0;
                    double right  = (i < m) ? phi_old[i+1][j] : 0.0;
                    double bottom = (j > 1) ? phi[i][j-1] : 0.0;
                    double top    = (j < n) ? phi_old[i][j+1] : 0.0;

                    double phi_gs = (q[i-1][j-1]
                                    + D*(left + right)/(delta*delta)
                                    + D*(bottom + top)/(gamma*gamma)) / diagCoef;

                    // SOR update
                    phi[i][j] = (1.0 - omega) * phi_old[i][j] + omega * phi_gs;
                }
            }

            finalError = calculateError(phi, phi_old);
            if (finalError < tolerance) {
                return true;
            }
        }
        return false;
    }

    //------------------------------------------------------------------
    //  Parallel Point-Jacobi using MPI (flag = 4)
    //------------------------------------------------------------------
    //  This method implements Items #1..#6 from your milestone instructions.
    //  Manager -> read & broadcast data; subdomain decomposition; halo exchanges;
    //  gather final solution on rank 0, etc.
    bool solveParallelPointJacobiMPI(std::vector<std::vector<double>>& globalPhi,
                                     int& iterations,
                                     double& finalError)
    {
        // MPI essentials
        int rank, size;
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        MPI_Comm_size(MPI_COMM_WORLD, &size);

        // For simplicity, we restrict to square mesh (m == n).
        // We also require sqrt(size) and n/sqrt(size) to be integers.
        if (m != n) {
            if (rank == 0) {
                std::cerr << "Error: For the parallel code (flag=4), we require m=n.\n";
                std::cerr << "Aborting.\n";
            }
            return false;
        }
        int nGlobal = n;  // since m = n

        // Check divisibility
        int rootP = (int)std::sqrt((double)size);
        if (rootP * rootP != size) {
            if (rank == 0) {
                std::cerr << "Error: sqrt(P) must be an integer.\n";
            }
            return false;
        }
        if (nGlobal % rootP != 0) {
            if (rank == 0) {
                std::cerr << "Error: n must be divisible by sqrt(P).\n";
            }
            return false;
        }

        // Determine local subdomain extents
        // We'll do a 2D block decomposition: rankX = rank % rootP, rankY = rank / rootP
        int rankX = rank % rootP;
        int rankY = rank / rootP;

        int blockSize = nGlobal / rootP; // subdomain dimension

        // local domain indices in [1..n], but we partition among ranks
        int iStart = rankX * blockSize + 1;
        int iEnd   = (rankX + 1) * blockSize;
        int jStart = rankY * blockSize + 1;
        int jEnd   = (rankY + 1) * blockSize;

        // Local size for convenience
        int localNX = iEnd - iStart + 1;
        int localNY = jEnd - jStart + 1;

        // Allocate local subdomain, plus a one-cell halo in each direction
        // We'll store these with indices from 0..localNX+1 and 0..localNY+1
        // but remember to map to global indices carefully.
        std::vector<std::vector<double>> phi_old(localNX+2, std::vector<double>(localNY+2, 0.0));
        std::vector<std::vector<double>> phi_new(localNX+2, std::vector<double>(localNY+2, 0.0));
        std::vector<std::vector<double>> q_local(localNX, std::vector<double>(localNY, 0.0));

        // Copy the local portion of q from global to q_local
        // Only rank 0 has the legitimate global "q" array,
        // so we broadcast slices of q to each rank. We can use MPI_Scatterv or
        // we can do a manual send.  For simplicity, we’ll do a manual approach:
        // Rank 0 sends each block of q to the appropriate rank.
        // However, since we already broadcasted input data in main.cpp,
        // each rank has "q" in the solver. We just have to copy the sub-block.

        // We assume 'q' is accessible on all ranks (because main broadcasted it),
        // so each rank picks out the sub-block it needs.
        for (int lx = 0; lx < localNX; lx++) {
            int globalX = iStart + lx - 1; // -1 because q is 0-based
            for (int ly = 0; ly < localNY; ly++) {
                int globalY = jStart + ly - 1;
                q_local[lx][ly] = q[globalX][globalY];
            }
        }

        // Diagonal coefficient
        double diagCoef = 2.0 * D/(delta*delta) + 2.0 * D/(gamma*gamma) + sigma_a;

        // Start iterative loop
        iterations = 0;
        finalError = 0.0;

        // We will do up to maxIterations iterations
        // On each iteration:
        //   - update local cells using old data
        //   - compute local max error
        //   - do MPI_Allreduce to get global max
        //   - check tolerance
        //   - do halo exchange before next iteration
        for (iterations = 0; iterations < maxIterations; ++iterations) {
            // Copy phi_new into phi_old (or vice versa) at iteration start
            // so we can do the Jacobi update
            phi_old = phi_new;

            // Perform local Jacobi updates
            double localMaxChange = 0.0;
            for (int lx = 1; lx <= localNX; lx++) {
                for (int ly = 1; ly <= localNY; ly++) {
                    // Global indices
                    int iGlobal = iStart + (lx - 1);
                    int jGlobal = jStart + (ly - 1);

                    // Grab neighbors from phi_old
                    double left   = (lx == 1)          ? ((iGlobal == 1)         ? 0.0 : phi_old[lx-1][ly]) : phi_old[lx-1][ly];
                    double right  = (lx == localNX)     ? ((iGlobal == nGlobal)   ? 0.0 : phi_old[lx+1][ly]) : phi_old[lx+1][ly];
                    double bottom = (ly == 1)          ? ((jGlobal == 1)         ? 0.0 : phi_old[lx][ly-1]) : phi_old[lx][ly-1];
                    double top    = (ly == localNY)     ? ((jGlobal == nGlobal)   ? 0.0 : phi_old[lx][ly+1]) : phi_old[lx][ly+1];

                    double numerator = q_local[lx-1][ly-1]
                                       + D*(left + right)/(delta*delta)
                                       + D*(bottom + top)/(gamma*gamma);

                    double oldVal = phi_old[lx][ly];
                    double newVal = numerator / diagCoef;
                    phi_new[lx][ly] = newVal;

                    // Compute local change (relative or absolute)
                    double denom = std::abs(newVal);
                    double diff = std::abs(newVal - oldVal);
                    double localChange = (denom > 1.0e-15) ? (diff / denom) : diff;
                    localMaxChange = std::max(localMaxChange, localChange);
                }
            }

            // Reduce to find global maximum change
            double globalMaxChange = 0.0;
            MPI_Allreduce(&localMaxChange, &globalMaxChange, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);

            // Check for convergence
            if (globalMaxChange < tolerance) {
                finalError = globalMaxChange;
                break; // converged
            }

            // Otherwise, do halo exchange before next iteration
            // We'll exchange boundaries in left-right direction, then top-bottom.

            // 1) Horizontal neighbors
            //    Send left column to left neighbor, receive from right, etc.
            // left neighbor rank
            int leftNeighbor  = (rankX > 0)          ? (rankY*rootP + (rankX - 1)) : MPI_PROC_NULL;
            int rightNeighbor = (rankX < rootP - 1)  ? (rankY*rootP + (rankX + 1)) : MPI_PROC_NULL;

            // buffers for send/recv
            std::vector<double> sendLeft(localNY), recvRight(localNY);
            std::vector<double> sendRight(localNY), recvLeft(localNY);

            // Pack sendLeft (column lx=1)
            for (int ly = 1; ly <= localNY; ly++) {
                sendLeft[ly-1]  = phi_new[1][ly];
                sendRight[ly-1] = phi_new[localNX][ly];
            }
            // Send to left, receive from right
            MPI_Sendrecv(sendLeft.data(),  localNY, MPI_DOUBLE, leftNeighbor,  101,
                         recvRight.data(), localNY, MPI_DOUBLE, rightNeighbor, 102,
                         MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            // Send to right, receive from left
            MPI_Sendrecv(sendRight.data(), localNY, MPI_DOUBLE, rightNeighbor, 103,
                         recvLeft.data(),  localNY, MPI_DOUBLE, leftNeighbor,  104,
                         MPI_COMM_WORLD, MPI_STATUS_IGNORE);

            // Unpack incoming buffers
            // receiving from left => place in column 0
            for (int ly = 1; ly <= localNY; ly++) {
                if (leftNeighbor != MPI_PROC_NULL) {
                    phi_new[0][ly] = recvLeft[ly-1];
                }
                if (rightNeighbor != MPI_PROC_NULL) {
                    phi_new[localNX+1][ly] = recvRight[ly-1];
                }
            }

            // 2) Vertical neighbors
            //    Send bottom row to bottom neighbor, receive from top, etc.
            int bottomNeighbor = (rankY > 0)          ? ((rankY - 1)*rootP + rankX) : MPI_PROC_NULL;
            int topNeighbor    = (rankY < rootP - 1)  ? ((rankY + 1)*rootP + rankX) : MPI_PROC_NULL;

            // buffers
            std::vector<double> sendBottom(localNX), sendTop(localNX);
            std::vector<double> recvTop(localNX), recvBottom(localNX);

            for (int lx = 1; lx <= localNX; lx++) {
                sendBottom[lx-1] = phi_new[lx][1];
                sendTop[lx-1]    = phi_new[lx][localNY];
            }
            // Send, receive
            MPI_Sendrecv(sendBottom.data(), localNX, MPI_DOUBLE, bottomNeighbor, 105,
                         recvTop.data(),    localNX, MPI_DOUBLE, topNeighbor,    106,
                         MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            MPI_Sendrecv(sendTop.data(),    localNX, MPI_DOUBLE, topNeighbor,    107,
                         recvBottom.data(), localNX, MPI_DOUBLE, bottomNeighbor, 108,
                         MPI_COMM_WORLD, MPI_STATUS_IGNORE);

            // Unpack
            for (int lx = 1; lx <= localNX; lx++) {
                if (bottomNeighbor != MPI_PROC_NULL) {
                    phi_new[lx][0] = recvBottom[lx-1];
                }
                if (topNeighbor != MPI_PROC_NULL) {
                    phi_new[lx][localNY+1] = recvTop[lx-1];
                }
            }
        } // end iteration loop

        // If we exited because iterations >= maxIterations, we are not converged
        bool converged = (iterations < maxIterations);

        // The last globalMaxChange (if we never converged) wasn’t stored, so finalError
        // might still be from a prior iteration or zero. We can store it here.
        if (!converged) {
            finalError = 0.0; // We can store best guess or 0
        }

        // Gather local phi_new subdomains into the globalPhi on rank 0.
        // globalPhi has dimension (m+2)x(n+2). We only fill [1..m][1..n].
        // We can gather with MPI_Gatherv in a 2D manner, or do it block by block.
        // Below is a simple manual approach: rank 0 receives each rank’s sub-block.
        // We’ll send the interior portion [1..localNX][1..localNY].
        // We do a loop over ranks from 0..size-1, each rank sends to manager.

        // Send local subdomain data to manager
        std::vector<double> sendbuf(localNX * localNY, 0.0);
        for (int lx = 1; lx <= localNX; lx++) {
            for (int ly = 1; ly <= localNY; ly++) {
                sendbuf[(lx-1)*localNY + (ly-1)] = phi_new[lx][ly];
            }
        }

        if (rank == 0) {
            // rank 0 will receive from all
            for (int src = 0; src < size; src++) {
                // figure out that src’s subdomain
                int rx = src % rootP;
                int ry = src / rootP;
                int iSt = rx*blockSize + 1;
                int iEd = (rx+1)*blockSize;
                int jSt = ry*blockSize + 1;
                int jEd = (ry+1)*blockSize;

                int locNX = iEd - iSt + 1;
                int locNY = jEd - jSt + 1;

                if (src == 0) {
                    // already have local data in sendbuf
                    for (int lx = 0; lx < locNX; lx++) {
                        for (int ly = 0; ly < locNY; ly++) {
                            globalPhi[iSt + lx][jSt + ly] = sendbuf[lx*locNY + ly];
                        }
                    }
                } else {
                    // receive from 'src'
                    std::vector<double> recvbuf(locNX * locNY, 0.0);
                    MPI_Recv(recvbuf.data(), locNX*locNY, MPI_DOUBLE, src, 200 + src,
                             MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                    for (int lx = 0; lx < locNX; lx++) {
                        for (int ly = 0; ly < locNY; ly++) {
                            globalPhi[iSt + lx][jSt + ly] = recvbuf[lx*locNY + ly];
                        }
                    }
                }
            }
        } else {
            // send local subdomain to rank 0
            MPI_Send(sendbuf.data(), localNX*localNY, MPI_DOUBLE, 0, 200 + rank,
                     MPI_COMM_WORLD);
        }

        return converged;
    }

    //------------------------------------------------------------------
    //  Direct solver (LUP) from Milestone 3 (used for comparison)
    //------------------------------------------------------------------
    std::vector<std::vector<double>> solveLUP() const {
        int size = m * n;
        std::vector<std::vector<double>> A(size, std::vector<double>(size, 0.0));
        std::vector<double> b(size, 0.0);

        // Fill matrix A and vector b
        for (int i = 1; i <= m; ++i) {
            for (int j = 1; j <= n; ++j) {
                int row = idx(i, j);

                double diag = 2.0*D/(delta*delta) + 2.0*D/(gamma*gamma) + sigma_a;
                A[row][row] = diag;

                if (i > 1) A[row][idx(i-1,j)] = -D/(delta*delta);
                if (i < m) A[row][idx(i+1,j)] = -D/(delta*delta);
                if (j > 1) A[row][idx(i,j-1)] = -D/(gamma*gamma);
                if (j < n) A[row][idx(i,j+1)] = -D/(gamma*gamma);

                b[row] = q[i-1][j-1];
            }
        }

        std::vector<std::vector<double>> L(size, std::vector<double>(size, 0.0));
        std::vector<std::vector<double>> U(size, std::vector<double>(size, 0.0));
        std::vector<std::vector<double>> P(size, std::vector<double>(size, 0.0));

        if (!lupFactorize(A, L, U, P)) {
            throw std::runtime_error("LUP factorization failed");
        }

        std::vector<double> Pb(size);
        applyPermutationMatrix(P, b, Pb);

        std::vector<double> y(size), x(size);
        forwardSubstitution(L, Pb, y);
        backSubstitution(U, y, x);

        // Convert back to 2D
        std::vector<std::vector<double>> phi(m+2, std::vector<double>(n+2, 0.0));
        for (int i = 1; i <= m; ++i) {
            for (int j = 1; j <= n; ++j) {
                phi[i][j] = x[idx(i, j)];
            }
        }

        return phi;
    }

public:
    // Read input (only rank=0 typically calls this in parallel scenario)
    bool readInput(const std::string& filename) {
        std::ifstream input(filename);
        if (!input.is_open()) {
          std::cerr << "Error: Cannot open input file " << filename << std::endl;
          return false;
        }

        input >> flag;
        switch (flag) {
            case 0: // LUP direct
                input >> maxIterations >> tolerance >> omega;
                // override for LUP
                maxIterations = 0;
                tolerance = 0.0;
                omega = 0.0;
                break;
            case 1: // Jacobi
            case 2: // Gauss-Seidel
                input >> maxIterations >> tolerance >> omega; // omega not used in 1 or 2
                omega = 0.0;
                break;
            case 3: // SOR
                input >> maxIterations >> tolerance >> omega;
                if (omega <= 0.0 || omega >= 2.0) {
                    std::cerr << "Error: SOR relaxation parameter must be in range (0,2)." << std::endl;
                    return false;
                }
                break;
            case 4: // Parallel Point-Jacobi
                // We will read maxIterations, tolerance, set others as needed
                input >> maxIterations >> tolerance >> omega;
                // For parallel Jacobi, omega is not used
                omega = 0.0;
                break;
            default:
                std::cerr << "Error: Invalid solution method flag. Must be 0-4." << std::endl;
                return false;
        }

        input >> a >> b;
        if (a <= 0 || b <= 0) {
            std::cerr << "Error: Rectangle dimensions invalid (must be > 0)." << std::endl;
            return false;
        }

        input >> m >> n;
        if (m <= 0 || n <= 0) {
            std::cerr << "Error: Grid dimensions must be positive." << std::endl;
            return false;
        }

        input >> D >> sigma_a;
        if (D <= 0) {
            std::cerr << "Error: D must be > 0." << std::endl;
            return false;
        }
        if (sigma_a < 0) {
            std::cerr << "Error: sigma_a must be >= 0." << std::endl;
            return false;
        }

        // read source
        q.resize(m, std::vector<double>(n, 0.0));
        for (int i = 0; i < m; i++) {
            for (int j = 0; j < n; j++) {
                input >> q[i][j];
                if (q[i][j] < 0) {
                    std::cerr << "Error: Negative source term at (" << i+1 << "," << j+1 << ")\n";
                    return false;
                }
            }
        }

        delta = a / (m + 1);
        gamma = b / (n + 1);

        return true;
    }

    // Alternatively set parameters by code
    void setParameters(int methodFlag, int maxIter, double tol, double w,
                       double rectA, double rectB, int gridM, int gridN,
                       double diffCoef, double removalXS,
                       const std::vector<std::vector<double>>& sourceTerms)
    {
        flag = methodFlag;
        maxIterations = maxIter;
        tolerance = tol;
        omega = w;

        a = rectA;
        b = rectB;
        m = gridM;
        n = gridN;

        D = diffCoef;
        sigma_a = removalXS;

        q = sourceTerms;
        delta = a / (m + 1);
        gamma = b / (n + 1);
    }

    // Main “solve” dispatcher (returns final flux, iteration count, error, etc.)
    std::vector<std::vector<double>> solve(int& iterations, double& finalError, bool& converged)
    {
        // Global flux array for all methods including parallel
        // We store boundary extra layers as zeros: dimension (m+2)x(n+2)
        std::vector<std::vector<double>> phi(m+2, std::vector<double>(n+2, 0.0));

        switch (flag) {
            case 0: { // LUP
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
            }
            case 1: { // Point Jacobi
                converged = solvePointJacobi(phi, iterations, finalError);
                break;
            }
            case 2: { // Gauss-Seidel
                converged = solveGaussSeidel(phi, iterations, finalError);
                break;
            }
            case 3: { // SOR
                converged = solveSOR(phi, iterations, finalError);
                break;
            }
            case 4: { // Parallel Point-Jacobi with MPI
                converged = solveParallelPointJacobiMPI(phi, iterations, finalError);
                break;
            }
            default:
                std::cerr << "Error: Invalid flag in solve()." << std::endl;
                converged = false;
                break;
        }
        return phi;
    }

    // Write output (rank=0 typically calls this; safe to do in serial if not method=4)
    void writeOutput(const std::vector<std::vector<double>>& phi,
                     const std::string& filename,
                     int iterations,
                     double finalError,
                     bool converged,
                     double executionTime)
    {
        std::ofstream output(filename);
        if (!output.is_open()) {
            std::cerr << "Error: unable to open output file " << filename << std::endl;
            return;
        }

        // Current date/time
        std::time_t now = std::time(nullptr);
        char timeStr[100];
        std::strftime(timeStr, sizeof(timeStr), "%b %d %Y %H:%M:%S", std::localtime(&now));

        output << "Steady State One-Speed Diffusion Equation Solver - Extended\n";
        output << "Version with MPI Parallel Option (flag=4)\n";
        output << "Author: Hasibul H. Rasheeq\n";
        output << "Date: " << timeStr << "\n";
        output << "----------------------------------------\n\n";

        output << "Solution Method: ";
        switch (flag) {
            case 0: output << "LUP Decomposition (Direct)"; break;
            case 1: output << "Point Jacobi (Serial)"; break;
            case 2: output << "Gauss-Seidel (Serial)"; break;
            case 3: output << "SOR (Serial), omega=" << omega; break;
            case 4: output << "Parallel Point-Jacobi (MPI)"; break;
            default: output << "Unknown"; break;
        }
        output << "\n\n";

        output << "Problem Parameters:\n";
        output << "  Grid dimensions: " << m << " x " << n << " (" << m*n << " unknowns)\n";
        output << "  Rectangle dimensions: " << a << " x " << b << " (cm)\n";
        output << "  Grid spacing: delta=" << delta << ", gamma=" << gamma << "\n";
        output << "  Physical params: D=" << D << " (cm), sigma_a=" << sigma_a << " (1/cm)\n\n";

        if (flag != 0) {
            output << "Iteration Info:\n";
            output << "  Max iterations allowed: " << maxIterations << "\n";
            output << "  Convergence tolerance: " << std::scientific << std::setprecision(6) << tolerance << "\n";
            if (converged) {
                output << "  Converged after " << iterations << " iterations\n";
                output << "  Final error: " << finalError << "\n";
            } else {
                output << "  NOT converged after " << maxIterations << " iterations\n";
                output << "  Last reported error: " << finalError << "\n";
            }
            output << "\n";
        }

        // Performance
        output << "Performance:\n";
        output << "  Execution time: " << std::fixed << std::setprecision(8) << executionTime << " s\n";
        double maxResidual = calculateMaxResidual(phi);
        output << "  Max absolute residual: " << std::scientific << std::setprecision(6) << maxResidual << "\n\n";

        // Print final flux
        output << "Scalar Flux (i,j,phi):\n";
        output << std::scientific << std::setprecision(6);
        for (int i = 0; i <= m+1; i++) {
            for (int j = 0; j <= n+1; j++) {
                output << "i=" << i << " j=" << j << " " << phi[i][j] << "\n";
            }
        }
        output.close();
    }

    // Compare solutions
    double compareSolutions(const std::vector<std::vector<double>>& phi1,
                            const std::vector<std::vector<double>>& phi2)
    {
        double maxDiff = 0.0;
        double eps = 1e-15;
        for (int i = 1; i <= m; ++i) {
            for (int j = 1; j <= n; ++j) {
                double denom = std::max(std::abs(phi1[i][j]), std::abs(phi2[i][j]));
                if (denom > eps) {
                    double relDiff = std::abs(phi1[i][j] - phi2[i][j]) / denom;
                    maxDiff = std::max(maxDiff, relDiff);
                } else {
                    double absDiff = std::abs(phi1[i][j] - phi2[i][j]);
                    maxDiff = std::max(maxDiff, absDiff);
                }
            }
        }
        return maxDiff;
    }

    // Simple getters
    int getFlag() const { return flag; }
    int getMaxIterations() const { return maxIterations; }
    double getTolerance() const { return tolerance; }
    double getOmega() const { return omega; }
    std::pair<int,int> getGridDimensions() const { return {m,n}; }

    double getA() const { return a; }
    double getB() const { return b; }
    double getD() const { return D; }
    double getSigmaA() const { return sigma_a; }
    const std::vector<std::vector<double>>& getQ() const { return q; }
};