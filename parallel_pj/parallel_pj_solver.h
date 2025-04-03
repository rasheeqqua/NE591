#ifndef PARALLEL_PJ_SOLVER_H
#define PARALLEL_PJ_SOLVER_H

#include <vector>
#include <mpi.h>

// Forward declaration of ProblemParameters struct
struct ProblemParameters;

/**
 * Solves the diffusion equation using parallel Point Jacobi method with MPI
 *
 * @param params Problem parameters including grid size, physical parameters, etc.
 * @param local_phi Output - Local solution grid including halo cells
 * @param iterations Output - Number of iterations performed
 * @param final_error Output - Final maximum relative error
 */
void parallelPointJacobiSolve(const ProblemParameters& params,
                             std::vector<std::vector<double>>& local_phi,
                             int& iterations,
                             double& final_error);

#endif // PARALLEL_PJ_SOLVER_H