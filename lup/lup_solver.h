#ifndef LUP_SOLVER_H
#define LUP_SOLVER_H

#include <vector>

// LUP decomposition of matrix A into L, U, and P
void lupDecomposition(const std::vector<std::vector<double>>& A,
                      std::vector<std::vector<double>>& L,
                      std::vector<std::vector<double>>& U,
                      std::vector<int>& P);

// Solve Ax = b using LUP decomposition (L, U, P)
std::vector<double> lupSolveSystem(const std::vector<std::vector<double>>& L,
                                  const std::vector<std::vector<double>>& U,
                                  const std::vector<int>& P,
                                  const std::vector<double>& b);

// Combined function to solve Ax = b using LUP decomposition
void lupSolve(const std::vector<std::vector<double>>& A,
             const std::vector<double>& b,
             std::vector<double>& x);

#endif // LUP_SOLVER_H