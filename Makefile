CXX       = mpicxx
CXXFLAGS  = -O2 -std=c++11
SRCS   = ./iterative/main.cpp ./iterative/DiffusionSolver.cpp ./LUP/LUPFactorization.cpp ./LUP/substitution.cpp ./LUP/ApplyPermutationMatrix.cpp
OBJS   = $(SRCS:.cpp=.o)

all: solver performance

solver: $(SRCS)
$(CXX) $(CXXFLAGS) -o solver $(SRCS)

performance: PerformanceAnalyzer.cpp ./iterative/DiffusionSolver.cpp ./LUP/LUPFactorization.cpp ./LUP/substitution.cpp ./LUP/ApplyPermutationMatrix.cpp
$(CXX) $(CXXFLAGS) -o performance PerformanceAnalyzer.cpp ./iteraive/DiffusionSolver.cpp ./LUP/LUPFactorization.cpp ./LUP/substitution.cpp ./LUP/ApplyPermutationMatrix.cpp

clean:
rm -f *.o solver performance *.log *.txt