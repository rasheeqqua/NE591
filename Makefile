CXX = g++
CXXFLAGS = -std=c++11 -Wall -O2

COMMON_SRCS = matrix_modules/matrix_operations.cpp \
              matrix_modules/verify_positive_definite_matrix.cpp \
              LUP/LUP_main.cpp \
              LUP/LUPFactorization.cpp \
              LUP/ApplyPermutationMatrix.cpp \
              LUP/substitution.cpp \
              SOR/SOR_main.cpp \
              CG/CG_main.cpp \
              PCG/PCG_main.cpp \
              PCG/jacobi_preconditioner.cpp \
              PI/power_iterations.cpp \
              II/inverse_iterations.cpp

MAIN_SRCS = main.cpp $(COMMON_SRCS)
MAIN_OBJS = $(MAIN_SRCS:.cpp=.o)

all: linear_solvers

linear_solvers: $(MAIN_OBJS)
	$(CXX) $(CXXFLAGS) -o $@ $(MAIN_OBJS)

%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

clean:
	rm -f $(MAIN_OBJS) $(PERF_OBJS) linear_solvers