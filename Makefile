CXX = g++
CXXFLAGS = -std=c++11 -Wall -O2

SRCS = main.cpp \
       matrix_modules/matrix_operations.cpp \
       matrix_modules/verify_positive_definite_matrix.cpp \
       LUP/LUP_main.cpp \
       LUP/LUPFactorization.cpp \
       LUP/ApplyPermutationMatrix.cpp \
       LUP/substitution.cpp \
       SOR/SOR_main.cpp \
       CG/CG_main.cpp

OBJS = $(SRCS:.cpp=.o)

all: linear_solvers

linear_solvers: $(OBJS)
	$(CXX) $(CXXFLAGS) -o $@ $(OBJS)

%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

clean:
	rm -f $(OBJS) linear_solvers