CXX = g++
CXXFLAGS = -std=c++11 -Wall -O2

all: cg_solver

cg_solver: main.o matrix_operations.o
	$(CXX) $(CXXFLAGS) -o cg_solver main.o matrix_operations.o

main.o: main.cpp matrix_operations.h
	$(CXX) $(CXXFLAGS) -c main.cpp

matrix_operations.o: matrix_operations.cpp matrix_operations.h
	$(CXX) $(CXXFLAGS) -c matrix_operations.cpp

clean:
	rm -f *.o cg_solver