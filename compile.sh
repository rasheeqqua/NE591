#!/bin/sh

# Clear any existing modules
module purge

# Load required compiler module
module load gcc

# Compile the main program
echo "Compiling S_N Solver..."
g++ -O2 -Wall main.cpp -o SNSolver
if [ $? -eq 0 ]; then
    echo "S_N Solver compilation successful."
else
    echo "S_N Solver compilation failed."
    exit 1
fi

# Compile the test suite
echo "Compiling S_N Solver Test Suite..."
g++ -O2 -Wall test.cpp -o SNSolverTest
if [ $? -eq 0 ]; then
    echo "S_N Solver Test Suite compilation successful."
else
    echo "S_N Solver Test Suite compilation failed."
    exit 1
fi