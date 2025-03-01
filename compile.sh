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