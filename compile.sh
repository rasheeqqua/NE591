#!/bin/bash

# Clear any existing modules
module purge

# Load required compiler module
module load gcc

# Compile the main program
echo "Compiling DiffusionSolver..."
g++ -O2 -Wall main.cpp -o DiffusionSolver
if [ $? -eq 0 ]; then
    echo "DiffusionSolver compilation successful."
else
    echo "DiffusionSolver compilation failed."
    exit 1
fi

# Compile PerformanceAnalyzer
echo "Compiling PerformanceAnalyzer..."
g++ -O2 -Wall PerformanceAnalyzer.cpp -o PerformanceAnalyzer
if [ $? -eq 0 ]; then
    echo "PerformanceAnalyzer compilation successful."
else
    echo "PerformanceAnalyzer compilation failed."
    exit 1
fi

# Compile VerificationTest
echo "Compiling VerificationTest..."
g++ -O2 -Wall VerificationTest.cpp -o VerificationTest
if [ $? -eq 0 ]; then
    echo "VerificationTest compilation successful."
else
    echo "VerificationTest compilation failed."
    exit 1
fi

echo "All compilations completed successfully."