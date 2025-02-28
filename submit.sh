#!/bin/tcsh
#BSUB -n 1                     # Request 1 core (single-threaded application)
#BSUB -R span[hosts=1]         # All processes on same node
#BSUB -W 30:00                 # 30-minutes wall clock limit
#BSUB -J diffusion_solver       # Job name
#BSUB -o diffusion_solver.%J.out # Output file name (%J is job ID)
#BSUB -e diffusion_solver.%J.err # Error file name (%J is job ID)

# Load the required module
module load gcc

# Run all three programs sequentially
echo "Running DiffusionSolver..."
./DiffusionSolver
echo "DiffusionSolver completed."

echo "Running PerformanceAnalyzer..."
./PerformanceAnalyzer
echo "PerformanceAnalyzer completed."

echo "Running VerificationTest..."
./VerificationTest
echo "VerificationTest completed."

echo "All programs executed successfully."