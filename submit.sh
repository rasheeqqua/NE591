#!/bin/tcsh
#BSUB -n 1                     # Request 1 core (single-threaded application)
#BSUB -R span[hosts=1]         # All processes on same node
#BSUB -W 30:00                 # 30-minutes wall clock limit
#BSUB -J sn_solver             # Job name
#BSUB -o sn_solver.%J.out      # Output file name (%J is job ID)
#BSUB -e sn_solver.%J.err      # Error file name (%J is job ID)

# Load the required module
module load gcc

# Run all three programs sequentially
echo "Running S_N Solver..."
./SNSolver
echo "S_N Solver completed."