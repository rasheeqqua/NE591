#!/bin/tcsh
#BSUB -n 8                     # Request 8 cores maximum
#BSUB -R span[hosts=1]         # All processes on same node
#BSUB -W 30:00                 # 30-minutes wall clock limit
#BSUB -J parallel_sn           # Job name
#BSUB -o parallel_sn.%J.out    # Output file name (%J is job ID)
#BSUB -e parallel_sn.%J.err    # Error file name (%J is job ID)

# Load the required modules
module load gcc
module load openmpi

# Parallel run for the S_N Solver (P=1)
echo "Running the S_N Solver with 2 process..."
mpirun -n 2 ./SNSolver > output_p2.log

# Prepare the input file for the test suite
echo "Creating the input file for testing purposes..."
mpirun -n 1 ./SNSolverTest

# Run the test suite
echo "Running test with 1 processes..."
mpirun -n 1 ./SNSolver > test_p1.log

# P=2 run
echo "Running test with 2 processes..."
mpirun -n 2 ./SNSolver > test_p2.log

# P=4 run
echo "Running test with 4 processes..."
mpirun -n 4 ./SNSolver > test_p4.log

# P=8 run
echo "Running test with 8 processes..."
mpirun -n 8 ./SNSolver > test_p8.log

echo "All runs completed."