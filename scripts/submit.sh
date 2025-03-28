#!/bin/tcsh
#BSUB -n 1                     # set the number of processes to 1
#BSUB -R span[hosts=1]         # All processes on same node
#BSUB -W 4                     # 5-minutes wall clock limit
#BSUB -q short
#BSUB -J outlab10              # Job name
#BSUB -o outlab10.%J.out       # Output file name (%J is job ID)
#BSUB -e outlab10.%J.err       # Error file name (%J is job ID)

# Load the required module
module load openmpi-gcc

# Run the tests with different process counts
cd ..
mpirun ./linear_solvers
mpirun ./performance_test