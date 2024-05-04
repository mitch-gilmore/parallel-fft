#!/bin/bash
# Attempt to load the necessary MPI module
module load openmpi 2>/dev/null

# Check if the module load was successful by checking if the MPI compiler is available
if ! command -v mpicxx &> /dev/null; then
    echo "Failed to load MPI module or MPI compiler not found."
    exit 1
fi

echo "MPI compiler loaded successfully."

# Compile the program
make clean
make

# Check if the executable exists
if [ ! -f "./fft_mpi" ]; then
    echo "Compilation failed or executable not found."
    exit 1
fi

# Run the executable with MPI
mpirun -np 4 ./fft_mpi > fft_results.txt

echo "Test completed. Results saved to fft_results.txt."
