#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <fftw3-mpi.h>

#define NX (int)pow(2, 5)
#define NY (int)pow(2, 5)
#define NZ (int)pow(2, 5)

#define SIZE NX * NY * NZ

int main(int argc, char **argv) {
    int rank, size;
    fftw_complex *in, *out;
    fftw_plan plan;

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    // Allocate memory for the input and output arrays
    in = fftw_alloc_complex(SIZE / size);
    out = fftw_alloc_complex(SIZE / size);

    // Create the FFTW plan
    plan = fftw_mpi_plan_dft_3d(NX, NY, NZ, in, out, MPI_COMM_WORLD, FFTW_FORWARD, FFTW_ESTIMATE);

    // Generate some example data for in array
    

    // Execute the FFTW plan
    fftw_execute(plan);

    // Gather results if necessary

    // Clean up
    fftw_destroy_plan(plan);
    fftw_free(in);
    fftw_free(out);

    MPI_Finalize();

    return 0;
}