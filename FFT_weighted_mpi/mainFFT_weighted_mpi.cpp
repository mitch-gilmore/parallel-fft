#include <mpi.h>
#include "fft_weighted.h"

int main(int argc, char** argv) {
    MPI_Init(&argc, &argv);

    int world_size;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);

    int world_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

    int N = Nfft;
    Complex *omega = new Complex[N];
    Complex *F = new Complex[N];
    Complex *Ftilde = new Complex[N];
    Complex *Fnew = new Complex[N];
    Complex *Fsave = new Complex[N];

    if (world_rank == 0) {
    makePhase(omega, N);
        for (int x = 0; x < N; x++) {
            // Use the same signal function as in the non-MPI version
            double real_part = 2 * sin(2.0 * PI * x / (double) N) + 4 * cos(2.0 * PI * 3.0 * x / (double) N) - 8 * cos(2.0 * PI * 4.0 * x / (double) N);
            F[x] = Complex(real_part, 0.0);  // Assume the imaginary part is 0 as the original signal is real
            Fsave[x] = F[x];
        }
    }

    // Broadcast omega to all processes with error checking
    int err = MPI_Bcast(omega, N, MPI_CXX_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD);
    if (err != MPI_SUCCESS) {
        std::cerr << "MPI_Bcast failed" << std::endl;
        MPI_Abort(MPI_COMM_WORLD, err);
    }

    // Scatter the initialized data to all processes
    Complex *local_F = new Complex[N / world_size];
    MPI_Scatter(F, N / world_size, MPI_CXX_DOUBLE_COMPLEX, local_F, N / world_size, MPI_CXX_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD);

    // Perform FFT on each piece
    Complex *local_Ftilde = new Complex[N / world_size];
    FFT(local_Ftilde, local_F, omega, N / world_size);

    // Gather results back to root
    MPI_Gather(local_Ftilde, N / world_size, MPI_CXX_DOUBLE_COMPLEX, Ftilde, N / world_size, MPI_CXX_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD);

    // Synchronize before printing results
    MPI_Barrier(MPI_COMM_WORLD);

    if (world_rank == 0) {
        std::cout << "FFT Results:" << std::endl;
        for (int i = 0; i < N; i++) {
            std::cout << "Ftilde[" << i << "] = " << Ftilde[i] << std::endl;
        }
    }

    delete[] omega;
    delete[] F;
    delete[] Ftilde;
    delete[] Fnew;
    delete[] Fsave;
    delete[] local_F;
    delete[] local_Ftilde;

    MPI_Finalize();
    return 0;
}
