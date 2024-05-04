#include "fft.h"
#include <mpi.h>  // Include MPI library

namespace fft
{
    namespace cpu
    {
        fftPlan::fftPlan(int N)
        {
            this->N = N;
            omega = new Complex[N];
            this->makePhase(omega);
        }

        fftPlan::~fftPlan()
        {
            delete[] omega;
        }

        double fftPlan::getWeights(unsigned int k)  // Added function
        {
            // Placeholder logic, adjust accordingly
            return 1.0; // Or any other weighting scheme
        }

        void fftPlan::makePhase(Complex *omega)
        {
            int rank, size;
            MPI_Comm_rank(MPI_COMM_WORLD, &rank);
            MPI_Comm_size(MPI_COMM_WORLD, &size);

            int local_n = N / size;
            int start = rank * local_n;
            int end = start + local_n;

            for (unsigned int k = start; k < end; k++) {
                omega[k] = exp(2.0 * M_PI * I * (double)k * getWeights(k) / (double)this->N);
            }

            // Ensure all phases are constructed
            MPI_Allgather(MPI_IN_PLACE, local_n, MPI_DOUBLE_COMPLEX, omega, local_n, MPI_DOUBLE_COMPLEX, MPI_COMM_WORLD);
        }
    }
}
