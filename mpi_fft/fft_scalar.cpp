#include "fft.h"
#include <mpi.h>
#include <iostream>

namespace fft {
    namespace cpu {
        namespace scalar {
            using fft::cpu::Complex;
            using fft::cpu::I;

            void BitRevArray(Complex *Frev, Complex *F, int N) {
                int rank, size;
                MPI_Comm_rank(MPI_COMM_WORLD, &rank);
                MPI_Comm_size(MPI_COMM_WORLD, &size);

                int local_n = N / size;
                Complex* local_Frev = new Complex[local_n];
                for (int i = 0; i < local_n; i++) {
                    int idx = reverseBits(rank * local_n + i, N);
                    local_Frev[i] = F[idx];
                }
                MPI_Allgather(local_Frev, local_n, MPI_DOUBLE_COMPLEX, Frev, local_n, MPI_DOUBLE_COMPLEX, MPI_COMM_WORLD);
                delete[] local_Frev;
            }

            unsigned int reverseBits(unsigned int x, unsigned int N) {
                unsigned int xrev = 0;
                unsigned int p = log2(N);
                unsigned int n;
                unsigned int power = N;
                for (unsigned int i = 0; i < p; i++) {
                    n = x % 2;
                    power /= 2;
                    xrev += n * power;
                    x /= 2;
                }
                return xrev;
            }

            void FFT(Complex *Fout, Complex *Fin, fft::cpu::fftPlan plan, bool forward) {
                int rank, size;
                MPI_Comm_rank(MPI_COMM_WORLD, &rank);
                MPI_Comm_size(MPI_COMM_WORLD, &size);

                int log2n = log2(plan.N);
                int local_n = plan.N / size;

                #pragma omp parallel for
                for (unsigned int x = 0; x < local_n; x++) {
                    Fout[reverseBits(rank * local_n + x, plan.N)] = Fin[rank * local_n + x];
                }

                MPI_Allgather(MPI_IN_PLACE, local_n, MPI_DOUBLE_COMPLEX, Fout, local_n, MPI_DOUBLE_COMPLEX, MPI_COMM_WORLD);

                for (int s = 0; s < log2n; s++) {
                    int Bsize = pow(2, s + 1);
                    int Nblocks = plan.N / Bsize;

                    for (int nb = 0; nb < Nblocks; nb++) {
                        for (int pairs = 0; pairs < Bsize / 2; pairs++) {
                            int k = pairs + nb * Bsize;
                            if (k / Bsize % size != rank) continue; // Only compute on the responsible rank

                            Complex phase = (forward) ? plan.omega[Nblocks * pairs] : conj(plan.omega[Nblocks * pairs]);
                            Complex low = Fout[k];
                            Complex high = phase * Fout[k + Bsize / 2];
                            Fout[k] = low + high;
                            Fout[k + Bsize / 2] = low - high;
                        }
                    }
                }
            }
        }
    }
}

// Main function outside any namespace
int main(int argc, char **argv) {
    MPI_Init(&argc, &argv);

    int N = 1024; // FFT size
    fft::cpu::fftPlan plan(N); // Ensure this is accessible

fft::cpu::Complex *input = new fft::cpu::Complex[N];
fft::cpu::Complex *output = new fft::cpu::Complex[N];

// Initialize input
for (int i = 0; i < N; i++) {
	input[i] = fft::cpu::Complex(cos(2 * M_PI * i / N), sin(2 * M_PI * i / N));
}

fft::cpu::scalar::FFT(output, input, plan, true); // Forward FFT

delete[] input;
delete[] output;
delete[] plan.omega;

MPI_Finalize();
    return 0;
}
