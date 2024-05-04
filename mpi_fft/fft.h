#pragma once

#include <cmath>
#include <complex>
#include <mpi.h> // Include MPI for any potential MPI data types or operations

namespace fft
{
    namespace cpu
    {
        typedef std::complex<double> Complex;

        const Complex I = Complex(0.0, 1.0);

        class fftPlan
        {
        public:
            Complex *omega;
            unsigned int N;
            fftPlan(int N);
            ~fftPlan();
        private:
            void makePhase(Complex *omega);
            double getWeights(unsigned int k);  // Function to add weights

        };

        namespace scalar
        {
            // Function declarations
            void FFT(Complex *Fout, Complex *Fin, fft::cpu::fftPlan plan, bool forward);
            void BitRevArray(Complex *Frev, Complex *F, int N);
            unsigned int reverseBits(unsigned int x, unsigned int N);
        }
    }
}
