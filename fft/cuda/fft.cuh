#pragma once

#include <cmath>
#include <cuComplex.h>
#include <cuda_runtime.h>


namespace fft
{
	namespace cuda
	{
		typedef cuDoubleComplex Complex;

		const Complex I = {0, 1};

		class fftPlan
		{
		public:
			Complex *omega;
			unsigned int N;
			fftPlan(unsigned int N);
			~fftPlan();
		private:
			void makePhase(Complex *omega, unsigned int N);
		};


		namespace base
		{
			// Function declarations
			void FFT(Complex *Ftilde, Complex *F, Complex *omega, int N);
			void BitRevArray(Complex *Frev, Complex *F, int N);
			unsigned int reverseBits(unsigned int x, unsigned int N);
			__global__ void makePhase(Complex *omega, unsigned int N);
		}
	}
}