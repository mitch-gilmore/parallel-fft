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
			int *indexShuffle;
			unsigned int N;
			fftPlan(unsigned int N);
			~fftPlan();
		private:
			void makePhase();
			void makeIndexShuffle();
		};


		namespace base
		{
			// Function declarations
			void FFT(Complex *Ftilde, Complex *F, Complex *omega, int N);
			__device__ void BitRevArray(Complex *Frev, Complex *F, int N);
			__device__ unsigned int reverseBits(unsigned int x, unsigned int N);
			__global__ void makePhase(Complex *omega, unsigned int N);
		}
	}
}