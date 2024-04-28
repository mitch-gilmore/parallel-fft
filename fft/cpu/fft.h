#pragma once

#include <cmath>
#include <complex>

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
		};


		namespace scalar
		{
			// Function declarations
			void FFT(Complex *Ftilde, Complex *F, Complex *omega, int N);
			void BitRevArray(Complex *Frev, Complex *F, int N);
			unsigned int reverseBits(unsigned int x, unsigned int N);
		}
	}
}