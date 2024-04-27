#pragma once

#include <cmath>
#include <complex>

#define PI M_PI
#define I Complex(0.0, 1.0)

// using namespace std;
typedef std::complex<double> Complex;

namespace fft
{
	namespace cpu
	{

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
			unsigned int reverseBits(unsigned int x, int N);
		}
	}
}