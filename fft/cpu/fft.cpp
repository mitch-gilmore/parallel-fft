#include "fft.h"

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

		void fftPlan::makePhase(Complex *omega)
		{
			for (unsigned int k = 0; k < N; k++)
				omega[k] = exp(2.0 * M_PI * I * (double)k / (double)this->N);
		}
	}
}