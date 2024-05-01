#pragma once

#include <cmath>
#include <math.h>
#include <cuComplex.h>


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
			unsigned int stages;
			unsigned int N;
			double phase_scale;
			fftPlan(unsigned int N, double phase_scale);
			~fftPlan();
		private:
			void makePhase();
			void makeIndexShuffle();
		};

		void FFT(Complex *Fout, Complex *Fin, fftPlan plan, bool forward);
	}
}