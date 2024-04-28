#include "fft.h"

namespace fft
{
	namespace cpu
	{
		namespace scalar
		{
			using fft::cpu::Complex;
			using fft::cpu::I;

			void BitRevArray(Complex *Frev, Complex *F, int N)
			{
				for (int x = 0; x < N; x++)
					F[reverseBits(x, N)] = F[x];
			}

			unsigned int reverseBits(unsigned int x, unsigned int N)
			{
				unsigned int xrev = 0;
				unsigned int p = log2(N); // p = 4
				unsigned int n;
				unsigned int power = N;
				for (unsigned int i = 0; i < p; i++)
				{
					n = x % 2; // find lowest bit
					power /= 2;
					xrev += n * power; //  add to highest 2^3
					x /= 2;
				}
				return xrev;
			}

			void FFT(Complex *Fout, Complex *Fin, fft::cpu::fftPlan plan, bool forward)
			{
				int log2n = log2(plan.N);
				int Bsize = 0;
				int Nblocks = 0;

				for (unsigned int x = 0; x < plan.N; x++)
					Fout[reverseBits(x, plan.N)] = Fin[x];

				for (int s = 0; s < log2n; s++)
				{
					Bsize = pow(2, s + 1);
					Nblocks = plan.N / Bsize;

					for (int nb = 0; nb < Nblocks; nb++) {
						for (int pairs = 0; pairs < Bsize / 2; pairs++)
						{
							int k = pairs + nb * Bsize;

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