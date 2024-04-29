#include "fft.cuh"

namespace fft
{
	namespace cuda
	{
		fftPlan::fftPlan(unsigned int N)
		{
			this->N = N;
			cudaMalloc((void**)&this->omega, N * sizeof(Complex));

			this->makePhase(this->omega, this->N);
		}

		fftPlan::~fftPlan()
		{
			delete[] omega;
		}

		void fftPlan::makePhase(Complex *omega, unsigned int N)
		{
			int threadsPerBlock = 32;
			int numBlocks = (this->N + threadsPerBlock - 1) / threadsPerBlock;

			base::makePhase<<<numBlocks, threadsPerBlock>>>(omega, N);
		}

		namespace base
		{
			using fft::cuda::Complex;
			using fft::cuda::I;
			using fft::cuda::fftPlan;

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

			void FFT(Complex *Fout, Complex *Fin, fftPlan plan, bool forward)
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

							Complex phase = (forward) ? plan.omega[Nblocks * pairs] : cuConj(plan.omega[Nblocks * pairs]);
							Complex low = Fout[k];
							Complex high = cuCmul(phase, Fout[k + Bsize / 2]);
							Fout[k] = cuCadd(low, high);
							Fout[k + Bsize / 2] = cuCsub(low, high);
						}
					}
				}
			}


			__global__ void makePhase(Complex *omega, unsigned int N)
			{
				int idx = threadIdx.x + blockIdx.x * blockDim.x;
	
				if (idx < N)
				{
					double phase = 2.0 * M_PI * (double)idx / (double)N;
					omega[idx].x = cos(phase);
					omega[idx].y = sin(phase);
				}
			}
		}
	}
}