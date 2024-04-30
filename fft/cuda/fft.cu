#include "fft.cuh"

namespace fft
{
	namespace cuda
	{

		// Constructor
		fftPlan::fftPlan(unsigned int N)
		{
			this->N = N;

			cudaMalloc((void**)&this->omega, N * sizeof(Complex));
			this->makePhase();

			cudaMalloc((void**)&this->indexShuffle, N * sizeof(int));
			this->makeIndexShuffle();
		}

		// Destructor
		fftPlan::~fftPlan()
		{
			cudaFree(this->omega);
			cudaFree(this->indexShuffle);
		}

		// Computes all the complex roots of unity
		void fftPlan::makePhase()
		{
			int threadsPerBlock = 32;
			int numBlocks = (this->N + threadsPerBlock - 1) / threadsPerBlock;

			base::makePhase<<<numBlocks, threadsPerBlock>>>(this->omega, this->N);
		}

		// Computes the index shuffle
		void fftPlan::makeIndexShuffle()
		{
			int threadsPerBlock = 32;
			int numBlocks = (this->N + threadsPerBlock - 1) / threadsPerBlock;

			base::makeIndexShuffle<<<numBlocks, threadsPerBlock>>>(this->indexShuffle, this->indexShuffle, this->N);
		}

		namespace base
		{
			using fft::cuda::Complex;
			using fft::cuda::I;
			using fft::cuda::fftPlan;
			
			__global__ static void FFT(Complex *Fout, Complex *Fin, int *indeces, unsigned int N, unsigned int stages, bool forward)
			{
				int idx = threadIdx.x + blockIdx.x * blockDim.x;

				complex twiddle = (forward) ? plan.omega[N / 2] : cuConj(plan.omega[N / 2]);

				if (idx < N/2)
				{
					// stockhams shuffle
					{
						int Bsize = 2;
						int Nblocks = (int)N / Bsize;

						int nb = idx / Bsize;
						int pairs = idx % Bsize / 2;
						int k = pairs + nb * Bsize;

						Complex phase = (forward) ? plan.omega[Nblocks * pairs] : cuConj(plan.omega[Nblocks * pairs]);
						Complex low = Fin[k];
						Complex high = cuCmul(phase, Fin[k + Bsize / 2]);

						Fout[indeces[k]] = cuCadd(low, high);
						Fout[indeces[k + Bsize / 2]] = cuCadd(low, cuCmul(twiddle, high));

						__syncthreads();

					}

					for (int s = 1; s < stages; s++)
					{
						int Bsize = pow(2, s + 1);
						int Nblocks = N / Bsize;

						int nb = idx / Bsize;
						int pairs = idx % Bsize / 2;
						int k = pairs + nb * Bsize;
								
						Complex phase = (forward) ? plan.omega[Nblocks * pairs] : cuConj(plan.omega[Nblocks * pairs]);
						Complex low = Fout[k];
						Complex high = cuCmul(phase, Fout[k + Bsize / 2]);
						Fout[k] = cuCadd(low, high);
						Fout[k + Bsize / 2] = cuCadd(low, cuCmul(twiddle, high));

						__syncthreads();
					}
				}
			}

			void FFT(Complex *Fout, Complex *Fin, fftPlan plan, bool forward)
			{
				int log2n = log2(plan.N);
				
				int threadsPerBlock = 32;
				int numBlocks = (plan.N + threadsPerBlock - 1) / (threadsPerBlock * 2);

				FFT<<<numBlocks, threadsPerBlock>>>(Fout, Fin, plan.indexShuffle, plan.N, log2n, true);

			}
			
			__global__ void BitRevArray(int *indeces, int N)
			{
				int idx = threadIdx.x + blockIdx.x * blockDim.x;
				if (idx < N)
					indeces[idx] = reverseBits(idx, N);
			}

			__forceinline __device__ unsigned int reverseBits(unsigned int x, unsigned int N)
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
			
			__global__ void makePhase(Complex *omega, unsigned int N)
			{
				int idx = threadIdx.x + blockIdx.x * blockDim.x;
	
				if (idx < N)
				{
					double phase = -1.0 * 2.0 * M_PI * (double)idx / (double)N;
					omega[idx].x = cos(phase);
					omega[idx].y = sin(phase);
				}
			}
		}
	}
}