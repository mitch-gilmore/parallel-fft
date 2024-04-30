#include "fft.cuh"

namespace fft
{
	namespace cuda
	{	
		namespace kernels
		{
			using fft::cuda::Complex;
			using fft::cuda::I;
			using fft::cuda::fftPlan;
			
			__global__ static void FFT(Complex *Fout, Complex *Fin, Complex *omega, int *indeces, unsigned int N, unsigned int stages, bool forward)
			{
				int idx = threadIdx.x + blockIdx.x * blockDim.x;
				
				Complex twiddle = (forward) ? omega[N / 2] : cuConj(omega[N / 2]);
				
				if (idx < N/2)
				{
					// stockhams shuffle
					{
						int Bsize = 2;
						int Nblocks = (int)N / Bsize;
						
						int nb = idx / Bsize;
						int pairs = idx % Bsize / 2;
						int k = pairs + nb * Bsize;
						
						Complex phase = (forward) ? omega[Nblocks * pairs] : cuConj(omega[Nblocks * pairs]);
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
						
						Complex phase = (forward) ? omega[Nblocks * pairs] : cuConj(omega[Nblocks * pairs]);
						Complex low = Fout[k];
						Complex high = cuCmul(phase, Fout[k + Bsize / 2]);
						Fout[k] = cuCadd(low, high);
						Fout[k + Bsize / 2] = cuCadd(low, cuCmul(twiddle, high));
						
						__syncthreads();
					}
				}
			}
			
			__forceinline__ __device__ static unsigned int reverseBits(unsigned int x, unsigned int stages, unsigned int N)
			{
				unsigned int xrev = 0;
				// unsigned int p = log2(N); // p = 4
				unsigned int n;
				unsigned int power = N;
				
				for (unsigned int i = 0; i < stages; i++)
				{
					n = x % 2; // find lowest bit
					power /= 2;
					xrev += n * power; //  add to highest 2^3
					x /= 2;
				}
				
				return xrev;
			}
			
			__global__ static void makeIndexShuffle(int *indeces, unsigned int stages, unsigned int N)
			{
				int idx = threadIdx.x + blockIdx.x * blockDim.x;
				if (idx < N)
					indeces[idx] = reverseBits(idx, stages, N);
			}
			
			__global__ void static makePhase(Complex *omega, unsigned int N)
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
		
		// Constructor
		fftPlan::fftPlan(unsigned int N)
		{
			this->N = N;
			this->stages = log2(N);
			
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

			kernels::makePhase<<<numBlocks, threadsPerBlock>>>(this->omega, this->N);
		}

		// Computes the index shuffle
		void fftPlan::makeIndexShuffle()
		{
			int threadsPerBlock = 32;
			int numBlocks = (this->N + threadsPerBlock - 1) / threadsPerBlock;

			kernels::makeIndexShuffle<<<numBlocks, threadsPerBlock>>>(this->indexShuffle, this->stages, this->N);
		}
		
		void FFT(Complex *Fout, Complex *Fin, fftPlan plan, bool forward)
		{
			
			int threadsPerBlock = 32;
			int numBlocks = (plan.N + threadsPerBlock - 1) / (threadsPerBlock * 2);
			
			kernels::FFT<<<numBlocks, threadsPerBlock>>>(Fout, Fin, plan.omega, plan.indexShuffle, plan.N, plan.stages, forward);
		}
	}
}