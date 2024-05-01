#include "cuda/fft.cuh"
#include <cuda_runtime.h>

#include <stdio.h>

float time_cuda(int N);
void print_cuda(int N);

int main(int argc, char **argv) {

	printf("N  | time (ms)\n");

	for (int i = 2; i <= 15; i++)
	{
		printf("%d %f\n", 1 << i, time_cuda(1 << i));
	}

	return 0;
}

float time_cuda(int N)
{
	double *h_in = new double[2*N];

	for (int i = 0; i < 2*N; i++) {
		h_in[i] = (!i) ? 1 : 0; // delta function
	}

	fft::cuda::Complex *d_in, *d_out;
	cudaMalloc(&d_in, N * sizeof(fft::cuda::Complex));
	cudaMalloc(&d_out, N * sizeof(fft::cuda::Complex));
	cudaMemcpy(d_in, h_in, N * sizeof(fft::cuda::Complex), cudaMemcpyHostToDevice);

	fft::cuda::fftPlan plan(N);

	cudaEvent_t start, stop;
	cudaEventCreate(&start);
	cudaEventCreate(&stop);

	cudaEventRecord(start);
	fft::cuda::FFT(d_in, d_out, plan, true);
	cudaEventRecord(stop);

	cudaEventSynchronize(stop);
	float milliseconds = 0;
	cudaEventElapsedTime(&milliseconds, start, stop);

	delete[] h_in;
	cudaFree(d_in);
	cudaFree(d_out);

	return milliseconds;
}