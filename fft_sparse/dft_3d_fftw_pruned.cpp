#include "fft_sparse.h"
#include "fftw3.h"

#define FFTW_TO_COMPLEX(num) Complex(num[0], num[1])
#define COMPLEX_TO_FFTW(cmplx, num) \
    num[0] = cmplx.real(); \
    num[1] = cmplx.imag();
    
double _fft_pruned_1d(int N, int K, Complex * f, Complex * f_tilde) {
  fftw_complex in[N], out[N], twids[(K-1)*(N/K-1)];

  // Precompute twiddle factors (since we usually want more than one FFT)
  for (int j = 1; j < N/K; ++j)
    for (int i = 1; i < K; ++i) {
      Complex cmp = cexp((I * FFTW_FORWARD * 2 * M_PI / N) * (i*j));
      twids[(j-1)*(K-1) + (i-1)][0] = cmp.real();
      twids[(j-1)*(K-1) + (i-1)][1] = cmp.imag();
    }

  // Plan N/K FFTs of size K.
  fftw_plan plan = fftw_plan_many_dft(1, &K, N/K, in, NULL, N/K, 1, out, NULL, 1, K, FFTW_FORWARD, FFTW_ESTIMATE);

  // Prepare input.
  for (int i = 0; i < N; i++) {
    in[i][0] = f[i].real();
    in[i][1] = f[i].imag();
  }

  // Execute.
  auto start = std::chrono::high_resolution_clock::now();

  fftw_execute(plan);

  auto end = std::chrono::high_resolution_clock::now();

  // Prepare output.
  for (int j = 1; j < N/K; ++j) {
    Complex res = FFTW_TO_COMPLEX(out[0]) + FFTW_TO_COMPLEX(out[j*K]);
    COMPLEX_TO_FFTW(res, out[0]);

    for (int i = 1; i < K; ++i) {
      Complex st = FFTW_TO_COMPLEX(out[i]) + (FFTW_TO_COMPLEX(out[i + j*K]) * FFTW_TO_COMPLEX(twids[(j-1)*(K-1) + (i-1)]));
      COMPLEX_TO_FFTW(st, out[i]);
    }
  }

  for (int i = 0; i < N; i++) {
    f_tilde[i] = Complex(out[i][0], out[i][1]);
  }

  fftw_destroy_plan(plan);

  // Return timing.
  double result;
  DURATION_MS(start, end, result);

  return result;
}