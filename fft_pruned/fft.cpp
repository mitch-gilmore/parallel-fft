#include "fft.h"

double ft_dft(int N, int K, Complex * f, Complex * f_tilde) {
  auto start = std::chrono::high_resolution_clock::now();

  for (int k = 0; k < N; ++k) {
    for (int n = 0; n < N; ++n) {
      f_tilde[k] += f[n] * exp(-2.0 * M_PI * Complex(0, 1) * static_cast<double>(k * n) / static_cast<double>(N));
    }
  }

  auto end = std::chrono::high_resolution_clock::now();

  double result;
  DURATION_MS(start, end, result);

  return result;
}

double ift_dft(int N, Complex * f, Complex * f_tilde) {
  auto start = std::chrono::high_resolution_clock::now();

  for (int k = 0; k < N; ++k) {
    f[k] = 0;
    for (int n = 0; n < N; ++n) {
      f[k] += f_tilde[n] * exp(2.0 * M_PI * Complex(0, 1) * static_cast<double>(k * n) / static_cast<double>(N)) * (1 / static_cast<double>(N));
    }
  }

  auto end = std::chrono::high_resolution_clock::now();

  double result;
  DURATION_MS(start, end, result);

  return result;
}

void _ft_fft_h(int N, Complex *f, Complex *f_tilde) {
  if(N <= 2) {
    f_tilde[0] = f[0] + f[1];
    f_tilde[1] = f[0] - f[1];
    return;
  }

  // Split even and odd.
  Complex *even = new Complex[N / 2];
  Complex *odd = new Complex[N / 2];

  Complex *even_tilde = new Complex[N / 2];
  Complex *odd_tilde = new Complex[N / 2];

  for (int i = 0; i < N / 2; ++i) {
    even[i] = f[2 * i];
    odd[i] = f[2 * i + 1];
  }

  _ft_fft_h(N / 2, even, even_tilde);
  _ft_fft_h(N / 2, odd, odd_tilde);

  double theta = M_PI * 2.0 / (double)N;

  for (int k = 0; k < N / 2; ++k) {
    Complex t = exp(theta * Complex(0, 1) *(double)k) * odd_tilde[k];
    f_tilde[k] = even_tilde[k] + t;
    f_tilde[k + N/2] = even_tilde[k] - t;
  }

  delete[] even;
  delete[] odd;

  delete[] even_tilde;
  delete[] odd_tilde;
}

double ft_fft(int N, int K, Complex * f, Complex * f_tidle) {
  auto start = std::chrono::high_resolution_clock::now();

  _ft_fft_h(N, f, f_tidle);

  auto end = std::chrono::high_resolution_clock::now();

  double result;
  DURATION_MS(start, end, result);

  return result;
}

double ft_fftw(int N, int K, Complex * f, Complex * f_tilde) {
  fftw_complex in[N], out[N];

  // Create plan before init input.
  fftw_plan p = fftw_plan_dft_1d(N, in, out, FFTW_FORWARD, FFTW_ESTIMATE);

  // Prepare input.
  for (int i = 0; i < N; i++) {
    in[i][0] = f[i].real();
    in[i][1] = f[i].imag();
  }

  // Execute.
  auto start = std::chrono::high_resolution_clock::now();

  fftw_execute(p);

  auto end = std::chrono::high_resolution_clock::now();

  // Prepare output.
  for (int i = 0; i < N; i++) {
    f_tilde[i] = Complex(out[i][0], out[i][1]);
  }

  fftw_destroy_plan(p);

  // Return timing.
  double result;
  DURATION_MS(start, end, result);

  return result;
}

double ft_pruned_fftw(int N, int K, Complex * f, Complex * f_tilde) {
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