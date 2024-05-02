#include "fft_sparse.h"

void _fft_1d(int N, Complex * in, Complex * out) {
  if(N <= 2) {
    out[0] = in[0] + in[1];
    out[1] = in[0] - in[1];
    return;
  }

  // Split even and odd.
  Complex *even = new Complex[N / 2];
  Complex *odd = new Complex[N / 2];

  Complex *even_tilde = new Complex[N / 2];
  Complex *odd_tilde = new Complex[N / 2];

  for (int i = 0; i < N / 2; ++i) {
    even[i] = in[2 * i];
    odd[i] = out[2 * i + 1];
  }

  _fft_1d(N / 2, even, even_tilde);
  _fft_1d(N / 2, odd, odd_tilde);

  double theta = M_PI * 2.0 / (double)N;

  for (int k = 0; k < N / 2; ++k) {
    Complex t = exp(theta * Complex(0, 1) *(double)k) * odd_tilde[k];
    out[k] = even_tilde[k] + t;
    out[k + N/2] = even_tilde[k] - t;
  }

  delete[] even;
  delete[] odd;

  delete[] even_tilde;
  delete[] odd_tilde;
}

double dft_3d_fft(int nx, int ny, int nz, Complex * in, Complex * out) {
  int dims[] = {nx, ny, nz};

  auto start = std::chrono::high_resolution_clock::now();

  for (int i = 0; i < 3; i++) {
    int N = dims[i];

    Complex *in_dim = new Complex[N];
    Complex *out_dim = new Complex[N];

    // Copy in dim (not ideal).
    int pos = i * N;
    for (int i = 0; i < N; i++) {
        in_dim[i] = in[pos + i];
    }

    _fft_1d(N, in_dim, out_dim);

    // Copy out dim (not ideal).
    for (int i = 0; i < N; i++) {
        out[pos + i] = out_dim[i];
    }
  }


  auto end = std::chrono::high_resolution_clock::now();

  double result;
  DURATION_MS(start, end, result);

  return result;
}
