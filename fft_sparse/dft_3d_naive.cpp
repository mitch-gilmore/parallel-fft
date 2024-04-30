#include "fft_sparse.h"

void dft_1d(int dim, int N, Complex * in, Complex * out) {
  int pos = dim * N;
  for (int k = 0; k < N; ++k) {
    for (int n = 0; n < N; ++n) {
      out[pos + k] += in[pos + n] * exp(-2.0 * M_PI * Complex(0, 1) * static_cast<double>(k * n) / static_cast<double>(N));
    }
  }
}

double dft_3d_naive(int nx, int ny, int nz, Complex * in, Complex * out) {
    auto start = std::chrono::high_resolution_clock::now();

    dft_1d(0, nx, in, out);
    dft_1d(1, ny, in, out);
    dft_1d(2, nz, in, out);

    auto end = std::chrono::high_resolution_clock::now();

    double result;
    DURATION_MS(start, end, result);

    return result;
}