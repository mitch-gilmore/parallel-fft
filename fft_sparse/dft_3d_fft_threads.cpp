#include <thread>
#include <cmath>
#include <chrono>

#include "fft_sparse.h"

void _fft_threads_1d(int N, Complex * in, Complex * out) {
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

  _fft_threads_1d(N / 2, even, even_tilde);
  _fft_threads_1d(N / 2, odd, odd_tilde);

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

double dft_3d_fft_threads(int nx, int ny, int nz, Complex * in, Complex * out) {
    int dims[] = {nx, ny, nz};

    auto start = std::chrono::high_resolution_clock::now();

    std::thread threads[3];

    for (int i = 0; i < 3; ++i) {
        int N = dims[i];

        threads[i] = std::thread([&, N, i] {
            Complex *in_dim = new Complex[N];
            Complex *out_dim = new Complex[N];

            int pos = i * N;
            for (int j = 0; j < N; ++j) {
                in_dim[j] = in[pos + j];
            }

            _fft_threads_1d(N, in_dim, out_dim);

            for (int j = 0; j < N; ++j) {
                out[pos + j] = out_dim[j];
            }

            delete[] in_dim;
            delete[] out_dim;
        });
    }

    for (int i = 0; i < 3; ++i) {
        threads[i].join();
    }

    auto end = std::chrono::high_resolution_clock::now();

    std::chrono::duration<double, std::milli> duration = end - start;

    return duration.count();
}