#include "fft_sparse.h"
#include "fftw3.h"

double dft_3d_fftw(int nx, int ny, int nz, Complex * in, Complex * out) {
  int size = nx*ny*nz;

  fftw_complex inw[size], outw[size];

  // Create plan before init input.
  fftw_plan plan = fftw_plan_dft_3d(nx, ny, nz, inw, outw, FFTW_FORWARD, FFTW_ESTIMATE);

  // Prepare input.
  for (int i = 0; i < size; i++) {
    inw[i][0] = in[i].real();
    inw[i][1] = in[i].imag();
  }

  // Execute.
  auto start = std::chrono::high_resolution_clock::now();

  fftw_execute(plan);

  auto end = std::chrono::high_resolution_clock::now();

  // Prepare output.
  for (int i = 0; i < size; i++) {
    out[i] = Complex(outw[i][0], outw[i][1]);
  }

  fftw_destroy_plan(plan);

  // Return timing.
  double result;
  DURATION_MS(start, end, result);

  return result;
}