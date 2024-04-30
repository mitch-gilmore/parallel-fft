#include <iostream>
#include <fstream>
#include <cmath>
#include <chrono>
#include <math.h>
#include <complex.h>

#include <fftw3.h>

#define DURATION_MS(start, end, timing) \
    std::chrono::duration<double, std::milli> duration = end - start; \
    timing = duration.count()

/// Shorthand Complex number.
typedef std::complex<double> Complex;

/// Naive 3D DFT (no parallelization)
double dft_3d_naive(int nx, int ny, int nz, Complex * in, Complex * out);

//// FFT (no parallelization)
double dft_3d_fft(int nx, int ny, int nz, Complex * in, Complex * out);

//// FFT with FFTW (no parallelization)
double dft_3d_fftw(int nx, int ny, int nz, Complex * in, Complex * out);