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

#define FFTW_TO_COMPLEX(num) Complex(num[0], num[1])
#define COMPLEX_TO_FFTW(cmplx, num) \
    num[0] = cmplx.real(); \
    num[1] = cmplx.imag();

/// Shorthand Complex number.
typedef std::complex<double> Complex;

/// Fourier transform solver method definition.
typedef double (*FTSolver)(int, int, Complex *, Complex *);

// Naive dft implementation.
double ft_dft(int N, int K, Complex * f, Complex * f_tilde);
double ift_dft(int N, Complex * f, Complex * f_tilde);

// Custom fft implementation.
double ft_fft(int N, int K, Complex * f, Complex * f_tidle);

// Using fftw.
double ft_fftw(int N, int K, Complex * f, Complex * f_tilde);

// Pruned variants.
double ft_pruned_fftw(int N, int K, Complex * f, Complex * f_tilde);

// Testing.
void run_method(
    FTSolver method, 
    std::string name,
    int N, 
    int K,
    Complex * f, 
    Complex * f_tilde,
    Complex * truth = NULL,
    bool set_truth = false,
    bool print_progress = true
);

void compare_all(
    int N, 
    int K, 
    Complex * f, 
    Complex * f_tilde, 
    bool print_progress = true
);
