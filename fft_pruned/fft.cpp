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
  if (N <= 1) {
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

  double theta = M_PI * 2 / N;
  for (int k = 0; k < N / 2; ++k) {
    Complex t = Complex(cos(k * theta), sin(k * theta)) * odd_tilde[k];
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

  // Prepare input.
  for (int i = 0; i < N; i++) {
    in[i][0] = f[i].real();
    in[i][1] = f[i].imag();
  }

  auto start = std::chrono::high_resolution_clock::now();

  fftw_plan p = fftw_plan_dft_1d(N, in, out, FFTW_FORWARD, FFTW_ESTIMATE);
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
  int i, j;
  fftw_complex in[N], out[N];

  // Prepare input.
  for (int i = 0; i < N; i++) {
    in[i][0] = f[i].real();
    in[i][1] = f[i].imag();
  }

  auto start = std::chrono::high_resolution_clock::now();

  // Plan N/K FFTs of size K.
  fftw_plan plan = fftw_plan_many_dft(1, &K, N/K, in, NULL, N/K, 1, out, NULL, 1, K, FFTW_FORWARD, FFTW_ESTIMATE);
  fftw_execute(plan);

  auto end = std::chrono::high_resolution_clock::now();

  // Prepare output.
  for (int i = 0; i < N; i++) {
    f_tilde[i] = Complex(out[i][0], out[i][1]);
  }

  fftw_destroy_plan(plan);

  // Return timing.
  double result;
  DURATION_MS(start, end, result);

  return result;
}

bool are_equal(double a, double b, double epsilon = 1e-9) {
    return std::abs(a - b) < epsilon;
}

void run_method(
    FTSolver method, 
    std::string name,
    int N, 
    int K,
    Complex * f, 
    Complex * f_tilde,
    Complex * truth = NULL,
    bool set_truth = false
  ) {

  // Execute FT solver.
  double timing = method(N, K, f, f_tilde);

  // Print results.
  printf("\n[%s] timing (ms): %.2f results:\n", name.c_str(), timing);

  // Check correctness.
  bool correct = true;
  for (int n = 0; n < N; ++n) {
    printf("[%s] freq: %3d %+9.5f %+9.5f I\n", name.c_str(), n, f_tilde[n].real() * 1./N, f_tilde[n].imag() * 1./N);

    if (truth == NULL) {
      continue;
    }

    if (set_truth) {
      truth[n] = f_tilde[n];
      continue;
    }

    bool real_eq = are_equal(truth[n].real(), f_tilde[n].real());
    bool imag_eq = are_equal(truth[n].imag(), f_tilde[n].imag());

    if (real_eq && imag_eq) {
      continue;
    } else {
      correct = false;
      break;
    }
  }

  if (truth != NULL) {
    printf("[%s] is correct: %d\n\n", name.c_str(), correct);
  } else {
    printf("\n");
  }
}

void compare_all(int N, int K, Complex * f, Complex * f_tilde) {
  FTSolver methods[] = {ft_dft, ft_fft, ft_fftw, ft_pruned_fftw};
  std::string names[] = {"Naive DFT", "Custom FFT", "FFTW FFT", "Pruned FFT FFTW"};

  Complex * truth = new Complex[N];

  for (int idx = 0; idx < 4; ++idx) {
    FTSolver method = methods[idx];
    std::string name = names[idx];

    run_method(method, name, N, K, f, f_tilde, truth, idx == 0);
  }

  delete [] truth;
}