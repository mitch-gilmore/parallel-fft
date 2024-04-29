#include "fft.h"

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
    Complex * truth,
    bool set_truth,
    bool print_progress
  ) {

  // Reset fftw.
  fftw_cleanup();

  // Execute FT solver.
  double timing = method(N, K, f, f_tilde);

  // Print results.
  printf("\n[%s] timing (ms): %.2f results:\n", name.c_str(), timing);

  // Check correctness.
  bool correct = true;
  for (int n = 0; n < N; ++n) {
    if (print_progress) {
        printf("[%s] freq: %3d %+9.5f %+9.5f I\n", name.c_str(), n, f_tilde[n].real() * 1./N, f_tilde[n].imag() * 1./N);
    }

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

void compare_all(int N, int K, Complex * f, Complex * f_tilde, bool print_progress) {
  FTSolver methods[] = {ft_dft, ft_fft, ft_fftw, ft_pruned_fftw};
  std::string names[] = {"Naive DFT", "Custom FFT", "FFTW FFT", "Pruned FFT FFTW"};

  Complex * truth = new Complex[N];

  for (int idx = 0; idx < 4; ++idx) {
    FTSolver method = methods[idx];
    std::string name = names[idx];

    run_method(method, name, N, K, f, f_tilde, truth, idx == 0, print_progress);
  }

  delete [] truth;
}