#include "fft.h"

#define N (int)pow(2, 6)

int main(void) {
    Complex * f = new Complex[N];
    Complex * f_tilde = new Complex[N];

    // Prepare a cosine wave.
    for (int i = 0; i < N; i++) {
        f[i] = Complex(cos(3 * 2 * M_PI * i / N), 0.0);
        f_tilde[i] = Complex(0.0, 0.0);
    }

    // Compare.
    compare_all(N, 16 * 2, f, f_tilde, false);
}