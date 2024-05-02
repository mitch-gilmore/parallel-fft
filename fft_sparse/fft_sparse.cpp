#include "fft_sparse.h"

#define NX (int)pow(2, 5)
#define NY (int)pow(2, 5)
#define NZ (int)pow(2, 5)

#define SIZE NX * NY * NZ

/// Generate cos data in 3 dimensions.
void generate_data(int nx, int ny, int nz, Complex * in) {
    int dims[] = {nx, ny, nz};

    for (int i = 0; i < 3; i++) {
        int dim = dims[i];
        int pos = i + dim;

        // Prepare a cosine wave.
        for (int i = 0; i < dim; i++) {
            in[pos + i] = Complex(cos(3 * 2 * M_PI * i / dim), 0.0);
        }
    }
}

void copy(int size, Complex * copy_from, Complex * copy_to) {
    for (int i = 0; i < size; i++) {
        copy_to[i] = copy_from[i];
    }
}

void find_largest(int size, Complex * out) {
    double magnitude = 0.0, largest = 0.0;
    int lg_index = 0;

    for (int i = 0; i < size; i++) {
        magnitude = sqrt(pow(out[i].real(),2)+pow(out[i].imag(),2));
        if (magnitude > largest) {
            largest = magnitude;
            lg_index = i;
        }
    }

    printf("Found magnitude %.8f at index %d\n\n", largest, lg_index);
}

int main(void) {
    double time;

    Complex * in = new Complex[SIZE];
    Complex * out = new Complex[SIZE];

    Complex * in_cp = new Complex[SIZE];
    Complex * out_cp = new Complex[SIZE];

    printf("\nInitializing 3D test data...\n");

    generate_data(NX, NY, NZ, in);
    copy(SIZE, in, in_cp);

    time = dft_3d_naive(NX, NY, NZ, in_cp, out_cp);
    printf("[DFT 3D NAIVE] time ms: %.8f \n", time);

    copy(SIZE, in, in_cp);

    time = dft_3d_fft(NX, NY, NZ, in_cp, out_cp);
    printf("[FFT 3D NAIVE] time ms: %.8f \n", time);

    copy(SIZE, in, in_cp);

    time = dft_3d_fft_threads(NX, NY, NZ, in_cp, out_cp);
    printf("[FFT 3D THREADS] time ms: %.8f \n", time);

    copy(SIZE, in, in_cp);

    time = dft_3d_fftw(NX, NY, NZ, in_cp, out_cp);
    printf("[FFT 3D FFTW] time ms: %.8f \n", time);

    printf("[FFT 3D SPARSE] time ms: IN PROGRESS (see fft_sparse.cpp and references)\n");

    printf("Finished...\n");

    delete [] in;
    delete [] out;

    delete [] in_cp;
    delete [] out_cp;
}
