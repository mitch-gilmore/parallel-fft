#include "fft_sparse.h"

#define NX (int)pow(2, 5)
#define NY (int)pow(2, 5)
#define NZ (int)pow(2, 5)

#define SIZE NX * NY * NZ

void copy(int size, Complex * copy_from, Complex * copy_to) {
    for (int i = 0; i < size; i++) {
        copy_to[i] = copy_from[i];
    }
}

int main(void) {
    int dims[] = {NX, NY, NZ};

    Complex * in = new Complex[SIZE];
    Complex * out = new Complex[SIZE];

    Complex * in_cp = new Complex[SIZE];
    Complex * out_cp = new Complex[SIZE];

    generate_data(NX, NY, NZ, in);

    double time;

    copy(SIZE, in, in_cp);
    copy(SIZE, out, out_cp);
    
    printf("\nInitializing...\n");

    time = dft_3d_naive(NX, NY, NZ, in_cp, out_cp);
    printf("[DFT 3D NAIVE] time ms: %.8f \n", time);

    copy(SIZE, in, in_cp);
    copy(SIZE, out, out_cp);

    time = dft_3d_fft(NX, NY, NZ, in_cp, out_cp);
    printf("[DFT 3D CUSTOM FFT] time ms: %.8f \n", time);

    copy(SIZE, in, in_cp);
    copy(SIZE, out, out_cp);

    time = dft_3d_fftw(NX, NY, NZ, in_cp, out_cp);
    printf("[DFT 3D FFTW] time ms: %.8f \n", time);

    printf("Finished...\n");

    delete [] in;
    delete [] out;

    delete [] in_cp;
    delete [] out_cp;
}