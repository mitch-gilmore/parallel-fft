#include <stdio.h>
#include <complex.h>
#include <fftw3.h>
#include <pthread.h>
#include <math.h>

// gcc -o dft_3d_fftw_threads.out dft_3d_fftw_threads.cpp -lfftw3 -lpthread

#define NX 32
#define NY 32
#define NZ 32

#define SIZE NX * NY * NZ

#define NUM_THREADS 4

typedef std::complex<double> Complex;

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

// Global variables
fftw_complex in[SIZE], out[SIZE];
fftw_plan plan;
pthread_barrier_t barrier;

// Worker function for FFT computation
void* compute_fft(void* arg) {
    int id = *(int*)arg;

    int start = id * (SIZE / NUM_THREADS);
    int end = (id + 1) * (SIZE / NUM_THREADS);

    // Perform FFT along the first dimension
    fftw_execute_dft(plan, &in[start * NX * NY], &out[start * NX * NY]);

    // Wait for all threads to finish FFT computation
    pthread_barrier_wait(&barrier);

    // Perform FFT along the second dimension
    for (int i = start; i < end; i++) {
        fftw_execute_dft(plan, &out[i * NX * NY], &out[i * NX * NY]);
    }

    // Wait for all threads to finish FFT computation
    pthread_barrier_wait(&barrier);

    // Perform FFT along the third dimension
    for (int i = start; i < end; i++) {
        fftw_execute_dft(plan, &out[i * NX * NY], &out[i * NX * NY]);
    }

    return NULL;
}

int main() {
    pthread_t threads[NUM_THREADS];
    int thread_ids[NUM_THREADS];

    // Initialize FFTW
    plan = fftw_plan_dft_3d(NX, NY, NZ, in, out, FFTW_FORWARD, FFTW_ESTIMATE);

    // Initialize barrier
    pthread_barrier_init(&barrier, NULL, NUM_THREADS);

    // Create threads
    for (int i = 0; i < NUM_THREADS; i++) {
        thread_ids[i] = i;
        pthread_create(&threads[i], NULL, compute_fft, &thread_ids[i]);
    }

    // Join threads
    for (int i = 0; i < NUM_THREADS; i++) {
        pthread_join(threads[i], NULL);
    }

    // Destroy barrier
    pthread_barrier_destroy(&barrier);

    // Clean up
    fftw_destroy_plan(plan);
    fftw_free(in);
    fftw_free(out);

    return 0;
}
