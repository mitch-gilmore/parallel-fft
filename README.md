# parallel-fft

## Running basic FFT functions & timing

Code exists in `./fft_sparse` directory.

Run the following from the repo root directory to execute:

```
cd ./fft_sparse && make && ./fft_sparse.out && cd ..
```

## Running DMSFT example algorithm

This is a modified example Sparse FFT algorithm designed by markiwen08. The parameters have been modified.

Code exists in `./fft_sparse` directory.

Run the following from the repo root directory to execute:

```
cd ./fft_sparse/dmsft && make && ./dmsft_ex.out && cd ../../
```

See existing research on sparse FFTs in the `./literature` directory.


## Cleaning the repo

Run the following from the repo root directory to clean the repo.

```
cd ./fft_sparse && make clean && cd ./dmsft && make clean && cd ../../
```

## Scaling FFT

The scaling FFT was implement in cuda and can be found in [here](fft/cuda/fft.cu). In order to compile the code first ensure you have the cuda tool kit installed along with GNU make and gcc. To compile the code into a test timing code simply run the following.

```bash
make
./fft_cuda_time
```
