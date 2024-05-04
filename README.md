# parallel-fft

## Running basic FFT functions & timing

Code exists in `./fft_sparse` directory.

Run the following from the repo root directory to execute:

```
cd ./fft_sparse && make && ./fft_sparse.out && cd ..
```

See existing research on sparse FFTs in the `./papers` directory.

## Running DMSFT example algorithm

This is a modified example Sparse FFT algorithm designed by markiwen08. The parameters have been modified.

Code exists in `./fft_sparse` directory.

Run the following from the repo root directory to execute:

```
cd ./fft_sparse/dmsft && make && ./dmsft_ex.out && cd ../../
```

## Cleaning the repo

Run the following from the repo root directory to clean the repo.

```
cd ./fft_sparse && make clean && cd ./dmsft && make clean && cd ../../
```