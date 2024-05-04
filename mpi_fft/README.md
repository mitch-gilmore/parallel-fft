# Fast Fourier Transform (FFT) with OpenACC

## Overview
This project implements a Fast Fourier Transform (FFT) using C++ and enhances its performance using OpenACC directives for parallel computing on accelerators like GPUs. It is designed to demonstrate the use of complex numbers and the efficiency improvements that can be achieved through parallel processing.

## Requirements
- PGI Compiler with OpenACC support, or another OpenACC-capable compiler like GCC.
- A system with an NVIDIA GPU is recommended for best performance when using OpenACC.

## File Descriptions
- `fft.h`: Header file containing the declaration of the `fftPlan` class and related functions.
- `fft.cpp`: Implementation of the `fftPlan` class, including methods for calculating the phase vector using FFT.
- `fft_scalar.cpp`: Contains the scalar implementation of the FFT and related utility functions.
- `Makefile`: Makefile for building the project with appropriate flags for OpenACC.
- `run_tests.sh`: Shell script for compiling and running the project, and saving the output.

## Setup
Clone the repository or download the source files to a local directory.

## Compilation
Use the provided Makefile to compile the project. Run the following command in the terminal:

```bash
make
```

This command compiles the code into an executable named `fft_acc`. It uses OpenACC flags to enable GPU acceleration.

## Running Tests
To run the FFT tests and save the results, execute the `run_tests.sh` script by following these steps:

1. Make the script executable:
   ```bash
   chmod +x run_tests.sh
   ```

2. Run the script:
   ```bash
   ./run_tests.sh
   ```

The script compiles the code (if not already compiled) and runs the executable, saving the output to `fft_results.txt`.

## Viewing Results
After running the tests, check the results in the `fft_results.txt` file:

```bash
cat fft_results.txt
```

This file will contain the output from the FFT calculations, demonstrating the transform's performance and results.

## Clean Up
To clean up the compiled objects and executable, run:

```bash
make clean
```

This will remove all generated files, leaving only the source code.
