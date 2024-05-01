# Makefile for FFT with weights
CXX = g++
CXXFLAGS = -std=c++11 -Wall
LDFLAGS = -lm  # Link against math library

GXX = nvcc

C_SOURCES = $(shell find fft/cpu -type f -name '*.cpp')
C_OBJECTS = $(C_SOURCES:.cpp=.o)

CUDA_SOURCES = $(shell find fft/cuda -type f -name '*.cu')
CUDA_OBJECTS = $(CUDA_SOURCES:.cu=.o)

OBJECTS = $(C_OBJECTS) $(CUDA_OBJECTS)

CUDA_TARGET = fft_cuda_time

TARGETS = $(CUDA_TARGET)

# Default target
all: $(OBJECTS) $(TARGETS)


# Compile each source file to an object file
%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

%.o: %.cu
	$(GXX) -c $< -o $@

fft_cuda_time: fft/main.cu $(CUDA_OBJECTS)
	$(GXX) -o $@ $^ $(LDFLAGS)

# Clean up build files
clean:
	rm -f $(shell find . -type f -name '*.o') $(TARGETS)

.PHONY: all clean
