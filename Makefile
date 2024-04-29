# Makefile for FFT with weights
CXX = g++
CXXFLAGS = -std=c++11 -Wall
LDFLAGS = -lm  # Link against math library

GXX = nvcc

C_SOURCES = $(shell find fft/cpu -type f -name '*.cpp')
C_OBJECTS = $(SOURCES:.cpp=.o)

CUDA_SOURCES = $(shell find fft/cuda -type f -name '*.cu')
CUDA_OBJECTS = $(CUDA_SOURCES:.cu=.o)

OBJECTS = $(C_OBJECTS) $(CUDA_OBJECTS)

# Default target
all: $(OBJECTS)

# $(TARGET): $(OBJECTS)
# 	$(CXX) $(CXXFLAGS) -o $@ $^ $(LDFLAGS)

# Compile each source file to an object file
%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

%.o: %.cu
	$(GXX) -c $< -o $@

# Clean up build files
clean:
	rm -f $(OBJECTS) $(TARGET)

.PHONY: all clean
