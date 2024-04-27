# Makefile for FFT with weights
CXX = g++
CXXFLAGS = -std=c++11 -Wall
LDFLAGS = -lm  # Link against math library

SOURCES = $(shell find fft -type f -name '*.cpp')
OBJECTS = $(SOURCES:.cpp=.o)

# Default target
all: $(OBJECTS)

# $(TARGET): $(OBJECTS)
# 	$(CXX) $(CXXFLAGS) -o $@ $^ $(LDFLAGS)

# Compile each source file to an object file
%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

# Clean up build files
clean:
	rm -f $(OBJECTS) $(TARGET)

.PHONY: all clean
