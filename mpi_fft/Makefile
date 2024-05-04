# Compiler settings
CC=mpicxx
CFLAGS=

# Target executable name
TARGET=fft_mpi

# Source files
SOURCES=fft.cpp fft_scalar.cpp
OBJECTS=$(SOURCES:.cpp=.o)

all: $(TARGET)

$(TARGET): $(OBJECTS)
	$(CC) $(CFLAGS) $^ -o $@

%.o: %.cpp
	$(CC) $(CFLAGS) -c $< -o $@

clean:
	rm -f $(OBJECTS) $(TARGET)

.PHONY: all clean
