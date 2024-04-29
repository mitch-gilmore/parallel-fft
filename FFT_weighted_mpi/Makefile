CXX = mpic++
CXXFLAGS = -std=c++11 -Wall
LDFLAGS = -lm

TARGET = fft_mpi_test
SOURCES = mainFFT_weighted_mpi.cpp fft_weighted.cpp
OBJECTS = $(SOURCES:.cpp=.o)

all: $(TARGET)

$(TARGET): $(OBJECTS)
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LDFLAGS)

%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

clean:
	rm -f $(OBJECTS) $(TARGET)

.PHONY: all clean
