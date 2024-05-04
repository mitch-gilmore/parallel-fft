LIBS = -lfftw3 -lm

CXXFLAGS = -g3 -std=c++11

all: dmsft_ex.out

dmsft.o: dmsft.cpp

dmsft_ex.out: dmsft.o
	$(CXX) $^ -O3 -o $@ $(LIBS)

clean:
	rm -f *.o *.out
