CC = g++
CFLAGS = -Wall -O2
LDFLAGS = -lfftw3

SRCS = fft_sparse.cpp dft_3d_fft.cpp dft_3d_fftw.cpp dft_3d_naive.cpp
OBJS = $(SRCS:.cpp=.o)
EXEC = fft_sparse.out

all: $(EXEC)

$(EXEC): $(OBJS)
	$(CC) $(CFLAGS) $(OBJS) -o $(EXEC) $(LDFLAGS)

fft_sparse.o: fft_sparse.cpp fft_sparse.h
	$(CC) $(CFLAGS) -c fft_sparse.cpp -o fft_sparse.o

clean:
	rm -f $(OBJS) $(EXEC)