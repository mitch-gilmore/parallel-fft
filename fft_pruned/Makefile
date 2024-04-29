CC = g++
CFLAGS = -Wall -O2
LDFLAGS = -lfftw3

SRCS = main.cpp fft.cpp compare.cpp
OBJS = $(SRCS:.cpp=.o)
EXEC = main

all: $(EXEC)

$(EXEC): $(OBJS)
	$(CC) $(CFLAGS) $(OBJS) -o $(EXEC) $(LDFLAGS)

main.o: main.cpp fft.h
	$(CC) $(CFLAGS) -c main.cpp -o main.o

clean:
	rm -f $(OBJS) $(EXEC)