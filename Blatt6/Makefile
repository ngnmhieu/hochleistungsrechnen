# Common definitions
CC = mpicc

# Compiler flags, paths and libraries
CFLAGS = -pedantic -Wall -Wextra -O3 -ggdb
LFLAGS = $(CFLAGS)

all: timempi timempi2

timempi: timempi.o
	$(CC) $(LFLAGS) -o $@ timempi.o

timempi2: timempi2.o
	$(CC) $(LFLAGS) -o $@ timempi2.o

%.o: %.c
	$(CC) -c $(CFLAGS) $*.c

clean:
	rm -f *.o 
	rm -f timempi timempi2
