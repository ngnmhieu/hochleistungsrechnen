# Common definitions
CC = mpicc

# Compiler flags, paths and libraries
CFLAGS = -pedantic -Wall -Wextra -O3 -ggdb
LFLAGS = $(CFLAGS)

all: timempi

timempi: timempi.o
	$(CC) $(LFLAGS) -o $@ timempi.o

timempi.o: 
	$(CC) -c $(CFLAGS) timempi.c

clean:
	rm timempi.o timempi
