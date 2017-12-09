# Common definitions
CC = mpicc

# Compiler flags, paths and libraries
CFLAGS = -pedantic -Wall -Wextra -O3 -ggdb
LFLAGS = $(CFLAGS)

all: circle

circle: circle.o
	$(CC) $(LFLAGS) -o $@ circle.o

%.o: %.c
	$(CC) -c $(CFLAGS) $*.c

clean:
	$(RM) *.o
	$(RM) circle
