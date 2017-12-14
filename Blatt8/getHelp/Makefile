# Common definitions
CC = mpicc

# Compiler flags, paths and libraries
CFLAGS = -std=c99 -pedantic -Wall -Wextra
LFLAGS = $(CFLAGS)
LIBS   = -lm

OBJS = partdiff-mpi.o askparams.o

# Rule to create *.o from *.c
.c.o:
	$(CC) -c $(CFLAGS) $*.c

# Targets ...
all: partdiff-mpi

partdiff-par: $(OBJS) Makefile
	$(CC) $(LFLAGS) -o $@ $(OBJS) $(LIBS)

vampir:
	vtcc -vt:cc $(CC) $(CFLAGS) $(OBCJS) -o circle -lm

clean:
	rm -f partdiff-mpi *.z *.thumb *.o *.otf

