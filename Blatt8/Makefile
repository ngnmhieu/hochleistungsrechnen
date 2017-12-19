# Common definitions
CC = mpicc

# Compiler flags, paths and libraries
CFLAGS = -std=c11 -pedantic -Wall -Wextra
LFLAGS = $(CFLAGS)
LIBS   = -lm

OBJS = partdiff-par.o askparams.o
OBJS_SEQ = partdiff-seq.o askparams.o

# Rule to create *.o from *.c
.c.o:
	$(CC) -c $(CFLAGS) $*.c

# Targets ...
all: partdiff-par partdiff-seq

partdiff-seq: $(OBJS_SEQ) Makefile
	gcc $(LFLAGS) -o $@ $(OBJS_SEQ) $(LIBS)

partdiff-par: $(OBJS) Makefile
	$(CC) $(LFLAGS) -o $@ $(OBJS) $(LIBS)

vampir:
	vtcc -vt:cc $(CC) $(CFLAGS) $(OBJS) -o circle -lm

clean:
	rm -f partdiff-seq partdiff-par *.z *.thumb *.o *.otf

test: partdiff-par partdiff-seq
	./test.sh
