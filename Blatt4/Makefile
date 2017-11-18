# Common definitions
CC = gcc

# Compiler flags, paths and libraries
CFLAGS = -std=c11 -pedantic -Wall -Wextra -O3 -ggdb
LFLAGS = $(CFLAGS)
LIBS   = -lm
OCFLAGS = $(CFLAGS) -fopenmp
OLFLAGS = $(OCFLAGS)

TGTS = partdiff-seq.c askparams.c
TGTS2 = partdiff-openmp.c askparams.c

OBJS = partdiff-seq.o askparams.o
OBJS2 = partdiff-openmp.o askparams.o

# Targets ...
all: partdiff-seq partdiff-openmp

add:
	$(CC) -c $(CFLAGS) $(TGTS)

partdiff-seq: add Makefile
	$(CC) $(LFLAGS) -o $@ $(OBJS) $(LIBS)

add2:
	$(CC) -c $(OCFLAGS) $(TGTS2)

partdiff-openmp: add2 Makefile
	$(CC) $(OLFLAGS) -o $@ $(OBJS2) $(LIBS)


partdiff-seq.o: partdiff-seq.c Makefile

askparams.o: askparams.c Makefile

# Rule to create *.o from *.c
%.o: %.c
	$(CC) -c $(CFLAGS) $*.c

clean:
	$(RM) $(OBJS)
	$(RM) $(TGTS)
