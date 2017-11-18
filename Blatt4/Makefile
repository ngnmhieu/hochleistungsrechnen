# Common definitions
CC = gcc

# Compiler flags, paths and libraries
CFLAGS = -std=c11 -pedantic -Wall -Wextra -O3 -ggdb
LFLAGS = $(CFLAGS)
LIBS   = -lm
OCFLAGS = $(CFLAGS) -fopenmp
OLFLAGS = $(OCFLAGS)

TGTS = partdiff-seq.c askparams.c
OTGTS = partdiff-openmp.c askparams.c

OBJS = partdiff-seq.o askparams.o
OOBJS = partdiff-openmp.o askparams.o

# Targets ...
all: partdiff-seq partdiff-openmp

link:
	$(CC) -c $(CFLAGS) $(TGTS)

partdiff-seq: link Makefile
	$(CC) $(LFLAGS) -o $@ $(OBJS) $(LIBS)

olink:
	$(CC) -c $(OCFLAGS) $(OTGTS)

partdiff-openmp: olink Makefile
	$(CC) $(OLFLAGS) -o $@ $(OOBJS) $(LIBS)


partdiff-seq.o: partdiff-seq.c Makefile

askparams.o: askparams.c Makefile

# Rule to create *.o from *.c
%.o: %.c
	$(CC) -c $(CFLAGS) $*.c

clean:
	$(RM) $(OBJS)
	$(RM) $(TGTS)
