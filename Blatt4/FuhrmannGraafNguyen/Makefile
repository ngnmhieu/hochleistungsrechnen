# Common definitions
CC = gcc

# Compiler flags, paths and libraries
CFLAGS = -std=c11 -pedantic -Wall -Wextra -O3 -ggdb
LFLAGS = $(CFLAGS)
OMP_CFLAGS = $(CFLAGS) -fopenmp
OMP_LFLAGS = $(OMP_CFLAGS)
LIBS   = -lm

TARGET_SEQ = partdiff-seq
TARGET_OMP = partdiff-openmp

OBJ_SEQ = partdiff-seq.o askparams.o
OBJ_OMP = partdiff-openmp.o askparams.o

all: partdiff-seq partdiff-openmp

askparams.o:
	$(CC) -c $(CFLAGS) askparams.c

partdiff-seq: $(OBJ_SEQ)
		$(CC) $(LFLAGS) -o $@ $(OBJ_SEQ) $(LIBS)

partdiff-seq.o:
	$(CC) -c $(CFLAGS) partdiff-seq.c

partdiff-openmp: $(OBJ_OMP)
		$(CC) $(OMP_LFLAGS) -o $@ $(OBJ_OMP) $(LIBS)

partdiff-openmp.o:
	$(CC) -c $(OMP_CFLAGS) partdiff-openmp.c 

# Rule to create *.o from *.c
%.o: %.c
	$(CC) -c $(CFLAGS) $*.c

clean:
	$(RM) $(OBJ_SEQ) $(OBJ_OMP)
	$(RM) $(TARGET_SEQ) $(TARGET_OMP)
