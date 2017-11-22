# Common definitions
CC = gcc

# Compiler flags, paths and libraries
CFLAGS = -std=c11 -pedantic -Wall -Wextra -O3 -ggdb
LFLAGS = $(CFLAGS)
PX_CFLAGS = $(CFLAGS) -pthread
PX_LFLAGS = $(PX_CFLAGS)
LIBS   = -lm

TARGET_SEQ = partdiff-seq
TARGET_PX = partdiff-posix

OBJ_SEQ = partdiff-seq.o askparams.o
OBJ_PX = partdiff-posix.o askparams.o

all: partdiff-seq partdiff-posix

askparams.o:
	$(CC) -c $(CFLAGS) askparams.c

partdiff-seq: $(OBJ_SEQ)
		$(CC) $(LFLAGS) -o $@ $(OBJ_SEQ) $(LIBS)

partdiff-seq.o:
	$(CC) -c $(CFLAGS) partdiff-seq.c

partdiff-posix: $(OBJ_PX)
		$(CC) $(PX_LFLAGS) -o $@ $(OBJ_PX) $(LIBS)

partdiff-posix.o:
	$(CC) -c $(PX_CFLAGS) partdiff-posix.c 

# Rule to create *.o from *.c
%.o: %.c
	$(CC) -c $(CFLAGS) $*.c

clean:
	$(RM) $(OBJ_SEQ) $(OBJ_PX)
	$(RM) $(TARGET_SEQ) $(TARGET_PX)
