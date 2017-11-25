#!/bin/sh

# Usage: ./partdiff-posix [num] [method] [lines] [func] [term] [prec/iter]
#
#   - num:       number of threads (1 .. 1024)
#   - method:    calculation method (1 .. 2)
#                  1: Gauß-Seidel
#                  2: Jacobi
#   - lines:     number of interlines (0 .. 10240)
#                  matrixsize = (interlines * 8) + 9
#   - func:      interference function (1 .. 2)
#                  1: f(x,y) = 0
#                  2: f(x,y) = 2 * pi^2 * sin(pi * x) * sin(pi * y)
#   - term:      termination condition ( 1.. 2)
#                  1: sufficient precision
#                  2: number of iterations
#   - prec/iter: depending on term:
#                  precision:  1e-4 .. 1e-20
#                  iterations:    1 .. 200000

SEQ_OUT=test-seq.txt
POSIX_OUT=test-posix.txt

NUM_THREAD=8
INTERLINES=0
FUNC=1
ITERATIONS=85
echo "\nTest 1: Threads = $NUM_THREAD ; Interlines = $INTERLINES; Iterations = $ITERATIONS"
./partdiff-seq   $NUM_THREAD 2 $INTERLINES $FUNC 2 $ITERATIONS > $SEQ_OUT
./partdiff-posix $NUM_THREAD 2 $INTERLINES $FUNC 2 $ITERATIONS > $POSIX_OUT
diff $SEQ_OUT $POSIX_OUT

NUM_THREAD=12
INTERLINES=64
FUNC=2
ITERATIONS=1000
echo "\nTest 2: Threads = $NUM_THREAD ; Interlines = $INTERLINES; Iterations = $ITERATIONS"
./partdiff-seq   $NUM_THREAD 2 $INTERLINES $FUNC 2 $ITERATIONS > $SEQ_OUT
./partdiff-posix $NUM_THREAD 2 $INTERLINES $FUNC 2 $ITERATIONS > $POSIX_OUT
diff $SEQ_OUT $POSIX_OUT

NUM_THREAD=12
INTERLINES=512
FUNC=2
ITERATIONS=1000
echo "\nTest 3: Threads = $NUM_THREAD ; Interlines = $INTERLINES; Iterations = $ITERATIONS"
./partdiff-seq   $NUM_THREAD 2 $INTERLINES $FUNC 2 $ITERATIONS > $SEQ_OUT
./partdiff-posix $NUM_THREAD 2 $INTERLINES $FUNC 2 $ITERATIONS > $POSIX_OUT
diff $SEQ_OUT $POSIX_OUT

rm $SEQ_OUT $POSIX_OUT
