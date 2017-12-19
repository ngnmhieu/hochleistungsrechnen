#!/bin/sh

# Usage: ./partdiff-par [num] [method] [lines] [func] [term] [prec/iter]
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

SINGLE_THREAD=1
GAUSS_SEIDEL=1
JACOBI=2
NACH_GENAUIGKEIT=1
NACH_ITERATION=2

#####################################################################

NUM_PROC=4
INTERLINES=0
ITERATIONS=82
FUNC=1
# SEQ_OUT=./referenz/Jacobi.f1
MPI_OUT=partdiff-par-1.txt
SEQ_OUT=partdiff-seq-1.txt
echo "\nTest 1: Prozesse = $NUM_PROC ; Interlines = $INTERLINES; Iterations = $ITERATIONS"
echo "Command: mpirun -np $NUM_PROC ./partdiff-par $SINGLE_THREAD $JACOBI $INTERLINES $FUNC $NACH_ITERATION $ITERATIONS > $MPI_OUT"
mpirun -np $NUM_PROC ./partdiff-par $SINGLE_THREAD $JACOBI $INTERLINES $FUNC $NACH_ITERATION $ITERATIONS > $MPI_OUT
./partdiff-seq $SINGLE_THREAD $JACOBI $INTERLINES $FUNC $NACH_ITERATION $ITERATIONS > $SEQ_OUT
diff $SEQ_OUT $MPI_OUT

#####################################################################

NUM_PROC=4
INTERLINES=0
ITERATIONS=85
FUNC=2
MPI_OUT=partdiff-par-2.txt
SEQ_OUT=partdiff-seq-2.txt
echo "\nTest 2: Prozesse = $NUM_PROC ; Interlines = $INTERLINES; Iterations = $ITERATIONS"
echo "Command: mpirun -np $NUM_PROC ./partdiff-par $SINGLE_THREAD $JACOBI $INTERLINES $FUNC $NACH_ITERATION $ITERATIONS > $MPI_OUT"
mpirun -np $NUM_PROC ./partdiff-par $SINGLE_THREAD $JACOBI $INTERLINES $FUNC $NACH_ITERATION $ITERATIONS > $MPI_OUT
./partdiff-seq $SINGLE_THREAD $JACOBI $INTERLINES $FUNC $NACH_ITERATION $ITERATIONS > $SEQ_OUT
diff $SEQ_OUT $MPI_OUT

#####################################################################

NUM_PROC=4
INTERLINES=100
ITERATIONS=500
FUNC=2
MPI_OUT=partdiff-par-3.txt
SEQ_OUT=partdiff-seq-3.txt
echo "\nTest 2: Prozesse = $NUM_PROC ; Interlines = $INTERLINES; Iterations = $ITERATIONS"
echo "Command: mpirun -np $NUM_PROC ./partdiff-par $SINGLE_THREAD $JACOBI $INTERLINES $FUNC $NACH_ITERATION $ITERATIONS > $MPI_OUT"
mpirun -np $NUM_PROC ./partdiff-par $SINGLE_THREAD $JACOBI $INTERLINES $FUNC $NACH_ITERATION $ITERATIONS > $MPI_OUT
./partdiff-seq $SINGLE_THREAD $JACOBI $INTERLINES $FUNC $NACH_ITERATION $ITERATIONS > $SEQ_OUT
diff $SEQ_OUT $MPI_OUT

echo "Done"
