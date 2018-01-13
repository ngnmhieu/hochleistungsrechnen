#!/bin/sh
# Usage: ./partdiff-par [num] [method] [lines] [func] [term] [prec/iter]
#
#   - num:       number of threads (1 .. 1024)
#   - method:    calculation method (1 .. 2)
#                  1: Gauß-Seidel
#                  2: Jacobi
#   - lines:     number of interlines (0 .. 10240)
#                  matrixsize = (interlines * 8) + 9 #   - func:      interference function (1 .. 2)
#                  1: f(x,y) = 0
#                  2: f(x,y) = 2 * pi^2 * sin(pi * x) * sin(pi * y)
#   - term:      termination condition ( 1.. 2)
#                  1: sufficient precision
#                  2: number of iterations
#   - prec/iter: depending on term:
#                  precision:  1e-4 .. 1e-20
#                  iterations:    1 .. 200000

# Verzeichnis für Ausgabe
mkdir -p output/

SINGLE_THREAD=1
GAUSS_SEIDEL=1
JACOBI=2
NACH_GENAUIGKEIT=1
NACH_ITERATION=2

#####################################################################

NUM_PROC=1
INTERLINES=0
ITERATIONS=82
FUNC=1
MPI_OUT=output/partdiff-par-0.txt
SEQ_OUT=output/partdiff-seq-0.txt
echo "\nTest 0: Prozesse = $NUM_PROC ; Interlines = $INTERLINES; Iterations = $ITERATIONS"
echo "Command: mpirun -np $NUM_PROC ./partdiff-par $SINGLE_THREAD $GAUSS_SEIDEL $INTERLINES $FUNC $NACH_ITERATION $ITERATIONS > $MPI_OUT"
mpirun -np $NUM_PROC ./partdiff-par $SINGLE_THREAD $GAUSS_SEIDEL $INTERLINES $FUNC $NACH_ITERATION $ITERATIONS > $MPI_OUT
./partdiff-seq $SINGLE_THREAD $GAUSS_SEIDEL $INTERLINES $FUNC $NACH_ITERATION $ITERATIONS > $SEQ_OUT
diff $SEQ_OUT $MPI_OUT

#####################################################################

NUM_PROC=1
INTERLINES=0
ITERATIONS=85
FUNC=2
MPI_OUT=output/partdiff-par-1.txt
SEQ_OUT=output/partdiff-seq-1.txt
echo "\nTest 1: Prozesse = $NUM_PROC ; Interlines = $INTERLINES; Iterations = $ITERATIONS"
echo "Command: mpirun -np $NUM_PROC ./partdiff-par $SINGLE_THREAD $GAUSS_SEIDEL $INTERLINES $FUNC $NACH_ITERATION $ITERATIONS > $MPI_OUT"
mpirun -np $NUM_PROC ./partdiff-par $SINGLE_THREAD $GAUSS_SEIDEL $INTERLINES $FUNC $NACH_ITERATION $ITERATIONS > $MPI_OUT
./partdiff-seq $SINGLE_THREAD $GAUSS_SEIDEL $INTERLINES $FUNC $NACH_ITERATION $ITERATIONS > $SEQ_OUT
diff $SEQ_OUT $MPI_OUT

#####################################################################

NUM_PROC=4
INTERLINES=0
ITERATIONS=82
FUNC=1
MPI_OUT=output/partdiff-par-2.txt
SEQ_OUT=output/partdiff-seq-2.txt
echo "\nTest 2: Prozesse = $NUM_PROC ; Interlines = $INTERLINES; Iterations = $ITERATIONS"
echo "Command: mpirun -np $NUM_PROC ./partdiff-par $SINGLE_THREAD $GAUSS_SEIDEL $INTERLINES $FUNC $NACH_ITERATION $ITERATIONS > $MPI_OUT"
mpirun -np $NUM_PROC ./partdiff-par $SINGLE_THREAD $GAUSS_SEIDEL $INTERLINES $FUNC $NACH_ITERATION $ITERATIONS > $MPI_OUT
./partdiff-seq $SINGLE_THREAD $GAUSS_SEIDEL $INTERLINES $FUNC $NACH_ITERATION $ITERATIONS > $SEQ_OUT
diff $SEQ_OUT $MPI_OUT

#####################################################################

NUM_PROC=4
INTERLINES=0
ITERATIONS=85
FUNC=2
MPI_OUT=output/partdiff-par-3.txt
SEQ_OUT=output/partdiff-seq-3.txt
echo "\nTest 3: Prozesse = $NUM_PROC ; Interlines = $INTERLINES; Iterations = $ITERATIONS"
echo "Command: mpirun -np $NUM_PROC ./partdiff-par $SINGLE_THREAD $GAUSS_SEIDEL $INTERLINES $FUNC $NACH_ITERATION $ITERATIONS > $MPI_OUT"
mpirun -np $NUM_PROC ./partdiff-par $SINGLE_THREAD $GAUSS_SEIDEL $INTERLINES $FUNC $NACH_ITERATION $ITERATIONS > $MPI_OUT
./partdiff-seq $SINGLE_THREAD $GAUSS_SEIDEL $INTERLINES $FUNC $NACH_ITERATION $ITERATIONS > $SEQ_OUT
diff $SEQ_OUT $MPI_OUT

#####################################################################

NUM_PROC=4
INTERLINES=50
ITERATIONS=300
FUNC=2
MPI_OUT=output/partdiff-par-4.txt
SEQ_OUT=output/partdiff-seq-4.txt
echo "\nTest 4: Prozesse = $NUM_PROC ; Interlines = $INTERLINES; Iterations = $ITERATIONS"
echo "Command: mpirun -np $NUM_PROC ./partdiff-par $SINGLE_THREAD $GAUSS_SEIDEL $INTERLINES $FUNC $NACH_ITERATION $ITERATIONS > $MPI_OUT"
mpirun -np $NUM_PROC ./partdiff-par $SINGLE_THREAD $GAUSS_SEIDEL $INTERLINES $FUNC $NACH_ITERATION $ITERATIONS > $MPI_OUT
./partdiff-seq $SINGLE_THREAD $GAUSS_SEIDEL $INTERLINES $FUNC $NACH_ITERATION $ITERATIONS > $SEQ_OUT
diff $SEQ_OUT $MPI_OUT

#####################################################################

NUM_PROC=5
INTERLINES=100
ITERATIONS=500
FUNC=2
MPI_OUT=output/partdiff-par-5.txt
SEQ_OUT=output/partdiff-seq-5.txt
echo "\nTest 5: Prozesse = $NUM_PROC ; Interlines = $INTERLINES; Iterations = $ITERATIONS"
echo "Command: mpirun -np $NUM_PROC ./partdiff-par $SINGLE_THREAD $GAUSS_SEIDEL $INTERLINES $FUNC $NACH_ITERATION $ITERATIONS > $MPI_OUT"
mpirun -np $NUM_PROC ./partdiff-par $SINGLE_THREAD $GAUSS_SEIDEL $INTERLINES $FUNC $NACH_ITERATION $ITERATIONS > $MPI_OUT
./partdiff-seq $SINGLE_THREAD $GAUSS_SEIDEL $INTERLINES $FUNC $NACH_ITERATION $ITERATIONS > $SEQ_OUT
diff $SEQ_OUT $MPI_OUT

#####################################################################

NUM_PROC=5
INTERLINES=100
PRECISION=1e-4
FUNC=1
MPI_OUT=output/partdiff-par-6.txt
SEQ_OUT=output/partdiff-seq-6.txt
echo "\nTest 6: Prozesse = $NUM_PROC ; Interlines = $INTERLINES; PRECISION = $PRECISION"
echo "Command: mpirun -np $NUM_PROC ./partdiff-par $SINGLE_THREAD $GAUSS_SEIDEL $INTERLINES $FUNC $NACH_GENAUIGKEIT $PRECISION > $MPI_OUT"
mpirun -np $NUM_PROC ./partdiff-par $SINGLE_THREAD $GAUSS_SEIDEL $INTERLINES $FUNC $NACH_GENAUIGKEIT $PRECISION > $MPI_OUT
./partdiff-seq $SINGLE_THREAD $GAUSS_SEIDEL $INTERLINES $FUNC $NACH_GENAUIGKEIT $PRECISION > $SEQ_OUT
diff $SEQ_OUT $MPI_OUT

#####################################################################

NUM_PROC=10
INTERLINES=20
PRECISION=1e-4
FUNC=2
MPI_OUT=output/partdiff-par-7.txt
SEQ_OUT=output/partdiff-seq-7.txt
echo "\nTest 7: Prozesse = $NUM_PROC ; Interlines = $INTERLINES; PRECISION = $PRECISION"
echo "Command: mpirun -np $NUM_PROC ./partdiff-par $SINGLE_THREAD $GAUSS_SEIDEL $INTERLINES $FUNC $NACH_GENAUIGKEIT $PRECISION > $MPI_OUT"
mpirun -np $NUM_PROC ./partdiff-par $SINGLE_THREAD $GAUSS_SEIDEL $INTERLINES $FUNC $NACH_GENAUIGKEIT $PRECISION > $MPI_OUT
./partdiff-seq $SINGLE_THREAD $GAUSS_SEIDEL $INTERLINES $FUNC $NACH_GENAUIGKEIT $PRECISION > $SEQ_OUT
diff $SEQ_OUT $MPI_OUT


echo "Done"
