/****************************************************************************/
/****************************************************************************/
/**                                                                        **/
/**                TU Muenchen - Institut fuer Informatik                  **/
/**                                                                        **/
/** Copyright: Prof. Dr. Thomas Ludwig                                     **/
/**            Andreas C. Schmidt                                          **/
/**                                                                        **/
/** File:      partdiff-seq.c                                              **/
/**                                                                        **/
/** Purpose:   Partial differential equation solver for Gauss-Seidel and   **/
/**            Jacobi method.                                              **/
/**                                                                        **/
/****************************************************************************/
/****************************************************************************/

/* ************************************************************************ */
/* Include standard header file.                                            */
/* ************************************************************************ */
#define _POSIX_C_SOURCE 200809L

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <inttypes.h>
#include <math.h>
#include <malloc.h>
#include <sys/time.h>

#include <mpi.h>
#include "partdiff-mpi.h"
#include <unistd.h>

struct calculation_arguments
{
    uint64_t  N;              /* number of spaces between lines (lines=N+1)     */
    uint64_t  num_matrices;   /* number of matrices                             */
    double    h;              /* length of a space between two lines            */
    double    ***Matrix;      /* index matrix used for addressing M             */
    double    *M;             /* two matrices with real values                  */
};

struct calculation_results
{
    uint64_t  m;
    uint64_t  stat_iteration; /* number of current iteration                    */
    double    stat_precision; /* actual precision of all slaves in iteration    */
};

/* ************************************************************************ */
/* Global variables                                                         */
/* ************************************************************************ */

/* time measurement variables */
struct timeval start_time;       /* time when program started                      */
struct timeval comp_time;        /* time when calculation completed                */


/* ************************************************************************ */
/* initVariables: Initializes some global variables                         */
/* ************************************************************************ */
static void initVariables (struct calculation_arguments* arguments, struct calculation_results* results, struct options const* options) {
    arguments->N = (options->interlines * 8) + 9 - 1;
    arguments->num_matrices = (options->method == METH_JACOBI) ? 2 : 1;
    arguments->h = 1.0 / arguments->N;

    results->m = 0;
    results->stat_iteration = 0;
    results->stat_precision = 0;
}

/* ************************************************************************ */
/* freeMatrices: frees memory for matrices                                  */
/* ************************************************************************ */
static void freeMatrices (struct calculation_arguments* arguments) {
    uint64_t i;

    for (i = 0; i < arguments->num_matrices; i++)
    {
        free(arguments->Matrix[i]);
    }

    free(arguments->Matrix);
    free(arguments->M);
}

/* ************************************************************************ */
/* allocateMemory ()                                                        */
/* allocates memory and quits if there was a memory allocation problem      */
/* ************************************************************************ */
static void* allocateMemory (size_t size) {
    void *p;

    if ((p = malloc(size)) == NULL)
    {
        printf("Speicherprobleme! (%" PRIu64 " Bytes)\n", size);
        /* exit program */
        exit(1);
    }

    return p;
}

/* ************************************************************************ */
/* allocateMatrices: allocates memory for matrices                          */
/* ************************************************************************ */
static void allocateMatrices (struct calculation_arguments* arguments, int length) {
    int i, j;
    int Num_matrices = (int) arguments->num_matrices;
    uint64_t const N = arguments->N;

    // Allokierung der kompletten Matrix
    arguments->M = (double*) allocateMemory(arguments->num_matrices * (N + 1) * (length) * sizeof(double));
    // Allokierung der Pointer zu den jeweiligen Matritzen
    arguments->Matrix = (double***) allocateMemory(arguments->num_matrices * sizeof(double**));

    for (i = 0; i < Num_matrices; i++)
    {
        // Allokierung der Pointer zu den jeweiligen Spalten
        arguments->Matrix[i] = (double**) allocateMemory((length) * sizeof(double*));

        for (j = 0; j < length; j++)
        {
            // Setzen der Pointer auf die jeweiligen Spalten
            arguments->Matrix[i][j] = arguments->M + (i * (N + 1) * (length)) + (j * (N + 1));
        }
    }
}

static void getPart(int* start, int* end, int* length, int rank, int nprocs, int N){
    N -= 1;
    // Ermitteln der Grundlaenge fuer jeden Prozess
    int base = N/nprocs;
    // Ermitteln des Restes 
    int rest = N%nprocs;
    // Aufsplitten des Restes unter den ersten Prozessen ( von 0 aufwaerts)
    if (rank < rest){
        *length = base+1;
    }
    else {
        *length = base;
    }
    // Hinzufuegen der Bufferlines fuer jeden Prozess. Die Bufferlines sind Fuer Rang 0 und nprocs - 1 equivalent zu den Randzeilen
    *length += 2;
    // Ermitteln der vorherigen aufgeteilten Reste, um die absolute Startzeile zu berechnen
    int restaddition;
    if (rank < rest) {
        restaddition = rank+1;
    }
    else {
        restaddition = rest;
    }
    // Berechnen von Start- und Endzeile
    *start = rank*base + restaddition + rank*2;
    *end = *start+ *length-1;
}


/* ************************************************************************ */
/* initMatrices: Initialize matrix/matrices and some global variables       */
/* ************************************************************************ */
static void initMatrices (struct calculation_arguments* arguments, struct options const* options, int rank, int length, int nprocs, int start) {
    int g, i, j;

    int const N = (int) arguments->N;
    double const h = arguments->h;
    double*** Matrix = arguments->Matrix;
    int const Num_Matrices = (int) arguments->num_matrices;

    /* initialize matrix/matrices with zeros */
    for (g = 0; g < Num_Matrices; g++)
    {
        for (i = 0; i < length; i++)
        {
            for (j = 0; j <= N; j++)
            {
                Matrix[g][i][j] = 0.0;
            }
        }
    }

    /* initialize borders, depending on function (function 2: nothing to do) */
    if (options->inf_func == FUNC_F0)
    {
        for (g = 0; g < Num_Matrices; g++)
        {
            // Initialisieren der Seiten
            for (i = 1; i < length-1; i++)
            {
                Matrix[g][i][0] = 1.0 - (h * (i+start));
                Matrix[g][i][N] = h * i;
            }

            // Initialisieren der obersten Reihe fuer rank 0
            if (rank == 0) {
                for (i = 0; i <= N; i++) {
                    Matrix[g][0][i] = 1.0 - (h * (i+start));
                }
            }

            // Initialisieren der untersten Reihe fuer rank nprocs -1
            if (rank == nprocs-1){
                for (i = 0; i <= N; i++) {
                    Matrix[g][length-1][i] = h * (i+start);
                }
            }

            // Initialisieren der Ecken auf 0
            if (rank == nprocs-1){
                Matrix[g][length-1][0] = 0.0;
            }
            if (rank == 0){
                Matrix[g][0][N] = 0.0;
            }
        }
    }
}

/* ************************************************************************ */
/* calculate: solves the equation                                           */
/* ************************************************************************ */
// Das hier ist die alte Funktion, welche unveraendert von Rank 0 fuer Gauss seidel aufgerufen wird.
static void calculate (struct calculation_arguments const* arguments, struct calculation_results *results, struct options const* options) {
    int i, j;                                   /* local variables for loops  */
    int m1, m2;                                 /* used as indices for old and new matrices       */
    double star;                                /* four times center value minus 4 neigh.b values */
    double residuum;                            /* residuum of current iteration                  */
    double maxresiduum;                         /* maximum residuum value of a slave in iteration */

    int const N = arguments->N;
    double const h = arguments->h;

    double pih = 0.0;
    double fpisin = 0.0;

    int term_iteration = options->term_iteration;

    /* initialize m1 and m2 depending on algorithm */
    if (options->method == METH_JACOBI)
    {
        m1 = 0;
        m2 = 1;
    }
    else
    {
        m1 = 0;
        m2 = 0;
    }

    if (options->inf_func == FUNC_FPISIN)
    {
        pih = PI * h;
        fpisin = 0.25 * TWO_PI_SQUARE * h * h;
    }

    while (term_iteration > 0)
    {
        double** Matrix_Out = arguments->Matrix[m1];
        double** Matrix_In  = arguments->Matrix[m2];

        maxresiduum = 0;

        /* over all rows */
        for (i = 1; i < N; i++)
        {
            double fpisin_i = 0.0;

            if (options->inf_func == FUNC_FPISIN)
            {
                fpisin_i = fpisin * sin(pih * (double)i);
            }

            /* over all columns */
            for (j = 1; j < N; j++)
            {
                star = 0.25 * (Matrix_In[i-1][j] + Matrix_In[i][j-1] + Matrix_In[i][j+1] + Matrix_In[i+1][j]);

                if (options->inf_func == FUNC_FPISIN)
                {
                    star += fpisin_i * sin(pih * (double)j);
                }

                if (options->termination == TERM_PREC || term_iteration == 1)
                {
                    residuum = Matrix_In[i][j] - star;
                    residuum = (residuum < 0) ? -residuum : residuum;
                    maxresiduum = (residuum < maxresiduum) ? maxresiduum : residuum;
                }

                Matrix_Out[i][j] = star;
            }
        }

        results->stat_iteration++;
        results->stat_precision = maxresiduum;

        /* exchange m1 and m2 */
        i = m1;
        m1 = m2;
        m2 = i;

        /* check for stopping calculation, depending on termination method */
        if (options->termination == TERM_PREC)
        {
            if (maxresiduum < options->term_precision)
            {
                term_iteration = 0;
            }
        }
        else if (options->termination == TERM_ITER)
        {
            term_iteration--;
        }
    }

    results->m = m2;
}

// Veraenderte Calculate fuer MPI Jacobi
static void calculate_mpi (struct calculation_arguments const* arguments, struct calculation_results *results, struct options const* options, int rank, int length, int nprocs) {
    int i, j;                                   /* local variables for loops  */
    int m1, m2;                                 /* used as indices for old and new matrices       */
    int maxindex = length -1;
    double star;                                /* four times center value minus 4 neigh.b values */
    double residuum;                            /* residuum of current iteration                  */
    double maxresiduum;                         /* maximum residuum value of a slave in iteration */

    int const N = arguments->N;
    double const h = arguments->h;

    double pih = 0.0;
    double fpisin = 0.0;

    int term_iteration = options->term_iteration;

    /* initialize m1 and m2 depending on algorithm */
    if (options->method == METH_JACOBI)
    {
        m1 = 0;
        m2 = 1;
    }
    else
    {
        m1 = 0;
        m2 = 0;
    }

    if (options->inf_func == FUNC_FPISIN)
    {
        pih = PI * h;
        fpisin = 0.25 * TWO_PI_SQUARE * h * h;
    }

    while (term_iteration > 0)
    {
        double** Matrix_Out = arguments->Matrix[m1];
        double** Matrix_In  = arguments->Matrix[m2];

        maxresiduum = 0;
        
        MPI_Request sdown, sup, rdown, rup;

        // Senden der hintersten Reihe an den naechsten Rang
        if (rank < nprocs-1) {
            MPI_Isend(Matrix_In[maxindex-1], N+1, MPI_DOUBLE, rank+1, 0, MPI_COMM_WORLD, &sdown);
        }
        // Senden der vordersten Reihe an den vorherigen Rang
        if (rank > 0) {
            MPI_Isend(Matrix_In[1], N+1, MPI_DOUBLE, rank-1, 0, MPI_COMM_WORLD, &sup);
        }

        // Empfangen der letzten Reihe des vorherigen Prozesses
        if (rank > 0) {
            MPI_Irecv(Matrix_In[0], N+1, MPI_DOUBLE, rank-1, 0, MPI_COMM_WORLD, &rdown);
        }
        // Empfangen der ersten Reihe des naechsten Prozesses
        if (rank < nprocs-1) {
            MPI_Irecv(Matrix_In[maxindex], N+1, MPI_DOUBLE, rank+1, 0, MPI_COMM_WORLD, &rup);
        }

        // Warten
        if (rank < nprocs-1) {
            MPI_Wait(&sdown, MPI_STATUS_IGNORE);
            MPI_Wait(&rup, MPI_STATUS_IGNORE);
        }

        if (rank > 0) {
            MPI_Wait(&sup, MPI_STATUS_IGNORE);
            MPI_Wait(&rdown, MPI_STATUS_IGNORE);
        }


        MPI_Barrier(MPI_COMM_WORLD);


        /* over all rows */
        for (i = 1; i < maxindex; i++)
        {
            double fpisin_i = 0.0;

            if (options->inf_func == FUNC_FPISIN)
            {
                fpisin_i = fpisin * sin(pih * (double)i);
            }

            /* over all columns */
            for (j = 1; j < N; j++)
            {
                star = 0.25 * (Matrix_In[i-1][j] + Matrix_In[i][j-1] + Matrix_In[i][j+1] + Matrix_In[i+1][j]);

                if (options->inf_func == FUNC_FPISIN)
                {
                    star += fpisin_i * sin(pih * (double)j);
                }

                if (options->termination == TERM_PREC || term_iteration == 1)
                {
                    residuum = Matrix_In[i][j] - star;
                    residuum = (residuum < 0) ? -residuum : residuum;
                    maxresiduum = (residuum < maxresiduum) ? maxresiduum : residuum;
                }

                Matrix_Out[i][j] = star;
            }
        }

        results->stat_iteration++;
        results->stat_precision = maxresiduum;

        /* exchange m1 and m2 */
        i = m1;
        m1 = m2;
        m2 = i;

        int residuums[nprocs];
        // Empfangen und Senden aller Residuum
        MPI_Allgather(&residuum, 1, MPI_DOUBLE, &residuums, 1, MPI_DOUBLE, MPI_COMM_WORLD);
        // Ermitteln des global hoechsten Residuums
        for (i=0; i<nprocs; i++) {
            maxresiduum = (residuums[i] < maxresiduum) ? maxresiduum : residuums[i];
        }

        /* check for stopping calculation, depending on termination method */
        if (options->termination == TERM_PREC)
        {
            if (maxresiduum < options->term_precision)
            {
                term_iteration = 0;
            }
        }
        else if (options->termination == TERM_ITER)
        {
            term_iteration--;
        }
    }

    results->m = m2;
}

/* ************************************************************************ */
/*  displayStatistics: displays some statistics about the calculation       */
/* ************************************************************************ */
static void displayStatistics (struct calculation_arguments const* arguments, struct calculation_results const* results, struct options const* options) {
    int N = arguments->N;
    double time = (comp_time.tv_sec - start_time.tv_sec) + (comp_time.tv_usec - start_time.tv_usec) * 1e-6;

    printf("Berechnungszeit:    %f s \n", time);
    printf("Speicherbedarf:     %f MiB\n", (N + 1) * (N + 1) * sizeof(double) * arguments->num_matrices / 1024.0 / 1024.0);
    printf("Berechnungsmethode: ");

    if (options->method == METH_GAUSS_SEIDEL)
    {
        printf("Gauss-Seidel");
    }
    else if (options->method == METH_JACOBI)
    {
        printf("Jacobi");
    }

    printf("\n");
    printf("Interlines:         %" PRIu64 "\n",options->interlines);
    printf("Stoerfunktion:      ");

    if (options->inf_func == FUNC_F0)
    {
        printf("f(x,y) = 0");
    }
    else if (options->inf_func == FUNC_FPISIN)
    {
        printf("f(x,y) = 2pi^2*sin(pi*x)sin(pi*y)");
    }

    printf("\n");
    printf("Terminierung:       ");

    if (options->termination == TERM_PREC)
    {
        printf("Hinreichende Genaugkeit");
    }
    else if (options->termination == TERM_ITER)
    {
        printf("Anzahl der Iterationen");
    }

    printf("\n");
    printf("Anzahl Iterationen: %" PRIu64 "\n", results->stat_iteration);
    printf("Norm des Fehlers:   %e\n", results->stat_precision);
    printf("\n");
}

/****************************************************************************/
/** Beschreibung der Funktion DisplayMatrix:                               **/
/**                                                                        **/
/**
 * rank and size are the MPI rank and size, respectively.
 * from and to denote the global(!) range of lines that this process is responsible for.
 *
 * Example with 9 matrix lines and 4 processes:
 * - rank 0 is responsible for 1-2, rank 1 for 3-4, rank 2 for 5-6 and rank 3 for 7.
 *   Lines 0 and 8 are not included because they are not calculated.
 * - Each process stores two halo lines in its matrix (except for ranks 0 and 3 that only store one).
 * - For instance: Rank 2 has four lines 0-3 but only calculates 1-2 because 0 and 3 are halo lines for other processes. It is responsible for (global) lines 5-6.
 */
/****************************************************************************/
static void DisplayMatrix_mpi (struct calculation_arguments* arguments, struct calculation_results* results, struct options* options, int rank, int size, int from, int to) {
    int const elements = 8 * options->interlines + 9;

    int x, y;
    double** Matrix = arguments->Matrix[results->m];
    MPI_Status status;

    /* first line belongs to rank 0 */
    if (rank == 0)
        from--;

    /* last line belongs to rank size - 1 */
    if (rank + 1 == size)
        to++;

    if (rank == 0)
        printf("Matrix:\n");

    for (y = 0; y < 9; y++)
    {
        int line = y * (options->interlines + 1);

        if (rank == 0)
        {
            /* check whether this line belongs to rank 0 */
            if (line < from || line > to)
            {
                /* use the tag to receive the lines in the correct order
                 * the line is stored in Matrix[0], because we do not need it anymore */
                MPI_Recv(Matrix[0], elements, MPI_DOUBLE, MPI_ANY_SOURCE, 42 + y, MPI_COMM_WORLD, &status);
            }
        }
        else
        {
            if (line >= from && line <= to)
            {
                /* if the line belongs to this process, send it to rank 0
                 * (line - from + 1) is used to calculate the correct local address */
                MPI_Send(Matrix[line - from + 1], elements, MPI_DOUBLE, 0, 42 + y, MPI_COMM_WORLD);
            }
        }

        if (rank == 0)
        {
            for (x = 0; x < 9; x++)
            {
                int col = x * (options->interlines + 1);

                if (line >= from && line <= to)
                {
                    /* this line belongs to rank 0 */
                    printf("%7.4f", Matrix[line][col]);
                }
                else
                {
                    /* this line belongs to another rank and was received above */
                    printf("%7.4f", Matrix[0][col]);
                }
            }

            printf("\n");
        }
    }

    fflush(stdout);
}

// Old DisplayMatrix
static void DisplayMatrix (struct calculation_arguments* arguments, struct calculation_results* results, struct options* options)
{
int x, y;

double** Matrix = arguments->Matrix[results->m];

int const interlines = options->interlines;

printf("Matrix:\n");

for (y = 0; y < 9; y++)
{
    for (x = 0; x < 9; x++)
    {
        printf ("%7.4f", Matrix[y * (interlines + 1)][x * (interlines + 1)]);
    }

    printf ("\n");
}

fflush (stdout);
}

/* ************************************************************************ */
/*  main                                                                    */
/* ************************************************************************ */
int main (int argc, char** argv) {
    struct options options;
    struct calculation_arguments arguments;
    struct calculation_results results;

    int rank, nprocs;
    /* get parameters */

    // Init functions for MPI
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

    // Get Parameter
    AskParams(&options, argc, argv, rank);
    // Initialize Variables
    initVariables(&arguments, &results, &options);

    // Falls Gauss_Seidel vorliegt, soll Rank 0 so rechnen, als waere er alleine vorhanden
    if (options.method == METH_GAUSS_SEIDEL && rank == 0) {
        nprocs = 1;
    }

    int start, end, length;
    // Get Segment of matrix per rank
    getPart(&start, &end, &length, rank, nprocs, arguments.N);
    // Allocate stuff
    allocateMatrices(&arguments, length);
    initMatrices(&arguments, &options, rank, length, nprocs, start);

    if (options.method == METH_GAUSS_SEIDEL) {
        // Nur Rank 0 soll Berechnungen fuer Gauss Seidel machen
        if (rank == 0 ){
            gettimeofday(&start_time, NULL);
            calculate(&arguments, &results, &options);
            gettimeofday(&comp_time, NULL);

            displayStatistics(&arguments, &results, &options);
            DisplayMatrix(&arguments, &results, &options);
        }
    } else {
        gettimeofday(&start_time, NULL);
        calculate_mpi(&arguments, &results, &options, rank, length, nprocs);
        gettimeofday(&comp_time, NULL);

        if (rank == 0 ) {
            displayStatistics(&arguments, &results, &options);
        }
        DisplayMatrix_mpi(&arguments, &results, &options, rank, length, start, end);
    }
    freeMatrices(&arguments);

    MPI_Finalize();
    return EXIT_SUCCESS;
}