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
#include <sys/time.h>
#include <mpi.h>
#include <stdbool.h>
#include <float.h>
#include <limits.h>

#include "split.h"
#include "partdiff-par.h"

#define MIN(a,b) ((a) < (b) ? (a) : (b))
#define MAX(a,b) ((a) > (b) ? (a) : (b))

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


static int rank;
static int comm_size;

static int start_line;
static int end_line; // inclusive
static int num_lines;

static bool is_first;
static bool is_last;

/* ************************************************************************ */
/* initVariables: Initializes some global variables                         */
/* ************************************************************************ */
static void initVariables (struct calculation_arguments* arguments, struct calculation_results* results, struct options const* options)
{
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
static void freeMatrices (struct calculation_arguments* arguments)
{
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
static void* allocateMemory (size_t size)
{
    void *p;

    if ((p = malloc(size)) == NULL)
    {
        printf("Speicherprobleme! (%lld Bytes)\n", (unsigned long long)size);
        /* exit program */
        exit(1);
    }

    return p;
}

/* ************************************************************************ */
/* allocateMatrices: allocates memory for matrices                          */
/* ************************************************************************ */
static void allocateMatrices (struct calculation_arguments* arguments)
{
    uint64_t i, j;

    uint64_t const N = arguments->N;

    arguments->M = allocateMemory(arguments->num_matrices * (N + 1) * (num_lines + 2) * sizeof(double));
    arguments->Matrix = allocateMemory(arguments->num_matrices * sizeof(double**));

    for (i = 0; i < arguments->num_matrices; i++)
    {
        arguments->Matrix[i] = allocateMemory((num_lines + 2) * sizeof(double*));

        for (j = 0; j < (uint64_t)num_lines + 2; j++)
        { //                                            | size of one matrix   |    |size of line|
            arguments->Matrix[i][j] = arguments->M + (i * (num_lines + 2) * (N + 1)) + (j * (N + 1));
        }
    }
}

static double *ml(double **m, int l) {
    return m[1 + l - start_line];
}

/* ************************************************************************ */
/* initMatrices: Initialize matrix/matrices and some global variables       */
/* ************************************************************************ */
static void initMatrices (struct calculation_arguments* arguments, struct options const* options)
{
    int64_t g, i, j;                                /*  local variables for loops   */

    int64_t const N = arguments->N;
    double const h = arguments->h;
    double*** Matrix = arguments->Matrix;

    /* initialize matrix/matrices with zeros */
    for (g = 0; g < (int64_t)arguments->num_matrices; g++)
    {
        for (i = 0; i < (int64_t)num_lines + 2; i++)
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
        for (g = 0; g < (int64_t)arguments->num_matrices; g++)
        {
            double **m = Matrix[g];
            for (i = 0; i <= N; i++)
            {
                if (is_first) ml(m, 0)[i] = 1.0 - (h * i);
                if (is_last)  ml(m, N)[i] = h * i;

                if (i >= start_line && i <= end_line) {
                    ml(m, i)[0] = 1.0 - (h * i);
                    ml(m, i)[N] = h * i;
                }
            }

            if (is_first) ml(m, 0)[N] = 0.0;
            if (is_last)  ml(m, N)[0] = 0.0;
        }
    }
}

/* ************************************************************************ */
/* Gauss Seidel Implementation der Calculate Funktion                       */
/* ************************************************************************ */
static void calculate_gauss_seidel (struct calculation_arguments const* arguments, struct calculation_results *results, struct options const* options, int nprocs)
{
    int i, j;/* local variables for loops  */
    int m1, m2;/* used as indices for old and new matrices       */
    double star;/* four times center value minus 4 neigh.b values */
    double residuum;/* residuum of current iteration                  */
    double maxresiduum;/* maximum residuum value of a slave in iteration */
    int iteration, stopping; // Current iteration count

    int const N = arguments->N;
    double const h = arguments->h;

    double pih = 0.0;
    double fpisin = 0.0;

    int term_iteration = options->term_iteration;
    MPI_Request drop1, drop2;

    m1 = 0;
    m2 = 0;
    maxresiduum = DBL_MAX;

    if (options->inf_func == FUNC_FPISIN)
    {
        pih = PI * h;
        fpisin = 0.25 * TWO_PI_SQUARE * h * h;
    }

    bool is_first_iteration = true;

    iteration = 0;
    stopping = -1;

    // Prozesse die noch nicht arbeiten können, werden beschäftigt.
    for (int k = 0; k < rank; k++)
    {
        MPI_Allreduce(MPI_IN_PLACE, &maxresiduum, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
        MPI_Allreduce(MPI_IN_PLACE, &stopping, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
    }

    while (term_iteration > 0 && (stopping == (-1) || iteration <= stopping))
    {
        //printf("%d: starte iteration\n", rank);
        /* Es wird nur eine Matrix benoetigt */
        double** mo = arguments->Matrix[m1];

        /* Receiving lines from previous and next process */
        if(rank != 0){
            MPI_Recv(ml(mo, start_line - 1), N+1, MPI_DOUBLE, rank-1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }
        if(rank != nprocs-1 && is_first_iteration != true){
            MPI_Recv(ml(mo, end_line + 1), N+1, MPI_DOUBLE, rank+1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }

        maxresiduum = 0;

        // Die entsprechenden Reduces zu denen for der while-Schleife
        MPI_Allreduce(MPI_IN_PLACE, &maxresiduum, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
        MPI_Allreduce(MPI_IN_PLACE, &stopping, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);

        // Abruch nach Genauigkeit. Hier geht was schief, aber wir wissen nicht, wieso.
        if (stopping != (-1) && iteration > stopping)
            break;

        is_first_iteration = false;
        maxresiduum = 0;

        // Das Herz der Berechnung. Hier gibt's nicht viel zu sehen.
        /* over all rows */
        for (i = MAX(start_line, 1); i <= MIN(end_line, N-1); i++)
        {
            double fpisin_i = 0.0;

            if (options->inf_func == FUNC_FPISIN)
            {
                fpisin_i = fpisin * sin(pih * (double)i);
            }

            /* over all columns */
            for (j = 1; j < N; j++)
            {
                double x1 = ml(mo, i-1)[j];
                double x2 = ml(mo, i)[j-1];
                double x3 = ml(mo, i)[j+1];
                double x4 = ml(mo, i+1)[j];
                star = 0.25 * (x1 + x2 + x3 + x4);

                if (options->inf_func == FUNC_FPISIN)
                {
                    star += fpisin_i * sin(pih * (double)j);
                }

                if (options->termination == TERM_PREC || term_iteration == 1)
                {
                    residuum = ml(mo, i)[j] - star;
                    residuum = (residuum < 0) ? -residuum : residuum;
                    maxresiduum = (residuum < maxresiduum) ? maxresiduum : residuum;
                }

                ml(mo, i)[j] = star;
            }
        }

        /* Sending lines to previous and next process */
        if (rank != 0) {
            MPI_Isend(ml(mo, start_line), N + 1, MPI_DOUBLE, rank - 1, 0, MPI_COMM_WORLD, &drop1);
        }
        if (rank != nprocs-1) {
            MPI_Isend(ml(mo, end_line), N + 1, MPI_DOUBLE, rank + 1, 0, MPI_COMM_WORLD, &drop2);
        }

        // Noch mehr reduces.
        MPI_Allreduce(MPI_IN_PLACE, &maxresiduum, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
        MPI_Allreduce(MPI_IN_PLACE, &stopping, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);

        //results->stat_iteration++;
        //results->stat_precision = maxresiduum;
        //iteration++;

        // Abbruchbedingungen. Abbruch nach Iterationsanzahl klappt wunderbar,
        // Abbruch nach Genauigkeit nicht. :-(
        if (options->termination == TERM_PREC)
        {

            if (maxresiduum < options->term_precision && stopping < 0)
            {
                stopping = rank + iteration;
            }
        }

        if (options->termination == TERM_ITER)
        {
                term_iteration--;
        }

        results->stat_iteration++;
        iteration++;
        //printf("%d: Iteration %d von %d, %e\n", rank, iteration, stopping, maxresiduum);
        
    }

    // Ähnlich wie for der while-Schleife müssen die fertigen Prozesse
    // beschäftigt werden.
    for (int k = 0; k < comm_size - 1 - rank; k++)
    {
        MPI_Allreduce(MPI_IN_PLACE, &maxresiduum, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
        MPI_Allreduce(MPI_IN_PLACE, &stopping, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
    }

    // Hier sammeln wir nochmal das maxresiduum ein. Um ehrlich zu sein sind wir uns
    // gar nicht mehr sicher, ob wir diese Zeile brauchen. Sie macht aber nichts kaputt.
    MPI_Allreduce(MPI_IN_PLACE, &maxresiduum, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);

    results->stat_precision = maxresiduum;

    results->m = m2;

}

/* ************************************************************************ */
/* calculate: solves the equation                                           */
/* ************************************************************************ */
static void calculate_jacobi (struct calculation_arguments const* arguments, struct calculation_results *results, struct options const* options)
{
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

    m1 = 0;
    m2 = 1;

    if (options->inf_func == FUNC_FPISIN)
    {
        pih = PI * h;
        fpisin = 0.25 * TWO_PI_SQUARE * h * h;
    }

    MPI_Request send_to_prev_req;
    MPI_Request send_to_next_req;

    bool is_first_iteration = true;

    while (term_iteration > 0)
    {
        //printf("%d: wheeeeeeeee\n", rank);

        double** mo = arguments->Matrix[m1];
        double** mi  = arguments->Matrix[m2];

        if (!is_first) {
            // wait for previous send to finish
            // if first iteration: don't wait
            if (!is_first_iteration) {
                MPI_Wait(&send_to_prev_req, MPI_STATUS_IGNORE);
            }

            // send first line to previous process
            MPI_Isend(ml(mi, start_line), N + 1, MPI_DOUBLE, rank - 1, 0, MPI_COMM_WORLD, &send_to_prev_req);
        }
        if (!is_last) {
            // wait for previous send to finish
            // if first iteration: don't wait
            if (!is_first_iteration) {
                MPI_Wait(&send_to_next_req, MPI_STATUS_IGNORE);
            }

            // send last line to next process
            MPI_Isend(ml(mi, end_line), N + 1, MPI_DOUBLE, rank + 1, 0, MPI_COMM_WORLD, &send_to_next_req);
        }

        if (!is_first) {
            // get the previous line from predecessor process
            MPI_Recv(ml(mi, start_line - 1), N + 1, MPI_DOUBLE, rank - 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }
        if (!is_last) {
            // get successive line from successor process
            MPI_Recv(ml(mi, end_line + 1), N + 1, MPI_DOUBLE, rank + 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }

        is_first_iteration = false;
        maxresiduum = 0;

        /* over all rows */
        for (i = MAX(start_line, 1); i <= MIN(end_line, N-1); i++)
        {
            double fpisin_i = 0.0;

            if (options->inf_func == FUNC_FPISIN)
            {
                fpisin_i = fpisin * sin(pih * (double)i);
            }

            /* over all columns */
            for (j = 1; j < N; j++)
            {
                double x1 = ml(mi, i-1)[j];
                double x2 = ml(mi, i)[j-1];
                double x3 = ml(mi, i)[j+1];
                double x4 = ml(mi, i+1)[j];
                star = 0.25 * (x1 + x2 + x3 + x4);

                if (options->inf_func == FUNC_FPISIN)
                {
                    star += fpisin_i * sin(pih * (double)j);
                }

                if (options->termination == TERM_PREC || term_iteration == 1)
                {
                    residuum = ml(mi, i)[j] - star;
                    residuum = (residuum < 0) ? -residuum : residuum;
                    maxresiduum = (residuum < maxresiduum) ? maxresiduum : residuum;
                }

                ml(mo, i)[j] = star;
            }
        }
        // Reduce maxresiduum over all processes to ensure correct termination by precision
        MPI_Allreduce(MPI_IN_PLACE, &maxresiduum, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);

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
    results->stat_precision = maxresiduum;

    results->m = m2;
}

/* ************************************************************************ */
/*  displayStatistics: displays some statistics about the calculation       */
/* ************************************************************************ */
static void displayStatistics (struct calculation_arguments const* arguments, struct calculation_results const* results, struct options const* options)
{
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
/** Die Funktion DisplayMatrix gibt eine Matrix                            **/
/** in einer "ubersichtlichen Art und Weise auf die Standardausgabe aus.   **/
/**                                                                        **/
/** Die "Ubersichtlichkeit wird erreicht, indem nur ein Teil der Matrix    **/
/** ausgegeben wird. Aus der Matrix werden die Randzeilen/-spalten sowie   **/
/** sieben Zwischenzeilen ausgegeben.                                      **/
/****************************************************************************/
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
static void DisplayMatrix (struct calculation_arguments* arguments, struct calculation_results* results, struct options* options, int rank, int size, int from, int to)
{
    int const elements = 8 * options->interlines + 9;

    int x, y;
    double** Matrix = arguments->Matrix[results->m];
    MPI_Status status;

    /* first line belongs to rank 0 */
    if (rank == 0){
        //from--;
    }

    /* last line belongs to rank size - 1 */
    if (rank + 1 == size) {
        //to++;
    }

    if (rank == 0) {
        printf("Matrix:\n");
    }

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
                MPI_Send(ml(Matrix, line), elements, MPI_DOUBLE, 0, 42 + y, MPI_COMM_WORLD);
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
                    printf("%7.4f", ml(Matrix, line)[col]);
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


/* ************************************************************************ */
/*  main                                                                    */
/* ************************************************************************ */
int main (int argc, char** argv)
{
    MPI_Init(&argc, &argv);

    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &comm_size);


    is_first = rank == 0;
    is_last = rank == comm_size - 1;

    struct options options;
    struct calculation_arguments arguments;
    struct calculation_results results;

    if (is_first) PrintHeader();

    /* get parameters */
    AskParams(&options, argc, argv);

    initVariables(&arguments, &results, &options);


    int **split_res = split(arguments.N + 1, comm_size);
    int *my_split = split_res[rank];
    start_line = my_split[0];
    end_line = my_split[1];
    num_lines = end_line + 1 - start_line;

    if (start_line != -1) {
        allocateMatrices(&arguments);
        initMatrices(&arguments, &options);

        gettimeofday(&start_time, NULL);
        if (options.method == METH_JACOBI) {
            calculate_jacobi(&arguments, &results, &options);
        }
        else {
            calculate_gauss_seidel(&arguments, &results, &options, comm_size);
        }
        gettimeofday(&comp_time, NULL);

        if(is_first) displayStatistics(&arguments, &results, &options);
        DisplayMatrix(&arguments, &results, &options, rank, comm_size, start_line, end_line);

        freeMatrices(&arguments);
    }

    MPI_Finalize();

    return 0;
}