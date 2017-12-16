/****************************************************************************/
/****************************************************************************/
/**                                                                        **/
/**                 TU München - Institut für Informatik                   **/
/**                                                                        **/
/** Copyright: Prof. Dr. Thomas Ludwig                                     **/
/**            Andreas C. Schmidt                                          **/
/**                                                                        **/
/** File:      partdiff-seq.c                                              **/
/**                                                                        **/
/** Purpose:   Partial differential equation solver for Gauß-Seidel and   **/
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
#include "partdiff-par.h"


int g_rank;					/* process rank */
int g_num_procs;		/* number of processes working */
int g_minMat;				/* lower index for matrix-section */
int g_maxMat;				/* upper index for matrix-section */
int g_size;				 /* number of matrix rows */
uint64_t g_alloc_size;

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
static
void
initVariables (struct calculation_arguments* arguments, struct calculation_results* results, struct options const* options)
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
static
void
freeMatrices (struct calculation_arguments* arguments)
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
static
void*
allocateMemory (size_t size)
{
	void *p;

	if ((p = malloc(size)) == NULL)
	{
		printf("Speicherprobleme! (%" PRIu64 " Bytes angefordert)\n", size);
		exit(1);
	}

	return p;
}


/* ************************************************************************ */
/* setLowAndHigh: sets Bounds of Sub-Matrix for the calling process         */
/* ************************************************************************ */
static
void
setLowAndHigh(int num_rows) {

	// Split all rows equally [base_size]
	int base_size = num_rows / g_num_procs;
	// and give rest to first n <= rest
	int num_big = num_rows % g_num_procs;
	/*
     * number of processes that receive one row more = rest => (0 --> rest - 1) => rank < num_big
     * calculate number of processes that receive one more, prior to current process
     * eg: num_big = 4 & rank = 3 ==> rank 0, 1, 2 have one more also --> #rank prior to this one
     * eg: num_big = 4 & rank = 4 ==> rank 0, 1, 2, 3 have one more --> num_big prior to this one
     */
	int num_big_pres = (g_rank < num_big) ? g_rank : num_big;
	/*
     * number of processes that do not receive one more, prior to current process
     * max(rank - num_big_pres, 0)
     * eg: num_big = 4 & rank = 3 ==> 3 bigger prior => 3 - 3 = 0; no small prior
     * eg: num_big = 4 & rank = 5 ==> 4 bigger prior (0,1,2,3) => 5 - 4 = 1; 1 smaller prior (4)
     */
	int num_sml_pres = (g_rank - num_big_pres > 0) ? g_rank - num_big_pres : 0;
	// number of bigger prior to this * bigger-size + number of smaller prior * base-size
	g_minMat = (num_big_pres * (base_size + 1)) + num_sml_pres * base_size;
	// if bigger: minIndex + bigSie; else: minIndex + smallSize
	g_size   = (g_rank < num_big) ? base_size + 1 : base_size;
	g_maxMat = g_minMat + g_size - 1;

	// extra rows for data from other processes
	// for the first and last rank only one extra row
	if (g_rank == 0 || g_rank == g_num_procs - 1) {
		g_alloc_size = g_size + 1;
	} else {
		g_alloc_size = g_size + 2;
	}
	printf("[Rank = %d] %d rows -> %d processes. Rank %d: %d -> %d (%d rows)\n", g_rank, num_rows, g_num_procs, g_rank, g_minMat, g_maxMat, g_size);
}


/* ************************************************************************ */
/* allocateMatrices: allocates memory for matrices                          */
/* ************************************************************************ */
static
void
allocateMatrices (struct calculation_arguments* arguments)
{
	uint64_t i, j;

	uint64_t const N = arguments->N;

	// set the appropriate g_minMat and g_maxMat
	setLowAndHigh(N);

	printf("[Rank = %d] Allocating %ld rows of memory.\n", g_rank, g_alloc_size);
	arguments->M = allocateMemory(arguments->num_matrices * g_alloc_size * (N + 1) * sizeof(double));
	arguments->Matrix = allocateMemory(arguments->num_matrices * sizeof(double**));

	for (i = 0; i < arguments->num_matrices; i++)
	{
		arguments->Matrix[i] = allocateMemory(g_alloc_size * sizeof(double*));

		for (j = 0; j < g_alloc_size; j++)
		{
			arguments->Matrix[i][j] = arguments->M + (i * g_alloc_size * (N + 1)) + (j * (N + 1));
		}
	}
}

/* ************************************************************************ */
/* initMatrices: Initialize matrix/matrices and some global variables       */
/* ************************************************************************ */
static
void
initMatrices (struct calculation_arguments* arguments, struct options const* options)
{
	uint64_t g, i, j;                                /*  local variables for loops   */

	uint64_t const N = arguments->N;
	double const h = arguments->h;
	double*** Matrix = arguments->Matrix;

	/* initialize matrix/matrices with zeros */
	for (g = 0; g < arguments->num_matrices; g++)
	{
		for (i = 0; i < g_alloc_size; i++)          // allocated rows + overlap (in [1,2])
		{
			for (j = 0; j <= N; j++)                // columns
			{
				Matrix[g][i][j] = 0.0;
			}
		}
	}

	/* initialize borders, depending on function (function 2: nothing to do) */
	if (options->inf_func == FUNC_F0) {
		for (g = 0; g < arguments->num_matrices; g++) {
			if (g_rank == 0) {
				for (i = 0; i <= N; i++)
				{
					Matrix[g][0][i] = 1.0 - (h * i);
				}
			} // if instead of else if -> can be only one
			if (g_rank == g_num_procs - 1) {
				for (i = 0; i <= N; i++)
				{
					Matrix[g][g_alloc_size-1][i] = h * i;
				}
			}
			// all
			for (i = 0; i < g_alloc_size; i++)
			{
				// emulate i over all rows
				j = g_minMat + i;
				Matrix[g][i][0] = 1.0 - (h * j);
				Matrix[g][i][N] = h * j;
			}

			if (g_rank == 0) {
				// first row of whole matrix
				Matrix[g][0][N] = 0.0;
			} // if instead of else if -> can be only one
			if (g_rank == g_num_procs - 1) {
				// last row of whole matrix
				Matrix[g][g_alloc_size-1][0] = 0.0;
			}
		}
	}
}

/* ************************************************************************ */
/* calculate: solves the equation                                           */
/* ************************************************************************ */
static
void
calculate (struct calculation_arguments const* arguments, struct calculation_results* results, struct options const* options)
{
	int i, j;                                   /* local variables for loops */
	int m1, m2;                                 /* used as indices for old and new matrices */
	double star;                                /* four times center value minus 4 neigh.b values */
	double residuum;                            /* residuum of current iteration */
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

		if (g_num_procs > 1)
		{
			// MSGs between processes
			MPI_Request send_to_next, send_to_prev, get_from_next, get_from_prev;
			// Send first to prev
			if (g_rank > 0){
				// Index 1 because 0 is buffer (dont send)
				MPI_Isend(Matrix_In[1], N+1, MPI_DOUBLE, g_rank-1, 0, MPI_COMM_WORLD, &send_to_prev);
			}
			// Send last to next
			if (g_rank < g_num_procs - 1){
				// Index g_alloc_size - 2 because last is buffer (dont send)
				MPI_Isend(Matrix_In[g_alloc_size-2], N+1, MPI_DOUBLE, g_rank + 1, 0, MPI_COMM_WORLD, &send_to_next);
			}
			// Get first from next
			if (g_rank < g_num_procs - 1){
				// Index g_alloc_size - 1 because last is buffer (receive)
				MPI_Irecv(Matrix_In[g_alloc_size - 1], N+1, MPI_DOUBLE, g_rank-1, 0, MPI_COMM_WORLD, &get_from_next);
			}
			// Get last from prev
			if (g_rank > 0){
				// Index 0 because 0 is buffer (receive)
				MPI_Irecv(Matrix_In[0], N+1, MPI_DOUBLE, g_rank+1, 0, MPI_COMM_WORLD, &get_from_prev);
			}

			// Wait for reqeusts to finish
			if (g_rank > 0){
				MPI_Wait(&send_to_prev, MPI_STATUS_IGNORE);
				MPI_Wait(&get_from_next, MPI_STATUS_IGNORE);
			}
			if (g_rank < g_num_procs - 1){
				MPI_Wait(&send_to_next, MPI_STATUS_IGNORE);
				MPI_Wait(&get_from_prev, MPI_STATUS_IGNORE);
			}

			// Synch
			MPI_Barrier(MPI_COMM_WORLD);
		}

		printf("[Rank = %d] Finish sending and receiving\n", g_rank);

		/* over all rows */
		for (i = 1; i < g_alloc_size-1; i++)
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

		// reduce maxresidium
		if (g_num_procs > 1) {
			MPI_Allreduce(&residuum, &maxresiduum, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
		}


		/* check for stopping calculation depending on termination method */
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
static
void
displayStatistics (struct calculation_arguments const* arguments, struct calculation_results const* results, struct options const* options)
{
	int N = arguments->N;
	double time = (comp_time.tv_sec - start_time.tv_sec) + (comp_time.tv_usec - start_time.tv_usec) * 1e-6;

	printf("Berechnungszeit:    %f s \n", time);
	printf("Speicherbedarf:     %f MiB\n", (N + 1) * (N + 1) * sizeof(double) * arguments->num_matrices / 1024.0 / 1024.0);
	printf("Berechnungsmethode: ");

	if (options->method == METH_GAUSS_SEIDEL)
	{
		printf("Gauß-Seidel");
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
static
void
DisplayMatrix (struct calculation_arguments* arguments, struct calculation_results* results, struct options* options, int rank, int size, int from, int to)
{
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

static
void
DisplayMatrix_single (struct calculation_arguments* arguments, struct calculation_results* results, struct options* options)
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
int
main (int argc, char** argv)
{
	struct options options;
	struct calculation_arguments arguments;
	struct calculation_results results;

	AskParams(&options, argc, argv);

	initVariables(&arguments, &results, &options);

	// MPI INIT
	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &g_rank);
	if (options.method == METH_JACOBI){
		MPI_Comm_size(MPI_COMM_WORLD, &g_num_procs);
	} else {
		// run METH_GAUSS_SEIDEL sequential
		g_num_procs = 1;
		// OR:
		// return -1;
	}

	allocateMatrices(&arguments);
	initMatrices(&arguments, &options);

	gettimeofday(&start_time, NULL);
	calculate(&arguments, &results, &options);
	gettimeofday(&comp_time, NULL);

	if (g_rank == 0) {
		displayStatistics(&arguments, &results, &options);
	}
	MPI_Barrier(MPI_COMM_WORLD);
	if (g_num_procs > 1){
		DisplayMatrix(&arguments, &results, &options, g_rank, g_size, g_minMat, g_maxMat);
	} else {
		DisplayMatrix_single(&arguments, &results, &options);
	}


	freeMatrices(&arguments);
	MPI_Finalize();
	return 0;
}

