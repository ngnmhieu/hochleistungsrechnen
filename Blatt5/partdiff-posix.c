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
#include <pthread.h>

#include "partdiff-seq.h"

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
/* allocateMatrices: allocates memory for matrices                          */
/* ************************************************************************ */
typedef struct {
  int tid;
  int N;
  int i;
  int rows_per_thread;
  double*** Matrix;
  double *M;
} allocate_matrix_data;

void* allocate_submatrix (void* arg)
{
  allocate_matrix_data *data = (allocate_matrix_data*) arg;

  int tid             = data->tid;
  int N               = data->N;
  int i               = data->i;
  int rows_per_thread = data->rows_per_thread;
  double*** Matrix    = data->Matrix;
  double *M           = data->M;

  int start_row = tid * rows_per_thread;
  int end_row = (tid + 1) * rows_per_thread;
  end_row = end_row < (N+1) ? end_row : (N+1);

  /* Debug message */
  /* printf("[TID = %d] Starting allocate matrix from row %d to %d\n", tid, start_row, end_row - 1); */

  for (int j = start_row; j < end_row; j++)
  {
    Matrix[i][j] = M + (i * (N + 1) * (N + 1)) + (j * (N + 1));
  }
  free(data);
  pthread_exit(NULL);
}

static
void
allocateMatrices (struct calculation_arguments* arguments, struct options const* options)
{
	uint64_t i;

	uint64_t const N = arguments->N;

  pthread_t threads[options->number];
  int rows_per_thread = (int) ceil(1.0 * (N + 1) / options->number);

	arguments->M = allocateMemory(arguments->num_matrices * (N + 1) * (N + 1) * sizeof(double));
	arguments->Matrix = allocateMemory(arguments->num_matrices * sizeof(double**));

	for (i = 0; i < arguments->num_matrices; i++)
	{
		arguments->Matrix[i] = allocateMemory((N + 1) * sizeof(double*));

    /* iterate over all threads */
		for (int j = 0; j < (int) options->number; j++)
		{
      allocate_matrix_data *data = (allocate_matrix_data*) malloc(sizeof(allocate_matrix_data));
      data->tid = j;
      data->N = N;
      data->i = i;
      data->rows_per_thread = rows_per_thread;
      data->Matrix = arguments->Matrix;
      data->M = arguments->M;
      pthread_create(&threads[j], NULL, allocate_submatrix, (void*) data);
		}

    /* joining all threads */
		for (int j = 0; j < (int) options->number; j++)
		{
      pthread_join(threads[j], NULL);
		}
	}
}

/* ************************************************************************ */
/* initMatrices: Initialize matrix/matrices and some global variables       */
/* ************************************************************************ */
typedef struct {
  int tid;
  int N;
  int g;
  double h;
  int rows_per_thread;
  double*** Matrix;
} init_matrix_data;

void* init_rows (void* arg)
{
  init_matrix_data *data = (init_matrix_data*) arg;

  int tid             = data->tid;
  int N               = data->N;
  int g               = data->g;
  int rows_per_thread = data->rows_per_thread;
  double*** Matrix    = data->Matrix;

  int start_row = tid * rows_per_thread;
  int end_row = (tid + 1) * rows_per_thread;
  end_row = end_row < (N+1) ? end_row : (N+1);

  /* Debug message */
  /* printf("[TID = %d] Starting init matrix from row %d to %d\n", tid, start_row, end_row - 1); */

  for (int i = start_row; i < end_row; i++)
  {
    for (int j = 0; j <= N; j++)
    {
      Matrix[g][i][j] = 0.0;
    }
  }
  free(data);
  pthread_exit(NULL);
}

void* init_borders (void* arg)
{
  init_matrix_data *data = (init_matrix_data*) arg;

  int tid             = data->tid;
  int N               = data->N;
  int g               = data->g;
  double h            = data->h;
  int rows_per_thread = data->rows_per_thread;
  double*** Matrix    = data->Matrix;

  int start_row = tid * rows_per_thread;
  int end_row = (tid + 1) * rows_per_thread;
  end_row = end_row < (N+1) ? end_row : (N+1);

  /* Debug message */
  /* printf("[TID = %d] Starting init matrix border from row %d to %d\n", tid, start_row, end_row - 1); */

  for (int i = start_row; i < end_row; i++)
  {
    Matrix[g][i][0] = 1.0 - (h * i);
    Matrix[g][i][N] = h * i;
    Matrix[g][0][i] = 1.0 - (h * i);
    Matrix[g][N][i] = h * i;
  }
  free(data);
  pthread_exit(NULL);
}

static
void
initMatrices (struct calculation_arguments* arguments, struct options const* options)
{
	uint64_t g, i;                                /*  local variables for loops   */

	uint64_t const N = arguments->N;
	double const h = arguments->h;
	double*** Matrix = arguments->Matrix;

  pthread_t threads[options->number];
  int rows_per_thread = (int) ceil(1.0 * (N + 1) / options->number);

	/* initialize matrix/matrices with zeros */
	for (g = 0; g < arguments->num_matrices; g++)
	{

    /* Debug message */
    /* printf("Matrix g = %ld ; N+1= %ld; rows_per_thread = %d \n", g, N+1, rows_per_thread); */
  
    /* iterate over all threads */
		for (i = 0; i < options->number; i++)
		{
      init_matrix_data *data = (init_matrix_data*) malloc(sizeof(init_matrix_data));
      data->tid = i;
      data->N = N;
      data->g = g;
      data->h = h;
      data->rows_per_thread = rows_per_thread;
      data->Matrix = Matrix;
      pthread_create(&threads[i], NULL, init_rows, (void*) data);
		}

    /* joining all threads */
		for (i = 0; i < options->number; i++)
		{
      pthread_join(threads[i], NULL);
		}
	}

	/* initialize borders, depending on function (function 2: nothing to do) */
	if (options->inf_func == FUNC_F0)
	{
		for (g = 0; g < arguments->num_matrices; g++)
		{
      /* Debug message */
      /* printf("Matrix g = %ld ; N+1= %ld; rows_per_thread = %d \n", g, N+1, rows_per_thread); */

      /* iterate over all threads */
      for (i = 0; i < options->number; i++)
			{
        init_matrix_data *data = (init_matrix_data*) malloc(sizeof(init_matrix_data));
        data->tid = i;
        data->N = N;
        data->g = g;
        data->h = h;
        data->rows_per_thread = rows_per_thread;
        data->Matrix = Matrix;
        pthread_create(&threads[i], NULL, init_borders, (void*) data);
			}

      /* joining all threads */
      for (i = 0; i < options->number; i++)
      {
        pthread_join(threads[i], NULL);
      }

			Matrix[g][N][0] = 0.0;
			Matrix[g][0][N] = 0.0;
		}
	}
}

/* ************************************************************************ */
/* calculate: solves the equation                                           */
/* ************************************************************************ */
typedef struct {
   int tid;
   double** Matrix_In;
   double** Matrix_Out;
   double fpisin;
   double pih;
   int term_iteration;
   double *maxresiduum;
   struct options const* options;
   int rows_per_thread;
   int N;
   pthread_mutex_t *mutex;
} calculate_row_data;

void* calculate_row (void* arg)
{
  calculate_row_data *data = (calculate_row_data*) arg;

	double star; /* four times center value minus 4 neigh.b values */
  double residuum; /* residuum of current iteration */

  int tid                       = data->tid;
  int N                         = data->N;
  int rows_per_thread           = data->rows_per_thread;
  double **Matrix_In            = data->Matrix_In;
  double **Matrix_Out           = data->Matrix_Out;
  double fpisin                 = data->fpisin;
  double pih                    = data->pih;
  int term_iteration            = data->term_iteration;
  struct options const* options = data->options;
  double *maxresiduum           = data->maxresiduum;
  pthread_mutex_t *mutex        = data->mutex;

  int start_row = tid * rows_per_thread + 1;
  int end_row = (tid + 1) * rows_per_thread + 1;
  end_row = end_row < N ? end_row : N;

  /* Debug message */
  /* printf("[TID = %d] Starting calculation from row %d to %d\n", tid, start_row, end_row - 1); */

  for (int i = start_row; i < end_row; i++)
  {
    double fpisin_i = 0.0;

    if (options->inf_func == FUNC_FPISIN)
    {
      fpisin_i = fpisin * sin(pih * (double)i);
    }

    /* over all columns */
    for (int j = 1; j < N; j++)
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
        pthread_mutex_lock (mutex);
        *maxresiduum = (residuum < *maxresiduum) ? *maxresiduum : residuum;
        pthread_mutex_unlock (mutex);
      }

      Matrix_Out[i][j] = star;
    }
  }
  free(data);
  pthread_exit(NULL);
}

static
void
calculate (struct calculation_arguments const* arguments, struct calculation_results* results, struct options const* options)
{
  int i;
	int m1, m2;                                 /* used as indices for old and new matrices */
	double maxresiduum;                         /* maximum residuum value of a slave in iteration */

	int const N = arguments->N;
	double const h = arguments->h;

	double pih = 0.0;
	double fpisin = 0.0;

  pthread_t threads[options->number];
  int rows_per_thread = (int) ceil(1.0 * N / options->number);

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

  /* Mutex for maxresiduum */
  pthread_mutex_t mutex = PTHREAD_MUTEX_INITIALIZER;

	while (term_iteration > 0)
	{
    /* Debug message */
    /* printf("Iteration %d ; N = %d; rows_per_thread = %d \n", term_iteration, N, rows_per_thread); */

		double** Matrix_Out = arguments->Matrix[m1];
		double** Matrix_In  = arguments->Matrix[m2];

		maxresiduum = 0;

    /* loop over all threads */
		for (i = 0; i < (int) options->number; i++)
		{
      calculate_row_data *data = (calculate_row_data*) malloc(sizeof(calculate_row_data));
      data->tid = i;
      data->Matrix_In = Matrix_In;
      data->Matrix_Out = Matrix_Out;
      data->fpisin = fpisin;
      data->pih = pih;
      data->term_iteration = term_iteration;
      data->maxresiduum = &maxresiduum;
      data->options = options;
      data->rows_per_thread = rows_per_thread;
      data->N = N;
      data->mutex = &mutex;
      pthread_create(&threads[i], NULL, calculate_row, (void*) data);
		}

    /* joining all threads */
		for (i = 0; i < (int) options->number; i++)
		{
      pthread_join(threads[i], NULL);
		}

		results->stat_iteration++;
		results->stat_precision = maxresiduum;

		/* exchange m1 and m2 */
		i = m1;
		m1 = m2;
		m2 = i;

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
DisplayMatrix (struct calculation_arguments* arguments, struct calculation_results* results, struct options* options)
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

	allocateMatrices(&arguments, &options);
	initMatrices(&arguments, &options);

	gettimeofday(&start_time, NULL);
	calculate(&arguments, &results, &options);
	gettimeofday(&comp_time, NULL);

	displayStatistics(&arguments, &results, &options);
	DisplayMatrix(&arguments, &results, &options);

	freeMatrices(&arguments);

	return 0;
}