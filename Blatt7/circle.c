#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <mpi.h>

#define ROOT_PID 0

int*
init (int N, int size)
{
  // Todo
  int* buf = malloc(sizeof(int) * size);

  srand(time(NULL));

  for (int i = 0; i < N; i++)
  {
    // Do not modify % 25
    buf[i] = rand() % 25;
  }

  return buf;
}

int*
circle (int* buf, int rank)
{
  // Todo
  return buf;
}

int
main (int argc, char** argv)
{
  char arg[256];
  int N;
  int* buf;
  int size;

  int rank, num_procs;

  if (argc < 2)
  {
    printf("Arguments error!\n");
    return EXIT_FAILURE;
  }

  sscanf(argv[1], "%s", arg);

  // Array length
  N = atoi(arg);

  //TODO check for input errors

  // Todo: myrank
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &num_procs);

  // ´N / num_procs´ splits evenly; if rest exists -> spread until no rest (<=> for rank < rest)
  size = (rank < N % num_procs) ? N / num_procs + 1 : N / num_procs;
  buf = init(N, size);

  /* NOT ALLOWED to receive all bufs and print -> collect one at a time and print*/
  printf("\nBEFORE\n");

  int *curr_buf;
  int curr_size;
  MPI_Status curr_status;
  if (rank == ROOT_PID)
  {
    // MASTER
    for (int i = 0; i < size; i++)
    {
      printf ("rank %d: %d\n", rank, buf[i]);
    }
    for (int j = 1; j < num_procs; j++)
    {
      /* Get size first */
      MPI_Recv(curr_size, 1, MPI_INT, j, 0, MPI_COMM_WORLD, &curr_status);
      /* Now get buf */
      MPI_Recv(curr_buf, curr_size, MPI_INT, j, 0, MPI_COMM_WORLD, &curr_status);
      for (int k = 0; k < size; k++)
      {
        printf ("rank %d: %d\n", j, buf[k]);
      }
    }
  } else {
    // SLAVES -- implicit barrier -> no barrier after init needed
    MPI_Send(size, 1, MPI_INT, ROOT_PID, 0, MPI_COMM_WORLD);
    MPI_Send(buf, size, MPI_INT, ROOT_PID, 0, MPI_COMM_WORLD);
  }

  //DEBUG TODO
  return EXIT_SUCCESS;

  circle(buf, rank);

  printf("\nAFTER\n");

  for (int j = 0; j < N; j++)
  {
    printf ("rank %d: %d\n", rank, buf[j]);
  }

  return EXIT_SUCCESS;
}
