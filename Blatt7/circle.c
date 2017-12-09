#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <mpi.h>

#define ROOT_PID 0

int*
init (int size, int rank)
{
  // Todo
  int* buf = malloc(sizeof(int) * size);

  //srand(time(NULL));
  srand(time(NULL)+rank+1);


  for (int i = 0; i < size; i++)
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

void printArray(int *buf, int size, int rank, int num_procs, int N){
  int curr_size[1];
  int curr_buf[N / num_procs + 1];	//maybe last not used
  MPI_Status curr_status;
  /* NOT ALLOWED to receive all bufs and print -> collect one at a time and print*/
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
      curr_size[0] = (j < N % num_procs) ? N / num_procs + 1 : N / num_procs;
      /* Now get buf */
      //int curr_buf = malloc(sizeof(int) * curr_size[0]);
      MPI_Recv(curr_buf, curr_size[0], MPI_INT, j, 0, MPI_COMM_WORLD, &curr_status);
      for (int k = 0; k < curr_size[0]; k++)
      {
        printf ("rank %d: %d\n", j, curr_buf[k]);
      }
    }
  } else {
    // SLAVES
    MPI_Send(buf, size, MPI_INT, ROOT_PID, 0, MPI_COMM_WORLD);
  }
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
  buf = init(size, rank);

  if (rank == ROOT_PID)
  {
    printf("\nBEFORE\n");
  }
  printArray(buf, size, rank, num_procs, N);

  //circle(buf, rank);

  if (rank == ROOT_PID)
  {
    printf("\nAFTER\n");
  }
  printArray(buf, size, rank, num_procs, N);


  /* wartet bis alle Ausgaben erfolgt sind */
  MPI_Barrier(MPI_COMM_WORLD);
  printf("Rang %d beendet jetzt!\n", rank);

  /* terminate processes */
  MPI_Finalize();
  return EXIT_SUCCESS;
}
