#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <mpi.h>

#define ROOT_PID 0

int*
init (int N, int size, int rank)
{
  // Todo
  int* buf = malloc(sizeof(int) * size);

  //srand(time(NULL));
  srand(time(NULL)+rank+1);


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
  buf = init(N, size, rank);

  /* NOT ALLOWED to receive all bufs and print -> collect one at a time and print*/
  MPI_Barrier(MPI_COMM_WORLD);

  int curr_size[1];
  int curr_buf[N / num_procs + 1];	//maybe last not used
  MPI_Status curr_status;
  if (rank == ROOT_PID)
  {
    printf("\nBEFORE\n");
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
      //int curr_buf = malloc(sizeof(int) * curr_size[0]);
      MPI_Recv(curr_buf, curr_size[0], MPI_INT, j, 0, MPI_COMM_WORLD, &curr_status);
      for (int k = 0; k < curr_size[0]; k++)
      {
        printf ("rank %d: %d\n", j, curr_buf[k]);
      }
    }
  } else {
    // SLAVES
    int size_buf[1];
    size_buf[0] = size;
    MPI_Send(size_buf, 1, MPI_INT, ROOT_PID, 0, MPI_COMM_WORLD);

    for (int k = 0; k < size; k++)
    {
      printf ("inside rank %d: %d\n", rank, buf[k]);
    }

    printf ("rank %d: mem @%p\n", rank, (void*)&buf);
    MPI_Send(buf, size, MPI_INT, ROOT_PID, 0, MPI_COMM_WORLD);
  }

  //DEBUG TODO
  /* wartet bis alle Ausgaben erfolgt sind */
  MPI_Barrier(MPI_COMM_WORLD);
  printf("Rang %d beendet jetzt!\n", rank);
  /* terminate processes */
  MPI_Finalize();
  return EXIT_SUCCESS;

  circle(buf, rank);

  printf("\nAFTER\n");

  for (int j = 0; j < N; j++)
  {
    printf ("rank %d: %d\n", rank, buf[j]);
  }

  return EXIT_SUCCESS;
}
