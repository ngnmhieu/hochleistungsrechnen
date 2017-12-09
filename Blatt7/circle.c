#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <mpi.h>

#define ROOT_PID 0

int getSize(int N, int num_procs, int rank, int itteration)
{
  // ´N / num_procs´ splits evenly; if rest exists -> spread until no rest (<=> for rank < rest)
  // return (rank < N % num_procs) ? N / num_procs + 1 : N / num_procs;
  // adjusted for current ittreration (<=> shift in "rank" for indexing)
  // added itteration * num_procs because modulo with negatives does not work
  //      printf ("adj rank %d: %d it(%d)\n", rank, (rank - itteration + ((itteration) * num_procs)) % num_procs, itteration);

  return (((rank - itteration + (itteration * num_procs)) % num_procs) < N % num_procs) ? N / num_procs + 1 : N / num_procs;
}

int*
init (int size, int rank)
{
  // Todo
  int* buf = malloc(sizeof(int) * size);

  srand(time(NULL)+rank+1);

  for (int i = 0; i < size; i++)
  {
    // Do not modify % 25
    buf[i] = rand() % 25;
  }

  return buf;
}

void printArray(int *buf, int size, int rank, int num_procs, int N, int itteration){
  int curr_size[1];
  int curr_buf[N / num_procs + 1];	//maybe last not used
  MPI_Status curr_status;
  /* NOT ALLOWED to receive all bufs and print -> collect one at a time and print*/
  if (rank == ROOT_PID)
  {
    // MASTER
    for (int i = 0; i < getSize(N, num_procs, ROOT_PID, itteration) ; i++)
    {
      printf ("rank %d: %d\n", rank, buf[i]);
    }
    for (int j = 1; j < num_procs; j++)
    {
      /* Get size first */
      curr_size[0] = getSize(N, num_procs, j, itteration);
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


/**
 * Gibt den Nachfolger zurueck
 * 0 wird zurueck gegeben falls rank == num_proc - 1
 */
int getNextRank(int rank, int num_proc) {
  return (rank + 1) % num_proc;
}

int getPrevRank(int rank, int num_proc) {
  return (rank == 0) ? num_proc - 1 : rank - 1;
}

int* circle (int* buf, int rank, int num_procs, int N)
{
  int firstVal;

  // den Wert fuer die Abbruchbedingung von P_0 -> P_N-1 senden
  if (rank == 0) {
    MPI_Send(buf, 1, MPI_INT, num_procs - 1, 0, MPI_COMM_WORLD); 
  }
  if (rank == num_procs - 1) {
    MPI_Recv(&firstVal, 1, MPI_INT, ROOT_PID, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
  }

  MPI_Barrier(MPI_COMM_WORLD);

  int running = 1;
  int nextRank, prevRank, bufSize, prevBufSize;
  for (int i = 0; running; i++) {

    // DEBUG
    if (rank == 0) {
      printf("===== Iteration %d ====\n", i);
    }
    MPI_Barrier(MPI_COMM_WORLD);
    // END DEBUG

    // Daten an Nachfolger verschicken
    nextRank = getNextRank(rank, num_procs);
    bufSize = getSize(N, num_procs, rank, i);
    printf("[PID = %d] Sending %d ints to process %d.\n", rank, bufSize, nextRank); 
    MPI_Send(buf, bufSize, MPI_INT, nextRank, 0, MPI_COMM_WORLD); 

    // Daten aus Vorgaenger empfangen
    prevRank = getPrevRank(rank, num_procs);
    prevBufSize = getSize(N, num_procs, prevRank, i);
    int* tempBuf = buf; // save pointer to buf
    buf = (int*) malloc(sizeof(int) * prevBufSize);
    printf("[PID = %d] Receiving %d ints from process %d.\n", rank, prevBufSize, prevRank); 
    MPI_Recv(buf, prevBufSize, MPI_INT, prevRank, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    free(tempBuf);

    MPI_Barrier(MPI_COMM_WORLD);

    // DEBUG
    printArray(buf, prevBufSize, rank, num_procs, N, 0);
    // END DEBUG

    if (rank == num_procs - 1) {
      if (buf[0] == firstVal){
        running = 0;
      }
    }

    MPI_Bcast(&running, 1, MPI_INT,num_procs-1, MPI_COMM_WORLD);

    MPI_Barrier(MPI_COMM_WORLD);
  }

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

  if (N <= 0)
  {
    printf("N must be positive!\n");
    return EXIT_FAILURE;
  }

  // Todo: myrank
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &num_procs);

  //size = (rank < N % num_procs) ? N / num_procs + 1 : N / num_procs;
  size = getSize(N, num_procs, rank, 0);
  buf = init(size, rank);

  if (rank == ROOT_PID) {
    printf("\nBEFORE\n");
  }
  printArray(buf, size, rank, num_procs, N, 0);
  MPI_Barrier(MPI_COMM_WORLD);

  circle(buf, rank, num_procs, N);

  MPI_Barrier(MPI_COMM_WORLD);
  if (rank == ROOT_PID) {
    printf("\nAFTER\n");
  }
  printArray(buf, size, rank, num_procs, N, 1);


  /* wartet bis alle Ausgaben erfolgt sind */
  MPI_Barrier(MPI_COMM_WORLD);

  printf("Rang %d beendet jetzt!\n", rank);

  /* terminate processes */
  MPI_Finalize();
  return EXIT_SUCCESS;
}
