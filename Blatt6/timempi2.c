#include <stdio.h>
#include <mpi.h>
#include <unistd.h>		// library for gethostname()
#include <limits.h>		// constants, eg for INT_MAX
#include <sys/time.h>		// library for gettimeofday()
#define TIMESTR_MAX_LEN 300
#define ROOT_PID 0

int main(int argc, char **argv)
{
  int pid, num_procs;
  int loc_usec = INT_MAX;
  int min_usec;

  /* Spawning slave processes, based on -np flag */
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &pid);		// get rank for every process
  MPI_Comm_size(MPI_COMM_WORLD, &num_procs);	// get number of processes

  /* master */
  if (pid == ROOT_PID)
  {
    char timestr[TIMESTR_MAX_LEN];
    MPI_Status status;

    for (int i = 1; i < num_procs; i++) {
      /* receive data from process #i -> ordered */
      MPI_Recv(timestr, TIMESTR_MAX_LEN, MPI_CHAR, i, 0, MPI_COMM_WORLD, &status);
      printf("%s\n", timestr);
    }
  }
  /* slaves */
  else
  {
    /* get hostname */
    char hostname[HOST_NAME_MAX];		// HOST_NAME_MAX defined with value 64 on Linux
    gethostname(hostname, HOST_NAME_MAX);

    /* get time */
    struct timeval tv;
    struct timezone tz;
    gettimeofday(&tv, &tz);
    loc_usec = tv.tv_usec;

    char timestr[TIMESTR_MAX_LEN];
    /* write to buffer */
    sprintf(timestr, "[pid = %d] %s: %ld.%ld", pid, hostname, tv.tv_sec, tv.tv_usec);

    /* send data to master */
    MPI_Send(timestr, TIMESTR_MAX_LEN, MPI_CHAR, ROOT_PID, 0, MPI_COMM_WORLD);
  }

  /* get min microsec */
  MPI_Reduce(&loc_usec, &min_usec, 1, MPI_INT, MPI_MIN, ROOT_PID, MPI_COMM_WORLD);
  if (pid == ROOT_PID) {
    printf("Kleinster Mikrosekunden-Anteil: %d\n", min_usec);
  }

  /* wartet bis alle Ausgaben erfolgt sind */
  MPI_Barrier(MPI_COMM_WORLD);

  printf("Rang %d beendet jetzt!\n", pid);

  /* terminate processes */
  MPI_Finalize();

  return 0;
}
