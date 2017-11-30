#include <stdio.h>
#include <mpi.h>
#include <unistd.h>
#include <limits.h>
#include <sys/time.h>
#define TIMESTR_MAX_LEN 300

int main(int argc, char **argv)
{
  int pid, num_procs;

  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &pid);
  MPI_Comm_size(MPI_COMM_WORLD, &num_procs);
  
  if (pid == 0)
  {
    char timestr[TIMESTR_MAX_LEN];
    MPI_Status status;
    for (int i = 1; i < num_procs; i++) {
      MPI_Recv(timestr, TIMESTR_MAX_LEN, MPI_CHAR, i, 0, MPI_COMM_WORLD, &status); 
      printf("%s\n", timestr);
    }
  }
  else
  {
    /* get hostname */
    char hostname[HOST_NAME_MAX];
    gethostname(hostname, HOST_NAME_MAX);

    /* get time */
    struct timeval tv;
    struct timezone tz;
    gettimeofday(&tv, &tz);

    char timestr[TIMESTR_MAX_LEN];
    sprintf(timestr, "[pid = %d] %s: %ld.%ld", pid, hostname, tv.tv_sec, tv.tv_usec);

    MPI_Send(timestr, TIMESTR_MAX_LEN, MPI_CHAR, 0, 0, MPI_COMM_WORLD);
  }

  /* wartet bis alle Ausgaben erfolgt sind */
  MPI_Barrier(MPI_COMM_WORLD);

  printf("Rang %d beendet jetzt!\n", pid);

  MPI_Finalize();

  return 0;
}
