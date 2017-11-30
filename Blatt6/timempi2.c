#include <stdio.h>
#include <mpi.h>
#include <unistd.h>
#include <limits.h>
#include <sys/time.h>
#define TIMESTR_MAX_LEN 300
#define ROOT_PID 0

int main(int argc, char **argv)
{
  int pid, num_procs;

  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &pid);
  MPI_Comm_size(MPI_COMM_WORLD, &num_procs);
  
  int local_usec = INT_MAX;
  int min_usec;

  if (pid == ROOT_PID)
  {
    MPI_Status status;
    char timestr[TIMESTR_MAX_LEN];

    for (int i = 1; i < num_procs; i++)
    {
      MPI_Recv(timestr, TIMESTR_MAX_LEN, MPI_CHAR, i, 0, MPI_COMM_WORLD, &status); 
      printf("%s\n", timestr);
    }
  }
  else
  {
    /* Hostname holen */
    char hostname[HOST_NAME_MAX];
    gethostname(hostname, HOST_NAME_MAX);

    /* Zeit holen */
    struct timeval tv;
    struct timezone tz;
    gettimeofday(&tv, &tz);
    local_usec = tv.tv_usec;

    char timestr[TIMESTR_MAX_LEN];
    sprintf(timestr, "[pid = %d] %s: %ld.%ld", pid, hostname, tv.tv_sec, tv.tv_usec);
  
    /* Ergebnis an dem Process 0 schicken */
    MPI_Send(timestr, TIMESTR_MAX_LEN, MPI_CHAR, 0, 0, MPI_COMM_WORLD);
  }

  /* den kleinsten Mikrosekunden-Anteil berechnen */
  MPI_Reduce(&local_usec, &min_usec, 1, MPI_INT, MPI_MIN, ROOT_PID, MPI_COMM_WORLD); 
  if (pid == ROOT_PID) {
    printf("Kleinster Mikrosekunden-Anteil: %d\n", min_usec);
  }

  /* wartet bis alle Ausgaben erfolgt sind */
  MPI_Barrier(MPI_COMM_WORLD);

  printf("Rang %d beendet jetzt!\n", pid);

  MPI_Finalize();

  return 0;
}
