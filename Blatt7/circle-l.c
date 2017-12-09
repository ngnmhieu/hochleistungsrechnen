#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>

#include <mpi.h>

// Function to get actual length of buffer depending on rank, step and mod
int getLength(int rank, int nprocs, int N, int step, int mod){
    int base = N/nprocs;
    int rest = N%nprocs;
    int length;
    int actual = (nprocs+rank-step)%nprocs + mod;
    // Get Last process if first process
    if (actual == -1) {
        actual = nprocs -1;
    }
    if (actual < rest){
        length=base+1;
    }
    else {
        length = base;
    }
    return length;
}


int* init (int rank, int nprocs, int N)
{
    // Initializing buffer with length_per_rank + 1 for maximal length
    int length_per_rank = N/nprocs;
    int* buf = (int*) malloc(sizeof(int) * length_per_rank+1);
    int actual_length = getLength(rank, nprocs, N, 0, 0);

    // Adding Rank do get unique random numbers
    srand(time(NULL)+rank+1);

    // Initializing array
    for (int i = 0; i < actual_length; i++) {
        buf[i] = rand() % 25; //do not modify %25
    }
    return buf;
}


int* circle (int* buf, int rank, int N, int nprocs, int* glob_step)
{
    int backup;
    int step = 0;
    // Sending target value to last Process
    if (rank == 0){
        MPI_Send(buf, 1, MPI_INT, nprocs-1, 0, MPI_COMM_WORLD);
    }
    else if (rank == nprocs-1){
        MPI_Recv(&backup, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }
    
    int length_per_rank = N/nprocs+1;

    int running = 1;
    while(running){
        int* temp_buf = (int*) malloc(sizeof(int) * length_per_rank);
        int actual_length_send = getLength(rank, nprocs, N, step, 0);
        int actual_length_resv = getLength(rank, nprocs, N, step, -1);
        // Sending buffer to next process
        if (rank == nprocs-1){
            MPI_Send(buf, actual_length_send, MPI_INT, 0, 0, MPI_COMM_WORLD);
        }
        else { 
            MPI_Send(buf, actual_length_send, MPI_INT, rank+1, 0, MPI_COMM_WORLD);
        }
        // Getting buffer from previous process
        if (rank == 0) {
            MPI_Recv(temp_buf, actual_length_resv, MPI_INT, nprocs-1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }
        else { 
            MPI_Recv(temp_buf, actual_length_resv, MPI_INT, rank-1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }

        // Wait for all processes to receive stuff
        MPI_Barrier(MPI_COMM_WORLD);
        // Copy temp buff to normal buffer
        memcpy(buf, temp_buf, sizeof(int) * actual_length_resv);
        // Free temp buff
        free(temp_buf);

        // Check if requirements for stopping are given, if true running is set to 0
        if (rank == nprocs-1) {
            if (buf[0] == backup){
                running = 0;
            }
        }
        // nprocs-1 sends his running state to all other processes.
        // This kills the while loop for all processes, if requirements are fulfilled
        MPI_Bcast(&running, 1, MPI_INT,nprocs-1, MPI_COMM_WORLD);
        
        MPI_Barrier(MPI_COMM_WORLD);
        // Incrementing step to get 
        step += 1;
    }
    *glob_step = step;
    return buf;
}

void print_sorted(int rank, int nprocs, int N, int* buf, int glob_step){
    // Allocating temp_buf with maximal size 
    int length_per_rank = N/nprocs + 1;
    // Rank 0 getting buffers from all processes
    if (rank == 0) {
        printf("Constellation: ");
        // Print buffer for rank 0
        int length_0 = getLength(rank, nprocs, N, glob_step, 0);
        for (int i = 0;  i<length_0; i++){
            printf("%i ", buf[i]);
        }
        for (int i=1 ; i < nprocs; i++) {
        int* temp_buf = (int*) malloc(sizeof(int) * length_per_rank);
            int length = getLength(i, nprocs, N, glob_step, 0);
            MPI_Recv(temp_buf, length, MPI_INT, i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            // Prints buffer for all other processes
            for (int i = 0; i<length; i++){
                printf("%i ", temp_buf[i]);
            }
        }
        printf("\n");
    }
    // All other processes send their buffer to rank 0 for printing
    else {
        int length = getLength(rank, nprocs, N, glob_step, 0);
        MPI_Send(buf, length, MPI_INT, 0, 0, MPI_COMM_WORLD);
    }
    MPI_Barrier(MPI_COMM_WORLD);
}

int main (int argc, char** argv)
{
    int rank, nprocs;
    // Init functions for MPI
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

    char arg[256];
    int N;
    int* buf;
    int glob_step = 0;

    if (argc < 2) {
        if (rank == 0){
            printf("Please specify the array length\n");
            return EXIT_FAILURE;
        }
    }

    sscanf(argv[1], "%s", arg);

    //array length
    N = atoi(arg);
    buf = init(rank, nprocs, N);

    // Check if there are more processes then array has length
    if (N < nprocs){
        if (rank == 0){
            printf("More processes than array length\n");
            return EXIT_FAILURE;
        }
    }

    print_sorted(rank, nprocs, N, buf, glob_step);
    
    circle (buf,rank,N,nprocs, &glob_step);

    print_sorted(rank, nprocs, N, buf, glob_step);

    MPI_Finalize();
    return EXIT_SUCCESS;
}