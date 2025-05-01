#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#define N 11
#define M 15

int a[M][N], b[M][N], c[M][N], d[M][N],tempv[M][N];

int main(int argc, char** argv) {

    int myrank, nprocs,is,ie,nlocal,nlast,globali;
    int is2, ie2, inext, iprev;

    MPI_Status status;
    MPI_Request req1, req2, req3,req4;

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

    if (myrank == 0)
    {
        for (int i = 0; i < M; i++)
        {
            for (int j = 0; j < N; j++)
            {
                a[i][j] = i + 10 * j;
            }
        }
        for (int i = 1; i < M - 1; i++)
        {
            for (int j = 1; j < N - 1; j++)
            {
                b[i][j] = a[i - 1][j] + a[i][j - 1] + a[i + 1][j] + a[i][j + 1];
                printf("%d ", b[i][j]);
            }
            printf("\n");
        }
    }

    nlocal = M / nprocs + 1;
    nlast = M - (nprocs - 1) / nlocal;

    is = myrank * nlocal;

    if (myrank == nprocs - 1)
    {
        ie = is + nlast - 1;
        ie2 = ie;
    }
    else
    {
        ie = is + nlocal - 1;
        ie2 = ie + 1;
    }

    if (myrank == 0)
    {
        is2 = is;
    }
    else
    {
        is2 = is - 1;
    }

    for (int i = is; i <= ie; i++)
    {
        for (int j = 0; j < N; j++)
        {
            c[i][j] = i + 10 * j;
        }
    }

    inext = myrank + 1;
    iprev = myrank - 1;

    //SYMETRICKA KOMUNIKACIA
    //pridaju sa fiktivne procesy pre tie co su na okrajoch
    if(myrank == 0)
    {
        iprev = MPI_PROC_NULL;
    }
    if (myrank == nprocs - 1)
    {
        inext = MPI_PROC_NULL;
    }

    MPI_Isend(&c[ie][0], N, MPI_INT, inext, 7, MPI_COMM_WORLD, &req1);
    MPI_Isend(&c[is][0], N, MPI_INT, iprev, 7, MPI_COMM_WORLD, &req2);
    MPI_Irecv(&c[ie2][0], N, MPI_INT, inext, 7, MPI_COMM_WORLD, &req3);
    MPI_Irecv(&c[is2][0], N, MPI_INT, iprev, 7, MPI_COMM_WORLD, &req4);

    MPI_Wait(&req1, &status);
    MPI_Wait(&req2, &status);
    MPI_Wait(&req3, &status);
    MPI_Wait(&req4, &status);


    for (int i = is; i <= ie; i++)
    {
        for (int j = 1; j < N-1; j++)
        {
            tempv[i][j] = c[i - 1][j] + c[i][j - 1] + c[i + 1][j] + c[i][j + 1];
        }
    }

    MPI_Allgather(&tempv[is][0], nlocal*N, MPI_INT, d, nlocal*N, MPI_INT, MPI_COMM_WORLD);

    if (myrank == 0)
    {
        printf("\n***********MPI_SEND_RECV**********\n\n");
        printf("\n******SERIOVE*******\n\n");
        for (int i = 1; i < M - 1; i++)
        {
            for (int j = 1; j < N - 1; j++)
            {
                printf("%d ", b[i][j]);
            }
            printf("\n");
        }

        printf("\n******PARALELNE*******\n\n");
        for (int i = 1; i < M - 1; i++)
        {
            for (int j = 1; j < N - 1; j++)
            {
                printf("%d ", d[i][j]);
            }
            printf("\n");
        }
    }


    MPI_Finalize();
}