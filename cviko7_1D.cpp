#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#define N 21

int a[N], b[N], c[N], d[N], tempv[N];
time_t t0, t1, t2;

int main(int argc, char** argv) {

    int myrank, nprocs,is,ie,nlocal,nlast;
    int is2, ie2, inext, iprev;

    MPI_Status status;
    MPI_Request req1, req2, req3,req4;

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

    if (myrank == 0)
    {
        for (int i = 0; i < N; i++)
        {
            b[i] = i + 1;
        }
        for (int i = 1; i < N - 1; i++)
        {
            a[i] = b[i - 1] + b[i + 1];
            printf("%d: \t %d\n", i, a[i]);
        }
    }

    nlocal = N / nprocs + 1;
    nlast = N - (nprocs - 1) / nlocal;

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
        c[i] = i + 1;
    }

    inext = myrank + 1;
    iprev = myrank - 1;


    //NESYMETRICKA KOMUNIKACIA
    //if (myrank == 0)
    //{
    //    MPI_Isend(c + ie, 1, MPI_INT, inext, 7, MPI_COMM_WORLD, &req1);
    //    MPI_Irecv(c + ie2, 1, MPI_INT, inext, 7, MPI_COMM_WORLD, &req2);
    //}
    //if ((myrank > 0) && (myrank < nprocs - 1))
    //{
    //    MPI_Isend(c + ie, 1, MPI_INT, inext, 7, MPI_COMM_WORLD, &req1);
    //    MPI_Isend(c + is, 1, MPI_INT, iprev, 7, MPI_COMM_WORLD, &req2);
    //    MPI_Irecv(c + ie2, 1, MPI_INT, inext, 7, MPI_COMM_WORLD, &req3);
    //    MPI_Irecv(c + is2, 1, MPI_INT, iprev, 7, MPI_COMM_WORLD, &req4);
    //}
    //if (myrank == nprocs-1)
    //{
    //    MPI_Isend(c + is, 1, MPI_INT, iprev, 7, MPI_COMM_WORLD, &req1);
    //    MPI_Irecv(c + is2, 1, MPI_INT, iprev, 7, MPI_COMM_WORLD, &req2);
    //}

    //MPI_Wait(&req1, &status);
    //MPI_Wait(&req2, &status);

    //if (myrank > 0 && myrank < nprocs - 1)
    //{
    //    MPI_Wait(&req3, &status);
    //    MPI_Wait(&req4, &status);
    //}

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

    MPI_Isend(c + ie, 1, MPI_INT, inext, 7, MPI_COMM_WORLD, &req1);
    MPI_Isend(c + is, 1, MPI_INT, iprev, 7, MPI_COMM_WORLD, &req2);
    MPI_Irecv(c + ie2, 1, MPI_INT, inext, 7, MPI_COMM_WORLD, &req3);
    MPI_Irecv(c + is2, 1, MPI_INT, iprev, 7, MPI_COMM_WORLD, &req4);

    MPI_Wait(&req1, &status);
    MPI_Wait(&req2, &status);
    MPI_Wait(&req3, &status);
    MPI_Wait(&req4, &status);


    for (int i = is; i <= ie; i++)
    {
        tempv[i] = c[i - 1] + c[i + 1];
    }

    MPI_Allgather(tempv + is, nlocal, MPI_INT, d, nlocal, MPI_INT, MPI_COMM_WORLD);

    if (myrank == 0)
    {
        printf("\n***********MPI_SEND_RECV**********\n\n");
        for (int i = 1; i < N-1; i++)
        {
            printf("%d\t%d\t%d\n", i, a[i], d[i]);
        }
    }


    MPI_Finalize();
}