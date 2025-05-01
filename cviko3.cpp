#include<iostream>
#include<stdlib.h>
#include<stdio.h>
#include<mpi.h>
#include <cmath>
#define N 20

//maticove nasobenie A[N][N]b[N]=c[N]
int main(int argc, char* argv[])
{
	int nprocs, myrank, istart, iend, nlocal, nlast;
	double a[N][N],b[N],c[N],d[N],temp[N];
	srand(time(0));

	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
	MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

	if (myrank == 0)
	{
		for (int i = 0; i < N; i++)
		{
			for (int j = 0; j < N; j++)
			{
				a[i][j] = 20*rand()/(RAND_MAX - 10.);
			}
			b[i] = 100*rand()/RAND_MAX;
		}

		for (int i = 0; i < N; i++)
		{
			c[i] = 0;
			for (int j = 0; j < N; j++)
			{
				c[i] += a[i][j] * b[j];
			}
		}
	}

	nlocal = (int)(N / nprocs+1 );
	nlast = N - nlocal * (nprocs - 1);
	MPI_Bcast(a, N*N, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Bcast(b, N, MPI_DOUBLE, 0, MPI_COMM_WORLD);


	istart = myrank * nlocal;
	if (myrank == (nprocs - 1))
	{
		iend = N - 1;
	}
	else
	{
		iend = istart + nlocal - 1;
	}

	//nasobenie
	for (int i = istart; i <= iend; i++)
	{
		temp[i] = 0;
		for (int j = 0; j < N; j++)
		{
			temp[i] += a[i][j] * b[j];
		}
	}
	MPI_Gather(&temp[istart], nlocal, MPI_DOUBLE, d, nlocal, MPI_DOUBLE, 0, MPI_COMM_WORLD);

	if (myrank == 0)
	{
		for (int i = 0; i < N; i++)
		{
			printf("pred: %f\t po: %f\n", c[i], d[i]);
		}
	}

	MPI_Finalize();
}


//hladanie minima, maxima a priemeru ako cez celok a pri paralelizacii- rozdelim vektor a na casti a ratam lokalne minima/maxima a potom priemer zo sum ktore boli zistene lokalne
//int main(int argc, char* argv[])
//{
//	int nprocs, myrank;
//	int a[N], sendmax[2],recvmax[2],sendmin[2],recvmin[2], istart, iend,nlocal,nlast;
//	int glob_max=INT_MIN, glob_max_loc, glob_min= INT_MAX, glob_min_loc,tmp,sum;
//	double glob_priemer=0.0, glob_std=0;
//
//
//	MPI_Init(&argc, &argv);
//	MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
//	MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
//
//	if (myrank == 0)
//	{
//		srand(time(0));
//		sum = 0;
//		for (int i = 0; i < N; i++)
//		{
//			a[i] = rand() % 20001-10000;
//			sum += a[i];
//			if (a[i] > glob_max)
//			{
//				glob_max = a[i];
//				glob_max_loc = i;
//			}
//			if (a[i] < glob_min)
//			{
//				glob_min = a[i];
//				glob_min_loc = i;
//			}
//		}
//		glob_priemer = (double)sum/ N;
//	}
//	MPI_Bcast(a, N, MPI_INT, 0, MPI_COMM_WORLD);
//
//
//	printf("\nJEDNOTLIVE STATISTIKY\n");
//	nlocal = (int)(N / nprocs + 1);
//	nlast = N - nlocal * (nprocs - 1);
//
//	istart = myrank * nlocal;
//	if (myrank == (nprocs - 1))
//	{
//		iend = N - 1;
//		printf("%d:\t %d-%d\t %d\n", myrank, istart, iend,nlocal);
//	}
//	else
//	{
//		iend = istart + nlocal - 1;
//		printf("%d:\t %d-%d\t %d\n", myrank, istart, iend,nlast);
//	}
//	
//	int lok_max=INT_MIN, lok_min=INT_MAX, lok_max_loc, lok_min_loc;
//	double lok_priemer;
//
//	sum = 0;
//	for (int i = istart; i <= iend; i++)
//	{
//		sum += a[i];
//		if (a[i] > lok_max)
//		{
//			lok_max = a[i];
//			lok_max_loc = i;
//		}
//		if (a[i] < lok_min)
//		{
//			lok_min = a[i];
//			lok_min_loc = i;
//		}
//	}
//
//	sendmax[0] = lok_max;
//	sendmax[1] = lok_max_loc;
//	sendmin[0] = lok_min;
//	sendmin[1] = lok_min_loc;
//
//	MPI_Reduce(sendmax, recvmax, 1, MPI_2INT, MPI_MAXLOC, 0, MPI_COMM_WORLD);
//	MPI_Reduce(sendmin, recvmin, 1, MPI_2INT, MPI_MINLOC, 0, MPI_COMM_WORLD);
//	MPI_Reduce(&sum, &tmp, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
//	
//	lok_priemer = (double)tmp / N;
//
//	printf("%d:\t Max: %d (%d)\t min: %d (%d)\n\n", myrank, lok_max, lok_max_loc, lok_min, lok_min_loc);
//
//
//	if (myrank == 0)
//	{
//		printf("\nGLOBALNE STATISTIKY PRED DELENIM\n");
//		printf("%d:\t Max: %d (%d)\t min: %d (%d)\t Priemer: %.20f\n\n", myrank, glob_max, glob_max_loc, glob_min, glob_min_loc, glob_priemer);
//		printf("\nREDUCE- STATISTIKY S DELENIM\n");
//		printf("%d:\t Max: %d (%d)\t min: %d (%d)\t priemer: %.20f\n\n", myrank, recvmax[0], recvmax[1], recvmin[0], recvmin[1],lok_priemer);
//	}
//
//	MPI_Finalize();
//}