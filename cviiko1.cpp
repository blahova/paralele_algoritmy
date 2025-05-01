#include<iostream>
#include<stdlib.h>
#include<stdio.h>
#include<mpi.h>


int main(int argc, char* argv[])
{
	int nprocs, myrank;

	//MPI BROADCAST
	// broadcastuje hodnoty z pola a do premennej (pola) s tym istym nazvom z procesu 0 do vsetkych procesov
	// int i, a[4];
	//MPI_Init(&argc, &argv);
	//MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
	//MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
	//for (i = 0; i < 4; i++)
	//{
	//	if (myrank == 0) a[i] = i;
	//	else a[i] = 0;
	//}
	//printf("%d-predMPI. %d %d %d %d \n", myrank, a[0], a[1], a[2], a[3]);
	//MPI_Bcast(&a[2], 2, MPI_INT, 2, MPI_COMM_WORLD);
	//printf("%d-poMPI. %d %d %d %d \n", myrank, a[0], a[1], a[2], a[3]);


	//MPI GATHER
	// zbiera udaje z premennej isend od ostatnych procesov. Da sa zadat do ktoreho procesu sa to bude zbierat.
	// je preddefinovane poradie, od nulteho, vsetky procesy musia poslat rovnaky pocet udajov
	// 
	// int i, isend[2], irecv[32];
	// 	MPI_Init(&argc, &argv);
	//MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
	//MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
	//isend[0] = myrank + 1;
	//isend[1] = myrank + 1;
	//MPI_Gather(isend, 2, MPI_INT, irecv, 2, MPI_INT, 0, MPI_COMM_WORLD); //(aka premenna sa posiela, kolko z tej premennej ci jeden prvok alebo viac, typ, kam sa uklada, kolko sa uklada, akeho typu, na ktory proces)
	//if (myrank == 0)
	//{
	//	for (i = 0; i < 2*nprocs; i++)
	//	{
	//		printf("%d\n", irecv[i]);
	//	}
	//	printf("\n");
	//}
	//MPI_Finalize();

	//MPI GATHERV
	//da sa specifikovat z ktoreho procesu pojde kolko prvkov aaj kam sa ulozia
	//int i, isend[5], irecv[14];
	//int  displs[4] = { 5,3, 10,0 };
	//int recvcounts[4] = { 5,2,4,3 };
	//int sendcount;
	//MPI_Init(&argc, &argv);
	//MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
	//MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
	//sendcount = recvcounts[myrank];
	//for (int j = 0; j < sendcount; j++)	isend[j] = myrank;
	//
	//MPI_Gatherv(isend, sendcount, MPI_INT, irecv, recvcounts, displs, MPI_INT, 0, MPI_COMM_WORLD);
	//if (myrank == 0)
	//{
	//	for (i = 0; i < 14; i++)
	//	{
	//		printf("%d\t", irecv[i]);
	//	}
	//	printf("\n");
	//}
	//MPI_Finalize();


	//MPI REDUCE SUM
	//rozdelime si pole na mensie casti, zratame sumy na kazdom a potom zratame sumu dokopy
	//int i, istart, iend;
	//int a[10], sum, tmp;
	//MPI_Init(&argc, &argv);
	//MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
	//MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
	//
	//istart = myrank * 3 + 1;
	//iend = istart + 2;
	//sum = 0;

	//for (i = istart; i <= iend; i++)
	//{
	//	a[i] = i;
	//	sum += a[i];
	//}

	//printf("suma na P%d: %d\n", myrank,sum);
	//MPI_Reduce(&sum, &tmp, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);

	//sum = tmp;
	//if (myrank == 0)
	//{
	//	printf("%d: celkova suma: %d\n", myrank,sum);
	//}

	//MPI_Finalize();



	//MPI REDUCE MAXLOC
	//rozdelime pole na viac casti, hladame v kazdom lokalne maximum, ulozim si ho spolu s poziciou, potom hladam globalne z lokalnych
	int i, istart, iend;
	int a[9] = { 12,15,2,20,8,3,7,24,52 };
	int max, loc, isend[2], irecv[2];
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
	MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

	istart = myrank * 3 + 1;
	iend = istart + 2;
	max = -999;

	for (i = istart; i <= iend; i++)
	{
		if (a[i] > max)
		{
			max = a[i];
			loc = i;
		}
	}
	isend[0] = max;
	isend[1] = loc;
	printf("max na P%d: %d na pozicii %d\n", myrank, max,loc);
	MPI_Reduce(isend, irecv, 1, MPI_2INT, MPI_MAXLOC, 0, MPI_COMM_WORLD);

	if (myrank == 0)
	{
		printf("%d: celkove max: %d na pozicii %d\n", myrank, irecv[0],irecv[1]);
	}

	MPI_Finalize();
}