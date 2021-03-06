#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <mpi.h>
#define h 0.1
#define k 1
#define T 0.1
#define L 1
#define EXCHANGE_TAG 33

int main(int argc, char* argv[])
{
	int myrank, size, begin, end;
 	MPI_Status Status;
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
	MPI_Request request;
	int N = L/h + 1;
	double  t0, total, tmp1, tmp2;
	int i = 0;
	int j = 0;
	int s = 1 / (h * size);
	//double *array[2];
	double* array;
	double* newarray;
	array = (double*)calloc(N, sizeof(double));
	newarray = (double*)calloc(N, sizeof(double));
	double t = h * h / k;
	int steps = T / t;
	int num = N / size;
	begin = num * myrank;
	end = begin + num;
	if (myrank == 0)
	{
		array[0] = 0;
		newarray[0] = 0;
		begin = 1;
	}
	if(myrank == size - 1)
	{
		array[s] = 0;
		newarray[s] = 0;
		end = N - 1;
	}	
	for (i = begin; i < end; i++)
	{
		array[i] = 1;
		newarray[i] = 1;
	}
	int tm = 0;
	double tou = 0.003;
	double co = (tou * k) / (h * h); 
 	for(i = 0; i < steps; i++)
	{
		if(myrank == 0)
		{
           		MPI_Isend(&array[end - 1], 1, MPI_DOUBLE, myrank + 1, EXCHANGE_TAG, MPI_COMM_WORLD, &request);
           		MPI_Irecv(&array[end], 1, MPI_DOUBLE, myrank + 1, EXCHANGE_TAG, MPI_COMM_WORLD, &request);
			MPI_Wait(&request, &Status);
		}
		else
		{
			if(myrank != size - 1)
			{
				MPI_Isend(&array[begin], 1, MPI_DOUBLE, myrank - 1, EXCHANGE_TAG, MPI_COMM_WORLD, &request); 
				MPI_Irecv(&array[begin - 1], 1, MPI_DOUBLE, myrank - 1, EXCHANGE_TAG, MPI_COMM_WORLD, &request);
		                MPI_Isend(&array[end - 1], 1, MPI_DOUBLE, myrank + 1, EXCHANGE_TAG, MPI_COMM_WORLD, &request);
				MPI_Irecv(&array[end], 1, MPI_DOUBLE, myrank + 1, EXCHANGE_TAG, MPI_COMM_WORLD, &request);
				MPI_Wait(&request, &Status);
			}
			else
			{
				MPI_Isend(&array[begin], 1, MPI_DOUBLE, myrank - 1, EXCHANGE_TAG, MPI_COMM_WORLD, &request); 
				MPI_Irecv(&array[begin - 1], 1, MPI_DOUBLE, myrank - 1, EXCHANGE_TAG, MPI_COMM_WORLD, &request);
				MPI_Wait(&request, &Status);
			}
		}
		for(j = begin; j < end; j++)
		{
			if(i % 2 == 0)
			{
			    newarray[j] = array[j] + co * (array[j + 1] - 2.0 * array[j] + array[j - 1]);
			}
			else
			{
			    array[j] = newarray[j] + co * (newarray[j + 1] - 2.0 * newarray[j] + newarray[j - 1]);
			}
		}
	}
	if (i % 2 == 0)
	{
	    memcpy(array, newarray, N);
	}
	if (myrank == 0) {
		for (i = 1; i < size - 1; i++) {			
			MPI_Recv(&array[num * i], num, MPI_DOUBLE, i, 0, MPI_COMM_WORLD, &Status);
		}
		MPI_Recv(&array[num * (size - 1)], num  +  N % size, MPI_DOUBLE, i, 0, MPI_COMM_WORLD, &Status);
		for (i = 0; i < N; i++) {
			printf("%f %f\n", h * i, array[i]);
		}
	} else {
		if(myrank == size - 1)
		    MPI_Send(&array[num * myrank], num + N % size, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
		else
		    MPI_Send(&array[num * myrank], num, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
	}
	printf("\n");
	free(array);
	free(newarray);
	MPI_Finalize();
	return 0;
}
