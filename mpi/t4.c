#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <mpi.h>
#define h 0.1
#define k 1
#define T 0.1
#define L 1
#define EXCHANGE_TAG 33
#define Yard 1000
#define PI 3.14159
#define E 2.71828

void precive(double x);
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
	int num = N/size;
	array = (double*)calloc(num + 2, sizeof(double));
	newarray = (double*)calloc(num + 2, sizeof(double));
	double t = h * h / k;
	int steps = T / t;
	//begin = num * myrank;
	//end = begin + num;
	for (i = 1; i < num; i++)
        {
                array[i] = 1;
                newarray[i] = 1;
        }
	if (myrank == 0)
	{
		array[0] = 0;
		newarray[0] = 0;
	}
	if(myrank == size - 1) 
	{
		array[num + 1] = 0;
		newarray[num + 1] = 0;
	}
	//end = N - 1;
	//end = end + 1;	
	int tm = 0;
	double tou = 0.003;
	double co = (tou * k) / (h * h); 
	for(i = 0; i < steps; i++)
	{
		if(myrank == 0)
		{
           		MPI_Isend(&array[num], 1, MPI_DOUBLE, myrank + 1, EXCHANGE_TAG, MPI_COMM_WORLD, &request);
           		MPI_Irecv(&array[num + 1], 1, MPI_DOUBLE, myrank + 1, EXCHANGE_TAG, MPI_COMM_WORLD, &request);
			MPI_Wait(&request, &Status);
		}
		else
		{
			if(myrank != size - 1)
			{
				MPI_Isend(&array[1], 1, MPI_DOUBLE, myrank - 1, EXCHANGE_TAG, MPI_COMM_WORLD, &request); 
				MPI_Irecv(&array[0], 1, MPI_DOUBLE, myrank - 1, EXCHANGE_TAG, MPI_COMM_WORLD, &request);
		                MPI_Isend(&array[num], 1, MPI_DOUBLE, myrank + 1, EXCHANGE_TAG, MPI_COMM_WORLD, &request);
				MPI_Irecv(&array[num + 1], 1, MPI_DOUBLE, myrank + 1, EXCHANGE_TAG, MPI_COMM_WORLD, &request);
				MPI_Wait(&request, &Status);
			}
			else
			{
				MPI_Isend(&array[1], 1, MPI_DOUBLE, myrank - 1, EXCHANGE_TAG, MPI_COMM_WORLD, &request); 
				MPI_Irecv(&array[0], 1, MPI_DOUBLE, myrank - 1, EXCHANGE_TAG, MPI_COMM_WORLD, &request);
				MPI_Wait(&request, &Status);
			}
		}
		for(j = 1; j <= num; j++)
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
	printf("\n%d ", myrank);
	for (j = 0; j <= num + 1; j++)
	{
		printf("%lf ", array[j]);
	}
	/*if (myrank == 0) {
		for (i = 1; i < size - 1; i++) {			
			MPI_Recv(&array[num * i], num, MPI_DOUBLE, i, 0, MPI_COMM_WORLD, &Status);
		}
		MPI_Recv(&array[num * (size - 1)], num  +  N % size, MPI_DOUBLE, i, 0, MPI_COMM_WORLD, &Status);
		for (i = 0; i < N; i++) {
			printf("%f %f\n", h * i, array[i]);
		}
		for (i = 0; i < N; i++)
		{
			precive(h * i);
		}
	} else {
		if(myrank == size - 1)
		    MPI_Send(&array[num * myrank], num + N % size, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
		else
		    MPI_Send(&array[num * myrank], num, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
	}*/
	printf("\n");
	free(array);
	free(newarray);
	MPI_Finalize();
	return 0;
}
void precive(double x)
{
	int i = 0;
	double sum = 0;
	for(i = 0; i < Yard; i++)
	{
		sum+=(4/PI) 
		* pow
		(E,
			 -(k * PI * PI * (2 * i + 1) * (2 * i + 1) * T / (L * L)
			/(2 * i + 1))
		) 
		* sin (PI * (2 * i + 1) * x / L);
	}
	printf(" \n%lf", sum);
}
