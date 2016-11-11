#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <mpi.h>
#define h 0.1
#define k 1
#define T 1000

int main(int argc, char* argv[])
{
	int myrank, size;
 	MPI_Status Status;
	MPI_Init(&argc, &argv);
	MPI_Request Request;
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
	double begin, end, total, tmp1, tmp2;
	int i = 0;
	double* array = calloc(size/h, sizeof(double));
	double* newarray = calloc(size/h, sizeof(double));
	double t = h * h / k;
	int s = int(size / h);
	for (i = 1; i < size/h; i++)
	{
		array[i] = 1;	
	}
	begin = MPI_WTIME();
	if (myrank == 0)
	{
		array[0] = 0;
	}
	else
	{
		if(myrank == size - 1)
		{
			array[size/h - 1] = 0;
		}	
	}
	while(total = begin - MPI_Wtime() < T)
	{
		if(myrank == 0)
		{
			tmp1 = 0;
                        MPI_Isendrecv(&array[size/h - 1], 1, MPI_DOUBLE, myrank + 1, myrank + 1, &tmp2, 1, MPI_DOUBLE, myrank + 1, myrank + 1, MPI_COMM_WORLD, &Status);
                        newarray[0] = array[0] + k*t/(h * h) * (array[1] - 2 * array[0] + tmp1);
                        newarray[size/h - 1] = array[size/h - 1] + k*t/(h * h) * (tmp2 - 2 * array[size/h - 1] + array[size/h - 2]);
		}
		else
		{
			if(myrank != size - 1)
			{
				MPI_Isendrecv(&array[0], 1, MPI_DOUBLE, myrank - 1, myrank - 1, &tmp1, 1, MPI_DOUBLE, myrank - 1, myrank - 1, MPI_COMM_WORLD, &Status);
                                MPI_Isendrecv(&array[size/h - 1], 1, MPI_DOUBLE, myrank + 1, myrank + 1, &tmp2, 1, MPI_DOUBLE, myrank + 1, myrank + 1, MPI_COMM_WORLD, &Status);
				newarray[0] = array[0] + k*t/(h * h) * (array[1] - 2 * array[0] + tmp1);
				newarray[size/h - 1] = array[size/h - 1] + k*t/(h * h) * (tmp2 - 2 * array[size/h - 1] + array[size/h - 2]);
			}
			else
			{
				tmp2 = 0;
                                MPI_Isendrecv(&array[0], 1, MPI_DOUBLE, myrank - 1, myrank - 1, &tmp1, 1, MPI_DOUBLE, myrank - 1, myrank - 1, MPI_COMM_WORLD, &Status);
                                newarray[0] = array[0] + k*t/(h * h) * (array[1] - 2 * array[0] + tmp1);
                                newarray[size - 1] = array[size/h - 1] + k*t/(h * h) * (tmp2 - 2 * array[size/h - 1] + array[size/h - 2]);

			}
		}
		for(i = 1; i < size/h - 2; i++)
		{
			newarray[i] = array[i] + k*t/(h * h) * (array[i + 1] - 2 * array[i] + array[i - 1]);
		}
		memcpy(array, newarray, size/h);
	}
}
