#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <mpi.h>
#define h 0.00005
#define k 1
#define T 0.00001
#define L 1
#define EXCHANGE_TAG 33
#define Yard 1000
#define PI 3.14159
#define E 2.71828

double once();
void precive(double x);
int main(int argc, char* argv[])
{
	int myrank, size, begin, end;
 	MPI_Status Status;
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
	MPI_Request request;
	double t0, t1;
	t0 = MPI_Wtick();
	t1 = MPI_Wtime();
	int N = L/h + 1;
	double   total, tmp1, tmp2;
	int i = 0;
	int j = 0;
	int s = 1 / (h * size);
	//double *array[2];
	double* array;
	double* newarray;
	array = (double*)calloc(2 * N, sizeof(double));
	newarray = (double*)calloc(2 * N, sizeof(double));
	double t = h * h / k;
	int steps = T / t;
	int num = N / size;
	begin = num * myrank;
	end = begin + num;
	if (myrank == 0)
	{
		 begin = 1;
		 t1 = MPI_Wtime();
	}
	if(myrank == size - 1) end = N - 1;
	for (i = begin - 1; i <= end; i++)
        {
                array[i] = 1;
                newarray[i] = 1;
        }
	array[0] = 0;
	newarray[0] = 0;
	array[N - 1] = 0;
	newarray[N - 1] = 0;
	//end = N - 1;
	//end = end + 1;	
	int tm = 0;
	double tou = 0.005;
	double co = (tou * k) / (h * h); 
	for(i = 0; i < steps; i++)
	{
		if (i % 2 == 0)
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
			newarray[j] = array[j] + co * (array[j + 1] - 2.0 * array[j] + array[j - 1]);
		  //	printf("%lf ", newarray[j]);
		}
		}
           	else
                {
                if(myrank == 0)
                {
                        MPI_Isend(&newarray[end - 1], 1, MPI_DOUBLE, myrank + 1, EXCHANGE_TAG, MPI_COMM_WORLD, &request);
                        MPI_Irecv(&newarray[end], 1, MPI_DOUBLE, myrank + 1, EXCHANGE_TAG, MPI_COMM_WORLD, &request);
                        MPI_Wait(&request, &Status);
                }
                else
                {
                        if(myrank != size - 1)
                        {
                                MPI_Isend(&newarray[begin], 1, MPI_DOUBLE, myrank - 1, EXCHANGE_TAG, MPI_COMM_WORLD, &request);
                                MPI_Irecv(&newarray[begin - 1], 1, MPI_DOUBLE, myrank - 1, EXCHANGE_TAG, MPI_COMM_WORLD, &request);
                                MPI_Isend(&newarray[end - 1], 1, MPI_DOUBLE, myrank + 1, EXCHANGE_TAG, MPI_COMM_WORLD, &request);
                                MPI_Irecv(&newarray[end], 1, MPI_DOUBLE, myrank + 1, EXCHANGE_TAG, MPI_COMM_WORLD, &request);
                                MPI_Wait(&request, &Status);
                        }
                        else
                        {
                                MPI_Isend(&newarray[begin], 1, MPI_DOUBLE, myrank - 1, EXCHANGE_TAG, MPI_COMM_WORLD, &request);
                                MPI_Irecv(&newarray[begin - 1], 1, MPI_DOUBLE, myrank - 1, EXCHANGE_TAG, MPI_COMM_WORLD, &request);
                                MPI_Wait(&request, &Status);
                        }
                }
                for(j = begin; j < end; j++)
                {
                        array[j] = newarray[j] + co * (newarray[j + 1] - 2.0 * newarray[j] + newarray[j - 1]);
                    //    printf("%lf ", array[j]);
                }
                }

		//printf(" %d\n", myrank);
	}
	/*if (i % 2 == 0)
	{
	    memcpy(array, newarray, N);
	}*/
	//printf("\n");
	if (myrank == 0) {
		for (i = 1; i < size - 1; i++) {			
			MPI_Recv(&array[num * i], num, MPI_DOUBLE, i, 0, MPI_COMM_WORLD, &Status);
		}
		MPI_Recv(&array[num * (size - 1)], num  +  N % size, MPI_DOUBLE, i, 0, MPI_COMM_WORLD, &Status);
		//for (i = 0; i < N; i++) {
		//	printf("%f %f\n", h * i, array[i]);
		//}
		t1 = MPI_Wtime() - t1;
		printf("%lf \n", once() / t1);
	} else {
		if(myrank == size - 1)
		    MPI_Send(&array[num * myrank], num + N % size, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
		else
		    MPI_Send(&array[num * myrank], num, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
	}
	//printf("\n");
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
double once()
{
        int i = 0;
        int j = 0;
	double t0, t1;
	t0 = MPI_Wtick();
	t1 = MPI_Wtime();
        double* array;
        double* newarray;
        int N = L/h + 1;
        array = (double*)calloc(2 * N, sizeof(double));
        newarray = (double*)calloc(2 * N, sizeof(double));
        double t = h * h / k;
        int steps = T / t;
        double tou = 0.005;
        double co = (tou * k) / (h * h);
        for(i = 1; i < N - 1; i++)
        {
                array[i] = 1;
                newarray[i] = 1;
        }
	for(i = 0; i < steps; i++)
        {
                if (i % 2 == 0)
                {
                        for(j = 1; j < N - 1; j++)
                        {
                                newarray[j] = array[j] + co * (array[j + 1] - 2.0 * array[j] + array[j - 1]);
                        }
                }
                else
                {
                        for(j = 1; j < N - 1; j++)
                        {
                                array[j] = newarray[j] + co * (newarray[j + 1] - 2.0 * newarray[j] + newarray[j - 1]);
                        }
                }
        }
	if (i % 2 == 0)
        {
                memcpy(array, newarray, N);
        }
	for (i = 0; i < N; i++)
        {
 //               printf("%lf %lf\n",h * i, array[i]);
        }
	free(array);
        free(newarray);
	return MPI_Wtime() - t1;
}

