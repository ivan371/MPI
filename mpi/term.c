#include <stdio.h>
#include <math.h>
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


}
