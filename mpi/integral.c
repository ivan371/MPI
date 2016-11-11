#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>

int main(int argc, char* argv[])
{
	int i;
	int array[10];
	int myrank, size;
	MPI_Status Status;
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
	if(myrank == 0)
	{
		for(i = 1; i < size; i++)
		{
			MPI_Send(&buf[i*N/size], N/size, MPI_INT, i, i, MPI_COMM_WORLD);
		}        }
        else
        {
                MPI_Recv (&buf[0], N/size, MPI_INT, 0, myrank, MPI_COMM_WORLD, &Status);
        }
        printf("I am %d of %d\n, myrank, size");
	MPI_Finalize();
        return 0;
}


