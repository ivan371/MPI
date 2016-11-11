#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#define split 100

int main(int argc, char* argv[])
{
  int i;
  double sum = 0;
  double val1 = 0;
  double val2 = 0;
  int myrank, size;
  MPI_Status Status;
  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
  if(myrank == 0)
  {
    int* array = (int*)calloc(2 * size, sizeof(int));
    for(i = 1; i < size - 1; i++)
    {
	array[2 * i] = (i * split)/size;      //begin
	array[2 * i + 1] = (i + 1) * split/size;            //end
    }
    array[2 * size - 2] = (size - 1)*split/size;
    array[2 * size - 1] = split;
    for (i = 1; i < size; i++)
    {
	MPI_Send(&array[2 * i], 2, MPI_INT, i, i, MPI_COMM_WORLD);
    }
    for(i = 0; i < array[2]; i++)
      {
	  val1 = 1/(1 + ((double)i / split) * ((double)i / split));
	  val2 = 1/(1 + ((double)(i + 1) / split) * ((double)(i + 1) / split));
	  sum += (val1 + val2) / (2 * split);
      }
      printf("I am %d begin %d end %d sum %lf\n", myrank, 0, array[2], sum);
  }
  else
  {
      int array[2];
      MPI_Recv(array, 2, MPI_INT, 0, myrank, MPI_COMM_WORLD, &Status);
      for(i = array[0]; i < array[1]; i++)
      {
	  val1 = 1/(1 + ((double)i / split) * ((double)i / split));
	  val2 = 1/(1 + ((double)(i + 1) / split) * ((double)(i + 1) / split));
	  sum += (val1 + val2) / (2 * split);
      }
      printf("I am %d begin %d end %d sum %lf\n", myrank, array[0], array[1], sum);
  }
  MPI_Finalize();
  return 0;
  
}
