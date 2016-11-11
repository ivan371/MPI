#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <math.h>
//#define split 100

int main(int argc, char* argv[])
{
  int i;
  int split = 0;
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
    for(split = 2; split < 150; split++)
    {
      sum = 0;
      for(i = 0; i < split; i++)
      {
	  val1 = 1/(1 + ((double)i / split) * ((double)i / split));
	  val2 = 1/(1 + ((double)(i + 1) / split) * ((double)(i + 1) / split));
	  sum += (val1 + val2) / (2 * split);
      }
      printf(" %d %lf %lf\n", split, sum, ((0.78539 - sum)/0.78539));
    }
  }
  MPI_Finalize();
  return 0;
  
}
