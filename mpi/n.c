#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mpi.h>

#define tay 0.005
#define h 0.1
#define k 1

double function(double t, double x){
	int m=0;
	double v=0, u=0;
	do {
	u += 4/M_PI*(exp(-k*(M_PI)*(M_PI)*(2*m+1)*(2*m+1)*t))/(2*m+1)*sin(M_PI*(2*m+1)*x);
	m++;
	v= u+ 4/M_PI*(exp(-k*(M_PI)*(M_PI)*(2*m+1)*(2*m+1)*t))/(2*m+1)*sin(M_PI*(2*m+1)*x);
	} while(fabs(u-v)>0.00001);
	return u;
	}

int main(int argc, char *argv[]){
int myrank, size;
MPI_Status Status;
MPI_Init(&argc, &argv);
MPI_Comm_size(MPI_COMM_WORLD, &size);
MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

int I = 1/h;
int N = 1/tay;
int ni = I/size;
int r = I+1-size*ni;

int T = 0.1/tay;
if (r>myrank){
	r=1;
	}
else r=0;
//double u [T+2][ni+2+r];
int j=1;
int i=1;
double *u[T+2];
for (i = 0; i < T+2; ++i)
{
u[i] = (double *)calloc(ni+2+r, sizeof(double));
}

for (j=1; j<ni+2+r; j++){
	u[0][j]=1;
	}
if (myrank == 0){
	for (j=0; j<T+1; j++){
	u[j][1]=0;
	}
	}
if (myrank==size-1){
	for (i=0; i<T+1; i++){
	u[i][ni+r]=0;
	}
	}
double x=0;
printf("I am %d, ni+r = %d+%d, T is %d\n", myrank, ni, r, T);

if(size==1){
	printf("The number of processes should be more than one\n");
	exit(-1);
	}
int begin = 0;
int finish = 0;
for (j=0; j<T+1; j++){
	x=myrank*ni*h;
	if (myrank == 0){
		MPI_Send(&u[j][ni+r], 1, MPI_DOUBLE, myrank+1, myrank+1, MPI_COMM_WORLD);
		MPI_Recv(&u[j][ni+1+r], 1, MPI_DOUBLE, myrank+1, myrank, MPI_COMM_WORLD, &Status);
		begin = 2;
		finish = ni+1+r;
		}
	if (myrank!=size-1 && myrank!=0){
		MPI_Send(&u[j][ni+r], 1, MPI_DOUBLE, myrank+1, myrank+1, MPI_COMM_WORLD);
		MPI_Send(&u[j][1], 1, MPI_DOUBLE, myrank-1, myrank-1, MPI_COMM_WORLD);
		MPI_Recv(&u[j][0], 1, MPI_DOUBLE, myrank-1, myrank, MPI_COMM_WORLD, &Status);
		MPI_Recv(&u[j][ni+1+r], 1, MPI_DOUBLE, myrank+1, myrank, MPI_COMM_WORLD, &Status);
		begin = 1;
		finish = ni+1+r;
		}
	if (myrank==size-1){
		MPI_Send(&u[j][1], 1, MPI_DOUBLE, myrank-1, myrank-1, MPI_COMM_WORLD);
		MPI_Recv(&u[j][0], 1, MPI_DOUBLE, myrank-1, myrank, MPI_COMM_WORLD, &Status);
		begin = 1;
		finish = ni+r;
		}
	for (i=begin; i<finish; i++){
		//printf("%d j is %d x %lf T %lf %d\n",myrank, j, x+h, u[j][i], T);
		u[j+1][i] = u[j][i]+k*tay/(h*h)*(u[j][i+1]-2*u[j][i]+u[j][i-1]);
		if (j==T){
			printf("x %lf T %lf\n", x+h, u[j][i]);
			} 
		x+=h;
		}
		
	}
if (myrank == 0){
double y[10];
x=h;
	for(j=0; j<9;j++){
		y[j] = function(0.1, x);
		x+=h;
		printf("%lf\n", y[j]);
	}
}
MPI_Finalize();

return 0;
}


