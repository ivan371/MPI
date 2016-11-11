#include<stdio.h>
#include<stdlib.h>
#include<pthread.h>
#include<math.h>
#include <string.h>
#include <time.h>
#include <semaphore.h>

typedef struct thread_data_struct
{
    long points;
}pthread_data;

void* pthread_func( void* arg_struct);

double calculate(long points);


double pi = 3.1415926;
sem_t sem;
double common_res;

int main(int argc, char *argv[])
{
    sem_init(&sem, 0, 1);
    struct timespec begin, end;
    double elapsed;
    clock_gettime(CLOCK_REALTIME, &begin); 

    int numb = atoi(argv[1]); //number of processes
    long points = atol(argv[2]);

    //long points = argv[2];
    //long numb = argv[1];
    pthread_t* pthr = (pthread_t*) malloc(numb*sizeof(pthread_t));
    pthread_data* data = (pthread_data*) malloc(numb*sizeof(pthread_data));
    long points_for_proc = (long) points/numb;
    int i = 0, rc;
    for( i = 0; i < numb; i++)
    {
        data[i].points = points_for_proc;
        rc = pthread_create(&pthr[i], NULL, (pthread_func), (void*) &data[i]);
        if (rc) printf("ERROR; return code from pthread_create() is %d \n", rc);
    }

    for(i = 0; i < numb; i++)
    {
        rc = pthread_join(pthr[i], NULL);
        if (rc) printf("ERROR; return code from pthread_join() is %d \n", rc);
    }

    double result = common_res/numb;
    printf("Result is %.4f\n", result);

    clock_gettime(CLOCK_REALTIME, &end);
    elapsed = end.tv_sec - begin.tv_sec;
    elapsed += (end.tv_nsec - begin.tv_nsec) / 1000000000.0;
    printf("Time is %.4f\n", elapsed);

    free(data);
    free(pthr);
    sem_destroy(&sem); 


    return 0;
}

void* pthread_func( void* arg_struct)
{
    pthread_data* data = (pthread_data*) arg_struct;
    long points = data->points;
    double res = calculate(points);

    sem_wait(&sem);
    common_res += res;
    sem_post(&sem);

    pthread_exit(NULL);
}

double calculate(long points)
{
    double res = 0.0;
    long i = 0;
    int x_k;
    double x, y;
    long goals = 0;
    for(i = 0; i< points; i ++ )
    {
        x = ((double)rand_r(&x_k) / RAND_MAX) * pi;
        y = (double)rand_r(&x_k) / RAND_MAX;
	//printf("%.2f %.2f\n", x, y);
        if (y <= sin(x))
        {
            goals ++;
            res += x*y;
        }
    }
    return 2*res/goals;
}

