#include <stdio.h>
#include <stdlib.h>
#include <pthread.h>
#include <semaphore.h>
#include <math.h>
#include <time.h>

#define c 1
//#define h 0.00001
#define h 0.1
#define tou 0.1
#define T 4
#define PI 3.14159
#define MAX_N 10
int MAX_TRANS = 10000000;
int TR;
int NUM_THREADS = 0;

sem_t* sem;
double* array;
double* new_array;

void *start_func(int param);
void* print_value();
int main(int argc, char* argv[])
{
    if (argc < 2)
    {
        printf("NOT 3 argc");
        exit(1);
    }
	NUM_THREADS = atoi(argv[1]);
	struct timespec begin, end;
	double elapsed, elapsed1;
	clock_gettime(CLOCK_REALTIME, &begin);
	once();
	clock_gettime(CLOCK_REALTIME, &end);
	elapsed1 = end.tv_sec - begin.tv_sec;
   	elapsed1 += (end.tv_nsec - begin.tv_nsec) / 1000000000.0;

    clock_gettime(CLOCK_REALTIME, &begin);

    int N = MAX_N / h;
    array = (double*)calloc(N, sizeof(double));
    new_array = (double*)calloc(N, sizeof(double));
	double* test_array = (double*)calloc(N, sizeof(double));    
	int param;
    int rc, i;
    void *arg;
    sem = (sem_t*)calloc(NUM_THREADS, sizeof(sem_t));
    for (i = 0; i < NUM_THREADS; i++)
    {
        sem_init(&sem[i], i, 0);
    }
    sem_post(&sem[0]);
	for(i = 0; i < N; i++)
    	{
    		if(h * i >= 0 && h * i <= 2 )
    		{
    	    		array[i] = h * i * (2 - h * i);
    	    		new_array[i] = h * i * (2 - h * i);
    		}
        	else
        	{
            		array[i] = 0;
           		 new_array[i] = 0;
        	}
    	}
	for(i = 0; i < N; i++)
	{
		if(h * i >= 4 && h * i <= 6 )
    	{
    	    test_array[i] = (h * i - 4) * (6 - h * i);
    	}
        else
        {
            test_array[i] = 0;
        }
	}
	pthread_t* pthr = calloc(NUM_THREADS, sizeof(pthread_t));
	if(NUM_THREADS == 1)
		once();
	else
	{
    		for (i = 0; i < NUM_THREADS; i++)
    		{
        		rc = pthread_create(&pthr[i], NULL, start_func, i);
    		}
		for (i = 0; i < NUM_THREADS; i++)
    		{
        		pthread_join(pthr[i], &arg);
    		}
	}
	free(pthr);
	for(i = 0; i < N; i++)
	{
		printf("%lf %lf\n", new_array[i], test_array[i]);
	}
	clock_gettime(CLOCK_REALTIME, &end);
    elapsed = end.tv_sec - begin.tv_sec;
    elapsed += (end.tv_nsec - begin.tv_nsec) / 1000000000.0;
//    printf("%lf\n", elapsed1/elapsed);
    
	for(i = 0; i < NUM_THREADS; i++)
    {
            sem_destroy(&sem[i]);
    }
	return 0;
}
void* print_value(local)
{
	int i = 0;
	int value;
	for(i = 0; i < NUM_THREADS; i++)
	{
		sem_getvalue(&sem[i], &value);
		printf("%d ", value);
	}
	printf("%d \n", local);
}
void* start_func(int param)
{
        int val;
        int local;
        local = (int*)param;
        int i, j;
        int steps = T / tou;
	int begin = MAX_N  * local / (h * NUM_THREADS);
	int end =  MAX_N * (local + 1) / (h * NUM_THREADS);
        for (i = 0; i < steps; i++)
        {
                sem_wait(&sem[local]);
               	for(j = begin; j <= end; j++)
                {
                        if(i % 2 == 0)
                       	{
                            new_array[j] = array[j] - c * tou * (array[j] - array[j - 1]) / h;
			           	}
                       	else
                       	{
                            array[j] = new_array[j] - c * tou  * (new_array[j] - new_array[j - 1]) / h;
                       	}
               	}
                if(local != NUM_THREADS - 1)
                {                
			sem_post(&sem[local + 1]);
		}
		if(local != 0)
			sem_post(&sem[local - 1]);
        }
        pthread_exit(NULL);
}

once()
{
	int i, j;
	int steps = T / tou;
	 int N = MAX_N / h;	
	 double* array1 = (double*)calloc(N, sizeof(double));
	 double* array2 = (double*)calloc(N, sizeof(double));
	 for (i = 0; i < steps; i++)
        {
                for(j = 1; j <= N; j++)
                {
                        if(i % 2 == 0)
                        {
                            array1[j] = array2[j] - c * tou * (array2[j] - array2[j - 1]) / h;
                        }
                        else
                        {
                            array2[j] = array1[j] - c * tou  * (array1[j] - array1[j - 1]) / h;
                        }
                }
        }

}
