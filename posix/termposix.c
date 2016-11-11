#include <stdio.h>
#include <stdlib.h>
#include <pthread.h>
#include <semaphore.h>
#include <math.h>
#include <time.h>

#define c 1
#define h 0.1
#define tou 0.05
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
int main(int argc, char* argv[])
{
        if (argc < 2)
        {
                printf("NOT 3 argc");
                exit(1);
        }
	NUM_THREADS = atoi(argv[1]);
	int N = MAX_N / h;
	array = (double*)calloc(N, sizeof(double));
	new_array = (double*)calloc(N, sizeof(double));
	int param;
        int rc, i;
        void *arg;
      	sem = (sem_t*)calloc(NUM_THREADS, sizeof(sem_t));
	for (i = 0; i < NUM_THREADS; i++)
	{
		sem_init(&sem[i], 0, 1);
	}
	for(i = 0; i < N; i++)
	{
		if(h * i >= 0 || h * i <= 2 )
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
        pthread_t* pthr = calloc(NUM_THREADS, sizeof(pthread_t));
        for (i = 0; i < NUM_THREADS; i++)
        {
                rc = pthread_create(&pthr[i], NULL, start_func, i);
        }
	for (i = 0; i < NUM_THREADS; i++)
        {
                pthread_join(pthr[i], &arg);
        }
        free(pthr);
	for(i = 0; i < NUM_THREADS; i++)
	{
		sem_destroy(&sem);
	}
        return 0;
}
void* start_func(int param)
{
        int val;
        int local;
        local = (int*)param;
        int i, j;
	int steps = T / tou;
	for (i = 0; i < steps; i++)
	{
		if(local == NUM_THREADS - 1)
		{
			sem_wait(&sem[0]);
		}
		else
		{
			sem_wait(&sem[local + 1]);
		}
		for(j = 1; j <= 1 / h; j++)
                {
                        if(i % 2 == 0)
                        {
                            new_array[j] = array[j] - c * tou * (array[j] - array[j - 1]);
                        }
                        else
                        {
                            array[j] = new_array[j] - c * tou  * (new_array[j] - new_array[j - 1]);
                        }
                }		
		if(local == NUM_THREADS - 1)
		{
			sem_post(&sem[NUM_THREADS - 1]);
		}
		else
		{
        		sem_post(&sem[local]);
		}
	}
        pthread_exit(NULL);
}

