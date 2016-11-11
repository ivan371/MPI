#include <stdio.h>
#include <stdlib.h>
#include <pthread.h>
#include <semaphore.h>
#include <math.h>

#define PI 3.14159
#define MAX_TRANS 1000000

double total_sum = 0;
sem_t sem;

void *start_func(void* param);
int main(int argc, char* argv[])
{
	if (argc < 2)
	{
		printf("NOT 3 argc");
		exit(1);
	}
	int NUM_THREADS = atoi(argv[1]);
	int param;
	int rc, i;
	void *arg;
	sem_init(&sem, 0, 1);
	pthread_t* pthr = calloc(NUM_THREADS, sizeof(pthread_t));
	for (i = 0; i < NUM_THREADS; i++)
	{
		sem_init(&sem, 0, 1);
		rc = pthread_create(&pthr[i], NULL, start_func, (void*)&param);
		pthread_join(pthr[i], &arg);
		sem_destroy(&sem);
	}
	printf("%lf\n", total_sum /(NUM_THREADS * MAX_TRANS));
	return 0;
}

void* start_func(void* param)
{	
	int val;
	sem_wait(&sem);
	int* local;
	local = (int*)param;
	int i;
	double sum = 0;
	double x = 0;
	double y = 0;
	for (i = 0; i < MAX_TRANS; i++)
	{
		x = (rand() % 4000)/1000;
		y = (rand() % 1000)/1000;
		if((x < PI) && (y < sin(x)))
		{
			sum++; 
		}
	}
	total_sum += sum;
//	printf("%d\n", total_sum);
	sem_post(&sem);
	sem_getvalue(&sem, &val);
	pthread_exit(NULL); 
}
