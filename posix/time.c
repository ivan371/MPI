#include <stdio.h>
#include <stdlib.h>
#include <pthread.h>
#include <semaphore.h>
#include <math.h>
#include <time.h>

#define PI 3.14159
int MAX_TRANS = 10000000;
int TR;

double total_sum = 0;
sem_t sem;
double tm2;

void *start_func(int param);
double once();
int newrand(int* nextp);

int main(int argc, char* argv[])
{
	if (argc < 2)
	{
		printf("NOT 3 argc");
		exit(1);
	}

	int param;
        int rc, i;
        void *arg;
	struct timespec begin, end, begin1, end1, end2, begin2;
	double elapsed, tme;
	pthread_t pthrd;
	sem_init(&sem, 0, 1);
	clock_gettime(CLOCK_REALTIME, &begin1);
	once();
	//rc = pthread_create(&pthrd, NULL, start_func, (void*)&param);
	//pthread_join(pthrd, &arg);
	clock_gettime(CLOCK_REALTIME, &end1);
	tme = end1.tv_sec - begin1.tv_sec;
	tme += (end1.tv_nsec - begin1.tv_nsec)/1000000000.0;
	total_sum  = 0;
	clock_gettime(CLOCK_REALTIME, &begin);

	int NUM_THREADS = atoi(argv[1]);
	MAX_TRANS = MAX_TRANS / NUM_THREADS;
	pthread_t* pthr = calloc(NUM_THREADS, sizeof(pthread_t));
	for (i = 0; i < NUM_THREADS; i++)
	{
		rc = pthread_create(&pthr[i], NULL, start_func, i);
	}
	for (i = 0; i < NUM_THREADS; i++)
	{
		pthread_join(pthr[i], &arg);
	}

	clock_gettime(CLOCK_REALTIME, &end);
	elapsed = end.tv_sec - begin.tv_sec;
	elapsed += (end.tv_nsec - begin.tv_nsec)/1000000000.0;
	printf("%lf %lf\n", total_sum /(NUM_THREADS), tme/elapsed);
	free(pthr);
	sem_destroy(&sem);
	return 0;
}

void* start_func(int param)
{	
	int val;
	int* local;
	local = (int*)param;
	int i;
	double sum = 0;
	int all_sum = 0;
	double x = 0;
	double y = 0;
	int x_k;
	for (i = 0; i < MAX_TRANS; i++)
	{
		x = ((double)rand_r(&x_k) / RAND_MAX) * PI;
		y = ((double)rand_r(&x_k) / RAND_MAX);
		if(y <= sin(x))
		{
			all_sum++; 
			sum += x * y;; 
		}
	}
	sem_wait(&sem);
	//printf("%lf %lf\n", all_sum, sum);
	total_sum += 2 * (sum / all_sum);
//	printf("%d\n", total_sum);
	sem_post(&sem);
//	sem_getvalue(&sem, &val);
	pthread_exit(NULL); 
}

double once()
{
	int i;
	double sum = 0;
	int all_sum = 0;
	double x = 0;
	double y = 0;
	int x_k;
	for (i = 0; i < MAX_TRANS; i++)
	{
		x = ((double)rand_r(&x_k) / RAND_MAX) * PI;
		y = ((double)rand_r(&x_k) / RAND_MAX);
		if(y <= sin(x))
                {
                        all_sum++;
                        sum += x * y;;
                }

	}
	return sum;
}


int newrand(int* nextp)
{
	*nextp = *nextp * 1103515245 + 12345;
	return (*nextp / 65536) % 32768;
}
