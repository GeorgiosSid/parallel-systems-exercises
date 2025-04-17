// SIDIROPOULOS GEORGIOS

#include <stdio.h>
#include <sys/time.h>
#include <stdlib.h>
#include <omp.h>

#define UPTO 10000000

long int count,      /* number of primes */
         lastprime;  /* the last prime found */


// timeit function accept a pointer to a function taking a long int and returning void
static
double timeit(void (*func)(long int), long int arg) {
    struct timeval start, end;
    gettimeofday(&start, NULL);
    func(arg);
    gettimeofday(&end, NULL);
    return (double) (end.tv_usec - start.tv_usec) / 1000000
        + (double) (end.tv_sec - start.tv_sec);
}

void serial_primes(long int n) {
	long int i, num, divisor, quotient, remainder;

	if (n < 2) return;
	count = 1;                         /* 2 is the first prime */
	lastprime = 2;

	for (i = 0; i < (n-1)/2; ++i) {    /* For every odd number */
		num = 2*i + 3;

		divisor = 1;
		do
		{
			divisor += 2;                  /* Divide by the next odd */
			quotient  = num / divisor;
			remainder = num % divisor;
		} while (remainder && divisor <= quotient);  /* Don't go past sqrt */

		if (remainder || divisor == num) /* num is prime */
		{
			count++;
			lastprime = num;
		}
	}
}


void openmp_primes(long int n) {
	long int i, num, divisor, quotient, remainder;

	if (n < 2) return;
	count = 1;                         /* 2 is the first prime */
	lastprime = 2;

    omp_set_dynamic(0);                /* nothing to do with the dynamic schedule! */
    omp_set_num_threads(4);            //schedule(runtime): for the different types of parallelism 
    #pragma omp parallel for  schedule(dynamic, 100000) private(i, num, divisor, quotient, remainder) reduction(+:count) reduction(max:lastprime) shared(n) 
        for (i = 0; i < (n-1)/2; ++i) {    /* For every odd number */
            num = 2*i + 3;

            divisor = 1;
            do
            {
                divisor += 2;                  /* Divide by the next odd */
                quotient  = num / divisor;
                remainder = num % divisor;
            } while (remainder && divisor <= quotient);  /* Don't go past sqrt */

            if (remainder || divisor == num) /* num is prime */
            {
                count++;
                lastprime = num;
            }
        }
}


int main()
{
    double exectime, exectimepar;

	printf("Serial and parallel prime number calculations:\n\n");


    /* serial_primes timing */
    exectime = timeit(serial_primes, UPTO);
	printf("[serial] count = %ld, last = %ld (time = %f)\n ------------------- \n", count, lastprime, exectime);

    /* openmp_primes timing */
    exectimepar = timeit(openmp_primes, UPTO);
	printf("[openmp] count = %ld, last = %ld (time = %f)\n", count, lastprime, exectimepar);

    //-- for the different types of parallelism --
    /* 
    int numThreads[] = {2, 3, 4};
    omp_sched_t schedules[] = {omp_sched_static, omp_sched_dynamic, omp_sched_guided};
    char* scheduleNames[] = {"static", "dynamic", "guided"};
    int chunkSizes[] = {100000, 250000, 500000, 1000000, 2500000}; // Example chunk sizes
    double sumExectime = 0;
    int runs = 4;

    for (int threadIndex = 0; threadIndex < 3; threadIndex++){
        for (int schedIndex = 0; schedIndex < 3; schedIndex++){
            for (int chunkIndex = 0; chunkIndex < 5; chunkIndex++){
                printf("Running with %d threads, schedule %s, chunk size %d four times.\n",
                        numThreads[threadIndex], scheduleNames[schedIndex], chunkSizes[chunkIndex]);
                for (int run = 0; run < runs; run++){
                    omp_set_dynamic(0);
                    omp_set_num_threads(numThreads[threadIndex]);
                    omp_set_schedule(schedules[schedIndex], chunkSizes[chunkIndex]);

                    count = 0;
                    lastprime = 0;

                    sumExectime += timeit(openmp_primes, UPTO);
                }
                printf("Average time is = %f\n\n", sumExectime / runs);
                sumExectime = 0;
            }
        }
    }
    */

	return 0;
}
