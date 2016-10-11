#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>


/* Function : generate_random_numbers
 * Generates n random numbers for both x and y on the host
 */
void generate_random_numbers(int n, float*x, float*y){
	srand(time(NULL));
	for (int i = 0; i < n; i++) {
		x[i] = ((float)rand() / RAND_MAX) - 0.5f;
		y[i] = ((float)rand() / RAND_MAX) - 0.5f;
	}
}

/* Function : calculate_pi_monte_carlo
 * Calculates pi on the host by using the monte carlo method, by using a set of
 * random points within a 2R square about point (0,0) we can calculate pi by
 * calculating the ratio of points within the a circle with radius R starting
 * from point (0,0) compared to that of the square. This is done on the host.
 */
float calculate_pi_monte_carlo(int n, float*x, float*y){
	int circle = 0;
	for(int i = 0; i < n; i++) {
		if( pow(x[i],2) + pow(y[i],2)  < pow(0.5f,2) ){
			circle++;
		}
	}
	return (4.0f*circle)/n;
}

int main(void)
{
    // Initiate variables
    int N = 1048576;
    int random_time, monte_carlo_time;
    float *x, *y;
    float pi;

    // Allocate memory for our random numbers
    x = (float*)malloc(N*sizeof(float));
    y = (float*)malloc(N*sizeof(float));

    //clock_t start = clock(), rt, mct; // Start Timer

    generate_random_numbers(N, x, y);

    //rt = clock() - start; // Time Random generation

    pi = calculate_pi_monte_carlo(N, x, y);

    //mct = clock() - rt; // Time Monte Carlo

    // Convert timing to a readable format (ms)
    //random_time = rt * 1000 / CLOCKS_PER_SEC;
    //monte_carlo_time = mct * 1000 / CLOCKS_PER_SEC;

    // Print out results
    printf("Random Time taken %d seconds %d milliseconds\n", random_time/1000, random_time%1000);
    printf("Monte Carlo Time taken %d seconds %d milliseconds\n", monte_carlo_time/1000, monte_carlo_time%1000);
    printf("Iterations : %d\n",N);
    printf("Pi is estimated to be %f\n", pi);

    // Free memory
    free(x);
    free(y);
}