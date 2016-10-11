#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>

#define PI 3.141592653
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
    int it;
    float *x, *y;
    double pi;
    pi=0.0;

    // Allocate memory for our random numbers
    x = (float*)malloc(N*sizeof(float));
    y = (float*)malloc(N*sizeof(float));


    for(int i=0; i<10024;i++){
        generate_random_numbers(N, x, y);
        pi += calculate_pi_monte_carlo(N,x,y);
        it = i+1;
        if (i%50==1){
            printf("Samples : %d, ", (N*it));
            printf("Pi Estimated : %f, ", pi/it);
            printf("Error : %f\n", (PI-(pi/it))/PI);
        } 
    }
    pi /= 1024;
    printf("Pi is estimated to be %f\n", pi);

    // Convert timing to a readable format (ms)
    //random_time = rt * 1000 / CLOCKS_PER_SEC;
    //monte_carlo_time = mct * 1000 / CLOCKS_PER_SEC;
    // Free memory
    free(x);
    free(y);
}