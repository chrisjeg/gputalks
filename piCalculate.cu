// nvcc -ccbin "D:\Program Files (x86)\Microsoft Visual Studio 11.0\VC\bin" piCalculate.cu -o piCalculate.exe

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <cuda.h>
#include <curand.h>
#include <curand_kernel.h>
#include <cuda_profiler_api.h>

#define MAX_CUDA_BLOCKS 65535
#define MAX_CUDA_THREADS 1024
#define PI 3.141592653

__global__ void cuda_calc_pi_step1(int n, int *circle, float *x, float *y)
{
    extern __shared__ int sdata[];
    float c;
    int t = 0;

    // each thread loads one element from global to shared mem
    unsigned int tid = threadIdx.x;
    unsigned int i = blockIdx.x * blockDim.x + threadIdx.x;

    c = (x[i] * x[i]) + (y[i] * y[i]);
    if (c < 0.25f && i<n) t = 1;
    sdata[tid] = t;
    __syncthreads();
    for (unsigned int s=blockDim.x/2; s>0; s>>=1) {
    if (tid < s) {
    sdata[tid] += sdata[tid + s];
    }
    __syncthreads();
    }
    // write result for this block to global mem
    if (tid == 0){
        circle[blockIdx.x] = sdata[0];
    };
}

__global__ void cuda_calc_pi_step2(int n, int *circle){
    extern __shared__ int sdata[];
    unsigned int tid = threadIdx.x;
    sdata[tid] = circle[tid];

    __syncthreads();

    // do reduction in shared mem
    for (unsigned int s = 1; s < blockDim.x; s *= 2){
        if (tid % (2 * s) == 0){
            sdata[tid] += sdata[tid + s];
        }
        __syncthreads();
    }

    // write result for this block to global mem
    if (tid == 0){
        circle[blockIdx.x] = sdata[0];
    };
}

/* Function : generate_random_numbers
 * Generates n random numbers for both x and y on the host
 */
void generate_random_numbers(int n, float *x, float *y)
{
    srand(time(NULL));
    for (int i = 0; i < n; i++)
    {
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
float calculate_pi_monte_carlo(int n, float *x, float *y)
{
    int circle = 0;
    for (int i = 0; i < n; i++)
    {
        if (pow(x[i], 2) + pow(y[i], 2) < pow(0.5f, 2))
        {
            circle++;
        }
    }
    return (4.0f * circle) / n;
}

double gpu_calc_pi_monte_carlo(int samples){
    int threads, blocks;

    //You can only do so many samples in an execution due to limitations of the card
    if(samples > (MAX_CUDA_THREADS*MAX_CUDA_THREADS)){
        printf("Too many samples\n");
        return 0.0f;
    }
    threads = (samples < MAX_CUDA_THREADS) ? samples : MAX_CUDA_THREADS;
    blocks = ((samples-1)/MAX_CUDA_THREADS)+1;

    int *d_circle_count, *circle_count;
    float *x, *y, *d_x, *d_y;
    double pi;

    // Allocate memory for our random numbers
    x = (float *)malloc(samples * sizeof(float));
    y = (float *)malloc(samples * sizeof(float));
    circle_count = (int *)malloc(blocks * sizeof(int));

    // Allocate memory on our GPU
    cudaMalloc(&d_x, samples * sizeof(float));
    cudaMalloc(&d_y, samples * sizeof(float));
    cudaMalloc(&d_circle_count,blocks * sizeof(int));

    generate_random_numbers(samples, x, y);

    cudaMemcpy(d_x, x, samples * sizeof(float), cudaMemcpyHostToDevice);
    cudaMemcpy(d_y, y, samples * sizeof(float), cudaMemcpyHostToDevice);

    cuda_calc_pi_step1<<<blocks,threads,threads*sizeof(int)>>>(samples, d_circle_count, d_x, d_y);
    cuda_calc_pi_step2<<<1,blocks,blocks*sizeof(int)>>>(blocks, d_circle_count);

    cudaMemcpy(circle_count, d_circle_count, blocks * sizeof(int), cudaMemcpyDeviceToHost);

    pi = 4.0 * (double)circle_count[0] / (double)samples;

    // Free memory
    free(x);
    free(y);
    free(circle_count);
    cudaFree(d_x);
    cudaFree(d_y);
    cudaFree(d_circle_count);

    return pi;
}

int main(void)
{
    // Initiate variables
    double pi;
    int it;
    long int N = 1048576;
    pi=0.0;

    for(int i=0; i<10024;i++){
        if(i==0)cudaProfilerStart();
        pi += gpu_calc_pi_monte_carlo(N);
        it = i+1;
        if (i%50==1){
            printf("Samples/1 mill : %d, ",it);
            printf("Pi Estimated : %f, ", pi/it);
            printf("Error : %f\n", (PI-(pi/it))/PI);
        } 
        if(i==100)cudaProfilerStop();
    }
    pi /= 10024;
    printf("Pi is estimated to be %f\n", pi);
}