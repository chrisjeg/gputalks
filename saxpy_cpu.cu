#include <stdio.h>

void cpu_saxpy(int n, float a, float*x, float *y)
{
  for(int i = 0; i < n; i++){
    y[i] = a*x[i] + y[i];
  }
}

int main(void)
{
  int N = 1<<20;
  float *x, *y;

  x = (float*)malloc(N*sizeof(float));
  y = (float*)malloc(N*sizeof(float));




  for (int i = 0; i < N; i++) {
    x[i] = 1.0f;
    y[i] = 2.0f;
  }

 


  // Perform SAXPY on 1M elements
  cpu_saxpy(N, 2.0f, x, y);



  free(x);
  free(y);


}