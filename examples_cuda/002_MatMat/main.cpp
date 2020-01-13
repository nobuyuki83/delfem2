#include <cstdio>
#include <cstdlib>
#include <cassert>
#include <ctime>
#include "cuda_runtime.h"

#include "delfem2/cuda/cu_matvec.h"

namespace dfm2 = delfem2;

void h_multiply(
    float *C,
    const float* A,
    const float *B,
    unsigned int WIDTH)
{
  for(int r=0;r<WIDTH;++r) {
    for(int c=0;c<WIDTH;++c){
      float tmp = 0.0;
      for(int i=0;i<WIDTH;++i){
        tmp += A[WIDTH*r+i] * B[WIDTH*i+c];
      }
      C[WIDTH*r+c] = tmp;
    }
  }
}

int main()
{
  const unsigned int WIDTH = 512;

  float h_A[WIDTH*WIDTH];
  float h_B[WIDTH*WIDTH];
  for(int i=0;i<WIDTH*WIDTH;++i){
    h_A[i] = (float)(i%WIDTH);
    h_B[i] = (float)(i%WIDTH);
  }

  clock_t time0 = clock();

  // ------------------------------------------------

  float h_C_cpu[WIDTH*WIDTH];
  h_multiply(h_C_cpu,
             h_A, h_B,WIDTH);

  clock_t time1 = clock();

  printf("on cpu it takes %.2f sec \n",(double)(time1-time0)/CLOCKS_PER_SEC);

  // ------------------------------------------------

  float h_C_gpu[WIDTH*WIDTH];
  dfm2::cuda::cuda_MatMat(h_C_gpu,
      h_A, h_B, WIDTH);

  clock_t time2 = clock();

  printf("on gpu it takes %.2f sec \n",(double)(time2-time1)/CLOCKS_PER_SEC);

  {
    float sumDiff = 0.0;
    for (int i = 0; i < WIDTH; ++i) {
      for (int j = 0; j < WIDTH; ++j) {
        float d0 = h_C_cpu[i * WIDTH + j] - h_C_gpu[i * WIDTH + j];
        sumDiff += d0*d0;
      }
    }
    printf("sumDiff: %f \n", sumDiff);
  }

}
