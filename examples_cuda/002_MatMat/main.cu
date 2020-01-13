#include <stdio.h>
#include <stdlib.h>
#include <cassert>
#include <time.h>
#include "cuda_runtime.h"

__global__
void d_multiply0(
    float *C,
    const float *A,
    const float *B,
    unsigned int WIDTH)
{
  const unsigned int r = blockIdx.y * blockDim.y + threadIdx.y; // row index
  const unsigned int c = blockIdx.x * blockDim.x + threadIdx.x; // column index
  float tmp = 0.0f;
  for(int k=0;k<WIDTH;++k){
    tmp += A[WIDTH*r+k] * B[WIDTH*k+c];
  }
  C[WIDTH*r+c] = tmp;
}


__global__
void kernel_MatMat_TPB16(
    float *C,
    const float *A,
    const float *B,
    unsigned int WIDTH)
{
  unsigned int r = blockDim.y * blockIdx.y + threadIdx.y;
  unsigned int c = blockDim.x * blockIdx.x + threadIdx.x;
  const unsigned int BLOCK = 16;
  assert(blockDim.x == BLOCK);
  assert(blockDim.y == BLOCK);
  __shared__ float s_A[BLOCK][BLOCK];
  __shared__ float s_B[BLOCK][BLOCK];
  float tmp = 0.0;
  for(int i=0;i<WIDTH;i+=BLOCK){
    s_A[threadIdx.y][threadIdx.x] = A[WIDTH*r+i+threadIdx.x];
    s_B[threadIdx.y][threadIdx.x] = B[WIDTH*(i+threadIdx.y) + c];
    __syncthreads(); // wait for copy is finished for all the thread in the block
    for(int j=0;j<BLOCK;++j) {
      tmp += s_A[threadIdx.y][j] * s_B[j][threadIdx.x];
    }
    __syncthreads();
  }
  C[WIDTH*r+c] = tmp;
}


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

  {
    const unsigned int BLOCK = 16;
    float *d_A, *d_B, *d_C;
    cudaMalloc((void **) &d_A, sizeof(float) * WIDTH * WIDTH);
    cudaMalloc((void **) &d_B, sizeof(float) * WIDTH * WIDTH);
    cudaMalloc((void **) &d_C, sizeof(float) * WIDTH * WIDTH);

    cudaMemcpy(d_A, h_A, sizeof(float) * WIDTH * WIDTH, cudaMemcpyHostToDevice);
    cudaMemcpy(d_B, h_B, sizeof(float) * WIDTH * WIDTH, cudaMemcpyHostToDevice);

    dim3 grid(WIDTH / BLOCK, WIDTH / BLOCK);
    dim3 block(BLOCK, BLOCK);

//    d_multiply0 << < grid, block >> > (d_C, d_A, d_B, WIDTH);
    kernel_MatMat_TPB16 <<< grid, block >>> (d_C, d_A,d_B, WIDTH);

    cudaMemcpy(h_C_gpu,
               d_C, sizeof(float) * WIDTH * WIDTH, cudaMemcpyDeviceToHost);

    cudaFree(d_A);
    cudaFree(d_B);
    cudaFree(d_C);
  }

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
