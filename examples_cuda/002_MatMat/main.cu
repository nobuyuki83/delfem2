#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "cuda_runtime.h"

const unsigned int WIDTH = 512;
const unsigned int BLOCK = 16;

__global__ void d_multiply0(float *A, float *B, float *C){
  unsigned int r = blockDim.y * blockIdx.y + threadIdx.y;
  unsigned int c = blockDim.x * blockIdx.x + threadIdx.x;
  float tmp = 0.0f;
  for(int i=0;i<WIDTH;++i){
    tmp += A[WIDTH*r+i] * B[WIDTH*i+c];
  }
  C[WIDTH*r+c] = tmp;
}

__global__ void d_multiply1(float *A, float *B, float *C){
  unsigned int r = blockDim.y * blockIdx.y + threadIdx.y;
  unsigned int c = blockDim.x * blockIdx.x + threadIdx.x;
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


void h_multiply(float *A, float* B, float *C){
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
  printf("hoge B\n");

  float h_A[WIDTH*WIDTH];
  float h_B[WIDTH*WIDTH];
  float h_C[WIDTH*WIDTH];
  for(int i=0;i<WIDTH*WIDTH;++i){
    h_A[i] = (float)i;
    h_B[i] = (float)i;
  }

  // ------------------------------------------------

  clock_t time0 = clock();

  float *d_A, *d_B, *d_C;
  cudaMalloc((void**)&d_A, sizeof(float)*WIDTH*WIDTH);
  cudaMalloc((void**)&d_B, sizeof(float)*WIDTH*WIDTH);
  cudaMalloc((void**)&d_C, sizeof(float)*WIDTH*WIDTH);

  cudaMemcpy(d_A, h_A, sizeof(float)*WIDTH*WIDTH, cudaMemcpyHostToDevice);
  cudaMemcpy(d_B, h_B, sizeof(float)*WIDTH*WIDTH, cudaMemcpyHostToDevice);

  dim3 grid(WIDTH/BLOCK,WIDTH/BLOCK);
  dim3 block(BLOCK,BLOCK);

  d_multiply0 <<< grid, block >>> (d_A,d_B,d_C);

  d_multiply1 <<< grid, block >>> (d_A,d_B,d_C);

  cudaMemcpy(h_C, d_C, sizeof(float)*WIDTH*WIDTH, cudaMemcpyDeviceToHost);

  cudaFree(d_A);
  cudaFree(d_B);
  cudaFree(d_C);

  clock_t time1 = clock();

  printf("on gpu it takes %.2f sec \n",(double)(time1-time0)/CLOCKS_PER_SEC);

  h_multiply(h_A, h_B, h_C);

  clock_t time2 = clock();

  printf("on cpu it takes %.2f sec \n",(double)(time2-time1)/CLOCKS_PER_SEC);
}
