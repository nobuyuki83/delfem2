#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "cuda_runtime.h"

__global__
void cuda_Dot_TPB64(float* d_res, const float* d_A, const float* d_B, int n) {
  const int idx = blockDim.x * blockIdx.x + threadIdx.x;
  if( idx >= n ){ return; }

  const int s_idx = threadIdx.x;

  const unsigned int TPB = 64;
  __shared__ float s_prod[TPB];
  s_prod[s_idx] = d_A[idx]*d_B[idx];
  __syncthreads();

  if( s_idx == 0 ) {
    float blockSum = 0;
    for(int j=0;j<blockDim.x;++j){
      blockSum += s_prod[j];
    }
    atomicAdd(d_res, blockSum);
  }
}


int main()
{
  const unsigned int n = 64*16;
  float h_A[n];
  float h_B[n];
  for(int i=0;i<n;++i){
    h_A[i] = (float)i;
    h_B[i] = (float)i;
  }

  // ------------------------------------------------

  clock_t time0 = clock();

  float h_res = 0.0;
  {
    float *d_A, *d_B, *d_res;
    cudaMalloc((void **) &d_A, sizeof(float) * n);
    cudaMalloc((void **) &d_B, sizeof(float) * n);
    cudaMalloc((void **) &d_res, sizeof(float));

    cudaMemset((void **) &d_res, 0.f, sizeof(float));
    cudaMemcpy(d_A, h_A, sizeof(float) * n, cudaMemcpyHostToDevice);
    cudaMemcpy(d_B, h_B, sizeof(float) * n, cudaMemcpyHostToDevice);

    const unsigned int BLOCK = 64;
    dim3 grid(n / BLOCK);
    dim3 block(BLOCK);

    cuda_Dot_TPB64 << < grid, block >> > (d_res, d_A, d_B, n);

    cudaMemcpy(&h_res, d_res, sizeof(float), cudaMemcpyDeviceToHost);

    cudaFree(d_A);
    cudaFree(d_B);
    cudaFree(d_res);
  }

  printf("gpu sum %.2f \n",h_res);

  clock_t time1 = clock();

  printf("on gpu it takes %.2f sec \n",(double)(time1-time0)/CLOCKS_PER_SEC);

  {
    float sum = 0;
    for(int i=0;i<n;++i){
      sum += h_A[i]*h_B[i];
    }
    printf("cpu sum %.2f \n",h_res);
  }
}
