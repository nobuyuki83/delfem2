#include <cstdio>
#include <cstdlib>
#include "cuda_runtime.h"

__global__
void kernel_VecScale(float *out, float *in, float scale, const int n)
{
  const unsigned int i = blockIdx.x * blockDim.x + threadIdx.x;
  if (i >= n) { return; }
  out[i] = in[i] * scale;
}

static void cuda_VecScale(float *hOut, float *hIn, float scale, const int n)
{
  float *dOut; cudaMalloc((void**)&dOut, sizeof(float)*n);
  float *dIn;  cudaMalloc((void**)&dIn,  sizeof(float)*n);
  cudaMemcpy(dIn, hIn, sizeof(float)*n, cudaMemcpyHostToDevice);

  unsigned int tpb = 64;
  unsigned int nblk = (unsigned int)((n-1)/tpb+1);
  kernel_VecScale<<<nblk, tpb>>>(dOut, dIn, scale, n);
  cudaDeviceSynchronize();

  cudaMemcpy(hOut, dOut, n * sizeof(float), cudaMemcpyDeviceToHost);
  cudaFree(dOut);
  cudaFree(dIn);
}

int main()
{
  printf("Hello\n");

  const int n = 4000;
  float *in = new float[n];
  float *out = new float[n];
  float *answer = new float[n];

  float scale = 2.0;
  for (int i = 0; i < n; i++) in[i] = rand() % 100;
  for (int i = 0; i < n; i++) answer[i] = in[i] * scale;

  cuda_VecScale(out, in, scale, n);

  for (int i = 0; i < n; i++) {
      if (answer[i] != out[i]) {
          printf("error at index = %d\n", i);
          break;
      }
  }
  printf("OK\n");

  delete[] in;
  delete[] out;
  delete[] answer;

  return 0;
}

