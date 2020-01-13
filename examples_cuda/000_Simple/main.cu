#include <cstdio>
#include <cstdlib>
#include "cuda_runtime.h"

static void vecDouble(int *in, int*out, const int n);

int main()
{
  printf("Hello\n");

  const int n = 10;
  int *in = new int[n];
  int *out = new int[n];
  int *answer = new int[n];

  for (int i = 0; i < n; i++) in[i] = rand() % 100;
  for (int i = 0; i < n; i++) answer[i] = in[i] * 2;

  vecDouble(in, out, n);

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

__global__
void cuda_VecDouble(int *in, int *out, const int n)
{
  const unsigned int i = threadIdx.x;
  if (i >= n) { return; }
  out[i] = in[i] * 2;
}

static void vecDouble(int *hIn, int *hOut, const int n)
{
  int *dIn;
  int *dOut;
  cudaMallocHost((void**)&dIn, n * sizeof(int));
  cudaMallocHost((void**)&dOut, n * sizeof(int));
  cudaMemcpy(dIn, hIn, n * sizeof(int), cudaMemcpyHostToDevice);

  cuda_VecDouble<<<1, n>>>(dIn, dOut, n);
  cudaDeviceSynchronize();

  cudaMemcpy(hOut, dOut, n * sizeof(int), cudaMemcpyDeviceToHost);
  cudaFree(dIn);
  cudaFree(dOut);
}
