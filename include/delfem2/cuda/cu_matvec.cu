#include <cstdio>
#include <cstdlib>
#include <cassert>
#include "cuda_runtime.h"

#include "cu_matvec.h"

__global__
void kernel_VecScale(float *out, float *in, float scale, const int n)
{
  const unsigned int i = blockIdx.x * blockDim.x + threadIdx.x;
  if (i >= n) { return; }
  out[i] = in[i] * scale;
}

void cuda_VecScale(float *hOut, float *hIn, float scale, const int n)
{
  float *dOut; cudaMalloc((void**)&dOut, sizeof(float)*n);
  float *dIn;  cudaMalloc((void**)&dIn,  sizeof(float)*n);
  cudaMemcpy(dIn, hIn, sizeof(float)*n, cudaMemcpyHostToDevice);

  const unsigned int tpb = 64;
  const unsigned int nblk = (unsigned int)((n-1)/tpb+1);
  kernel_VecScale<<<nblk, tpb>>>(dOut, dIn, scale, n);
  cudaDeviceSynchronize();

  cudaMemcpy(hOut, dOut, n * sizeof(float), cudaMemcpyDeviceToHost);
  cudaFree(dOut);
  cudaFree(dIn);
}


// -------------------------------------------------------------------

/**
 * @brief dot product of two vectors
 */
__global__
void kernel_Dot_TPB64(
    float* d_res,
    const float* d_A,
    const float* d_B,
    int n)
{
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


float cuda_Dot(
    const float* h_A,
    const float* h_B,
    unsigned int n)
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

  kernel_Dot_TPB64 << < grid, block >> > (d_res, d_A, d_B, n);

  float h_res;
  cudaMemcpy(&h_res, d_res, sizeof(float), cudaMemcpyDeviceToHost);

  cudaFree(d_A);
  cudaFree(d_B);
  cudaFree(d_res);
  return h_res;
}

// ------------------------------------------------------------------

__global__
void kernel_MatMat_TPB16(
    float *C,
    const float *A,
    const float *B,
    unsigned int N)
{
  unsigned int r = blockDim.y * blockIdx.y + threadIdx.y;
  unsigned int c = blockDim.x * blockIdx.x + threadIdx.x;
  const unsigned int BLOCK = 16;
  assert(blockDim.x == BLOCK);
  assert(blockDim.y == BLOCK);
  __shared__ float s_A[BLOCK][BLOCK];
  __shared__ float s_B[BLOCK][BLOCK];
  float tmp = 0.0;
  for(int i=0;i<N;i+=BLOCK){
    s_A[threadIdx.y][threadIdx.x] = A[N*r+i+threadIdx.x];
    s_B[threadIdx.y][threadIdx.x] = B[N*(i+threadIdx.y) + c];
    __syncthreads(); // wait for copy is finished for all the thread in the block
    for(int j=0;j<BLOCK;++j) {
      tmp += s_A[threadIdx.y][j] * s_B[j][threadIdx.x];
    }
    __syncthreads();
  }
  C[N*r+c] = tmp;
}

void cuda_MatMat(
    float* h_C_gpu,
    const float* h_A,
    const float* h_B,
    unsigned int WIDTH)
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