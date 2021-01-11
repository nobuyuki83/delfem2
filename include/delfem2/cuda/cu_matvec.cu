/*
 * Copyright (c) 2019 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */


#include <cstdio>
#include <cstdlib>
#include <cassert>
#include <cstdint>
#include <thrust/device_vector.h>
#include <thrust/copy.h>
#include <thrust/sort.h>
#include "cuda_runtime.h"

#include "cu_matvec.h"

namespace dfm2 = delfem2;

// -------------------------------------------------------------------

// -------------------------------------------------------------------

__global__
void kernel_VecScale(
    float *out,
    const float *in,
    float scale,
    const int n)
{
  const unsigned int i = blockIdx.x * blockDim.x + threadIdx.x;
  if (i >= n) { return; }
  out[i] = in[i] * scale;
}

void dfm2::cuda::cuda_VecScale(
    float *hOut,
    const float *hIn,
    float scale,
    const unsigned int n)
{
  const thrust::device_vector<float> dIn(hIn,hIn+n);
  thrust::device_vector<float> dOut(n);
  {
    const unsigned int tpb = 64;
    const auto nblk = (unsigned int) ((n - 1) / tpb + 1);
    kernel_VecScale <<< nblk, tpb >>> (
        thrust::raw_pointer_cast(dOut.data()),
        thrust::raw_pointer_cast(dIn.data()),
        scale, n);
  }
  thrust::copy(dOut.begin(),dOut.end(), hOut);
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
    unsigned int n)
{
  const unsigned int idx = blockDim.x * blockIdx.x + threadIdx.x;
  if( idx >= n ){ return; }

  const unsigned int s_idx = threadIdx.x;

  const unsigned int TPB = 64;
  __shared__ float s_prod[TPB];
  s_prod[s_idx] = d_A[idx]*d_B[idx];
  __syncthreads();

  if( s_idx == 0 ) {
    int ns = TPB;
    if( blockDim.x * (blockIdx.x+1) > n ) {
      ns = n - blockDim.x * blockIdx.x;
    }
    float blockSum = 0;
    for(int j=0;j<ns;++j){
      blockSum += s_prod[j];
    }
    atomicAdd(d_res, blockSum);
  }
}


float dfm2::cuda::cuda_Dot(
    const float* h_A,
    const float* h_B,
    unsigned int n)
{
  const thrust::device_vector<float> d_A(h_A, h_A+n);
  const thrust::device_vector<float> d_B(h_B, h_B+n);
  thrust::device_vector<float> d_R(1, 0.f);
  {
    const unsigned int BLOCK = 64;
    dim3 grid((n - 1) / BLOCK + 1);
    dim3 block(BLOCK);
    kernel_Dot_TPB64 <<< grid, block >>> (
        thrust::raw_pointer_cast(d_R.data()),
        thrust::raw_pointer_cast(d_A.data()),
        thrust::raw_pointer_cast(d_B.data()),
        n);
  }
  float h_res;
  thrust::copy(d_R.begin(), d_R.end(), &h_res);
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
  const unsigned int r = blockDim.y * blockIdx.y + threadIdx.y;
  const unsigned int c = blockDim.x * blockIdx.x + threadIdx.x;
  const unsigned int BLOCK = 16;
  assert(blockDim.x == BLOCK && blockDim.y == BLOCK);
  __shared__ float s_A[BLOCK][BLOCK];
  __shared__ float s_B[BLOCK][BLOCK];
  float tmp = 0.0;
  for(int i=0;i<N;i+=BLOCK){
    if( i+threadIdx.x < N && r < N ) {
      s_A[threadIdx.y][threadIdx.x] = A[N * r + i + threadIdx.x];
    }
    if( i+threadIdx.y < N && c < N ) {
      s_B[threadIdx.y][threadIdx.x] = B[N * (i + threadIdx.y) + c];
    }
    __syncthreads(); // wait for copy is finished for all the thread in the block
    int ns = BLOCK;
    if( i+BLOCK >= N ) { ns = N-i; }
    for(int j=0;j<ns;++j) {
      tmp += s_A[threadIdx.y][j] * s_B[j][threadIdx.x];
    }
    __syncthreads();
  }
  if( r >= N || c >= N ){ return; }
  C[N*r+c] = tmp;
}

void dfm2::cuda::cuda_MatMat(
    float *h_C_gpu,
    const float *h_A,
    const float *h_B,
    unsigned int WIDTH)
{
  const thrust::device_vector<float> d_A(h_A,h_A+WIDTH*WIDTH);
  const thrust::device_vector<float> d_B(h_B,h_B+WIDTH*WIDTH);
  thrust::device_vector<float> d_C(WIDTH*WIDTH);
  {
    const unsigned int BLOCK = 16;
    dim3 grid((WIDTH - 1) / BLOCK + 1, (WIDTH - 1) / BLOCK + 1);
    dim3 block(BLOCK, BLOCK);
    kernel_MatMat_TPB16 << < grid, block >> > (
        thrust::raw_pointer_cast(d_C.data()),
//        d_C1,
        thrust::raw_pointer_cast(d_A.data()),
        thrust::raw_pointer_cast(d_B.data()),
        WIDTH);

  }
  thrust::copy(d_C.begin(), d_C.end(), h_C_gpu);
}
