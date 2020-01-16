#include <cstdio>
#include <cstdlib>
#include <cassert>
#include "cuda_runtime.h"

#include "cu_matvec.h"

namespace dfm2 = delfem2;

// -------------------------------------------------------------------

__device__
void atomicMaxFloat(float * const address, const float value)
{
  if ( *address >= value ) { return; }
  int * const address_as_i = (int *)address;
  int old = * address_as_i, assumed;
  do {
    assumed = old;
    if (__int_as_float(assumed) >= value) { break; }
    old = atomicCAS(address_as_i, assumed, __float_as_int(value));
  } while (assumed != old);
}

__device__
void atomicMinFloat(float * const address, const float value)
{
  if ( *address <= value ) { return; }
  int * const address_as_i = (int *)address;
  int old = * address_as_i, assumed;
  do {
    assumed = old;
    if (__int_as_float(assumed) <= value) { break; }
    old = atomicCAS(address_as_i, assumed, __float_as_int(value));
  } while (assumed != old);
}


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
    const int n)
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
  float *d_A, *d_B, *d_res;
  cudaMalloc((void **) &d_A, sizeof(float) * n);
  cudaMalloc((void **) &d_B, sizeof(float) * n);
  cudaMalloc((void **) &d_res, sizeof(float));

  cudaMemset((void **) &d_res, 0.f, sizeof(float));
  cudaMemcpy(d_A, h_A, sizeof(float) * n, cudaMemcpyHostToDevice);
  cudaMemcpy(d_B, h_B, sizeof(float) * n, cudaMemcpyHostToDevice);

  const unsigned int BLOCK = 64;
  dim3 grid( (n-1)/BLOCK + 1);
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
  float *d_A, *d_B, *d_C;
  cudaMalloc((void **) &d_A, sizeof(float) * WIDTH * WIDTH);
  cudaMalloc((void **) &d_B, sizeof(float) * WIDTH * WIDTH);
  cudaMalloc((void **) &d_C, sizeof(float) * WIDTH * WIDTH);

  cudaMemcpy(d_A, h_A, sizeof(float) * WIDTH * WIDTH, cudaMemcpyHostToDevice);
  cudaMemcpy(d_B, h_B, sizeof(float) * WIDTH * WIDTH, cudaMemcpyHostToDevice);

  const unsigned int BLOCK = 16;
  dim3 grid( (WIDTH-1)/BLOCK+1, (WIDTH-1)/BLOCK+1);
  dim3 block(BLOCK, BLOCK);

  kernel_MatMat_TPB16 << < grid, block >> > (d_C, d_A, d_B, WIDTH);

  cudaMemcpy(h_C_gpu,
             d_C, sizeof(float) * WIDTH * WIDTH, cudaMemcpyDeviceToHost);

  cudaFree(d_A);
  cudaFree(d_B);
  cudaFree(d_C);
}

// ------------------------------------------------------------------------

__global__
void kernel_MinMax_TPB256(
    float *d_minmax,
    const float *d_XYZ,
    unsigned int np)
{
  const unsigned int idx = blockDim.x * blockIdx.x + threadIdx.x;
  const unsigned int s_idx = threadIdx.x;
  assert( blockDim.y == 3 && blockIdx.y == 0 );
  unsigned int idy = threadIdx.y;
  if( idx >= np ){ return; }
  // ---------------
  const unsigned int BLOCK = 256;
  assert(blockDim.x == BLOCK);
  __shared__ float s_XYZ[BLOCK][3];
  s_XYZ[s_idx][idy] = d_XYZ[idx*3+idy];
  __syncthreads();
  if( s_idx == 0 ) {
    float vmin = s_XYZ[0][idy];
    float vmax = s_XYZ[0][idy];
    int ns = BLOCK;
    if( blockDim.x * (blockIdx.x+1) > np ) {
      ns = np - blockDim.x * blockIdx.x;
    }
    for(int is=0;is<ns;++is){
      if( s_XYZ[is][idy] < vmin ){ vmin = s_XYZ[is][idy]; }
      if( s_XYZ[is][idy] > vmax ){ vmax = s_XYZ[is][idy]; }
    }
    atomicMinFloat(d_minmax+idy+0,vmin);
    atomicMaxFloat(d_minmax+idy+3,vmax);
  }
}

void dfm2::cuda::cuda_MinMax_Point3D(
    float *h_minmax,
    const float *h_XYZ,
    unsigned int np)
{
  h_minmax[0] = h_minmax[3] = h_XYZ[0];
  h_minmax[1] = h_minmax[4] = h_XYZ[1];
  h_minmax[2] = h_minmax[5] = h_XYZ[2];
  // --------------------------------------
  float *d_minmax, *d_XYZ;
  cudaMalloc((void **) &d_minmax, sizeof(float) * 6);
  cudaMalloc((void **) &d_XYZ, sizeof(float) * np * 3);
  cudaMemcpy(d_XYZ,
             h_XYZ, sizeof(float) * np * 3, cudaMemcpyHostToDevice);
  cudaMemcpy(d_minmax,
             h_minmax, sizeof(float) * 6, cudaMemcpyHostToDevice);

  {
    const unsigned int BLOCK = 256;
    dim3 grid((np - 1) / BLOCK + 1);
    dim3 block(BLOCK, 3);
    kernel_MinMax_TPB256 <<< grid, block >>> (d_minmax, d_XYZ, np);
  }

  cudaMemcpy(h_minmax,
             d_minmax, sizeof(float) * 6, cudaMemcpyDeviceToHost);

  cudaFree(d_minmax);
  cudaFree(d_XYZ);
}

// ---------------------------------------------------------------------------------------------------------------------

__device__
void kernel_dist3(
    float *d,
    const float p0[3],
    const float p1[3])
{
  float v = (p1[0]-p0[0])*(p1[0]-p0[0]) + (p1[1]-p0[1])*(p1[1]-p0[1]) + (p1[2]-p0[2])*(p1[2]-p0[2]);
  *d = sqrtf(v);
}

__global__
void kernel_CentRad_MeshTri3D_TPB256(
    float *dXYZ_c,
    float *dRad,
    const float *dXYZ,
    const unsigned int nXYZ,
    const unsigned int *dTri,
    const unsigned int nTri)
{
  const unsigned int idx = blockDim.x * blockIdx.x + threadIdx.x;
  if( idx >= nTri ) return;
  // ----------------------------
  const unsigned int itri = idx;
  const unsigned int i0 = dTri[itri*3+0];
  const unsigned int i1 = dTri[itri*3+1];
  const unsigned int i2 = dTri[itri*3+2];
  const float p0[3] = {dXYZ[i0*3+0],dXYZ[i0*3+1],dXYZ[i0*3+2]};
  const float p1[3] = {dXYZ[i1*3+0],dXYZ[i1*3+1],dXYZ[i1*3+2]};
  const float p2[3] = {dXYZ[i2*3+0],dXYZ[i2*3+1],dXYZ[i2*3+2]};
  const float pc[3] = {
      (p0[0]+p1[0]+p2[0])/3.f,
      (p0[1]+p1[1]+p2[1])/3.f,
      (p0[2]+p1[2]+p2[2])/3.f };
  dXYZ_c[itri*3+0] = pc[0];
  dXYZ_c[itri*3+1] = pc[1];
  dXYZ_c[itri*3+2] = pc[2];
  // ---------------------
  float l0,l1,l2;
  kernel_dist3(&l0, pc, p0);
  kernel_dist3(&l1, pc, p1);
  kernel_dist3(&l2, pc, p2);
  if( l0 > l1 && l0 > l2 ){ dRad[itri] = l0; return; }
  if( l1 > l0 && l1 > l2 ){ dRad[itri] = l1; return; }
  dRad[itri] = l2;
}

void dfm2::cuda::cuda_CentRad_MeshTri3D(
    float* hXYZ_c,
    float* hRad,
    const float *hXYZ,
    const unsigned int nXYZ,
    const unsigned int *hTri,
    const unsigned int nTri)
{
  float *dXYZ, *dXYZ_c, *dRad;
  unsigned int *dTri;
  cudaMalloc((void **) &dXYZ, sizeof(float) * nXYZ * 3);
  cudaMalloc((void **) &dTri, sizeof(unsigned int) * nTri * 3);
  cudaMalloc((void **) &dXYZ_c, sizeof(float) * nTri * 3);
  cudaMalloc((void **) &dRad, sizeof(float) * nTri);
  cudaMemcpy(dXYZ,
             hXYZ, sizeof(float) * nXYZ * 3, cudaMemcpyHostToDevice);
  cudaMemcpy(dTri,
             hTri, sizeof(unsigned int) * nTri * 3, cudaMemcpyHostToDevice);

  {
    const unsigned int BLOCK = 64;
    dim3 grid( (nTri-1)/BLOCK + 1 );
    dim3 block( BLOCK );
    kernel_CentRad_MeshTri3D_TPB256 <<< grid, block >>> (dXYZ_c, dRad,
        dXYZ, nXYZ,
        dTri, nTri);
  }

  cudaMemcpy(hXYZ_c,
             dXYZ_c, sizeof(float) * nTri * 3, cudaMemcpyDeviceToHost);
  cudaMemcpy(hRad,
             dRad, sizeof(float) * nTri, cudaMemcpyDeviceToHost);

  cudaFree(dTri);
  cudaFree(dXYZ);
  cudaFree(dXYZ_c);
  cudaFree(dRad);
}
