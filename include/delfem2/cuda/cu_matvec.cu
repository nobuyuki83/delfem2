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

__device__
void device_AtomicMaxFloat(float * const address, const float value)
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
void device_AtomicMinFloat(float * const address, const float value)
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

__device__
float kernel_dist3(
    const float p0[3],
    const float p1[3])
{
  float v = (p1[0]-p0[0])*(p1[0]-p0[0]) + (p1[1]-p0[1])*(p1[1]-p0[1]) + (p1[2]-p0[2])*(p1[2]-p0[2]);
  return sqrtf(v);
}


__device__
unsigned int device_ExpandBits(unsigned int v)
{
  v = (v * 0x00010001u) & 0xFF0000FFu;
  v = (v * 0x00000101u) & 0x0F00F00Fu;
  v = (v * 0x00000011u) & 0xC30C30C3u;
  v = (v * 0x00000005u) & 0x49249249u;
  return v;
}

__device__
unsigned int device_MortonCode(float x, float y, float z)
{
  auto ix = (unsigned int)fmin(fmax(x * 1024.0f, 0.0f), 1023.0f);
  auto iy = (unsigned int)fmin(fmax(y * 1024.0f, 0.0f), 1023.0f);
  auto iz = (unsigned int)fmin(fmax(z * 1024.0f, 0.0f), 1023.0f);
  //  std::cout << std::bitset<10>(ix) << " " << std::bitset<10>(iy) << " " << std::bitset<10>(iz) << std::endl;
  ix = device_ExpandBits(ix);
  iy = device_ExpandBits(iy);
  iz = device_ExpandBits(iz);
  //  std::cout << std::bitset<30>(ix) << " " << std::bitset<30>(iy) << " " << std::bitset<30>(iz) << std::endl;
  return ix * 4 + iy * 2 + iz;
}

__device__
int delta(unsigned int* sortedMortonCodes, int x, int y, int numObjects)
{
  if (x >= 0 && x <= numObjects - 1 && y >= 0 && y <= numObjects - 1)
  {
//         return __clz(sortedMortonCodes[x] ^ sortedMortonCodes[y]);
#if HASH_64
    return __clzll(sortedMortonCodes[x] ^ sortedMortonCodes[y]);
#else
    return __clz(sortedMortonCodes[x] ^ sortedMortonCodes[y]);
#endif
  }
  return -1;
}

__device__
int2 determineRange(unsigned int* sortedMortonCodes, int numObjects, int idx)
{
  int d = delta(sortedMortonCodes, idx, idx + 1, numObjects) - delta(sortedMortonCodes, idx, idx - 1, numObjects);
  d = d > 0 ? 1 : -1;

  //compute the upper bound for the length of the range
  int dmin = delta(sortedMortonCodes, idx, idx - d, numObjects);
  int lmax = 2;
  while (delta(sortedMortonCodes, idx, idx + lmax * d, numObjects) > dmin){
    lmax = lmax * 2;
  }

  //find the other end using binary search
  int l = 0;
  for (int t = lmax / 2; t >= 1; t /= 2)
  {
    if (delta(sortedMortonCodes, idx, idx + (l + t)*d, numObjects) > dmin)
    l += t;
  }
  int j = idx + l*d;

  int2 range;
  range.x = min(idx, j);
  range.y = max(idx, j);
  if (idx == 38043 || idx == 38044 || idx == 38045 || idx == 38046 || idx == 38047 || idx == 38048)
  printf("idx %d range :%d - %d j: %d morton: %d\n", idx, range.x, range.y, j, sortedMortonCodes[idx]);
  return range;
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
    const unsigned int n)
{
  const thrust::device_vector<float> dIn(hIn,hIn+n);
  thrust::device_vector<float> dOut(n);
  {
    const unsigned int tpb = 64;
    const unsigned int nblk = (unsigned int) ((n - 1) / tpb + 1);
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
    device_AtomicMinFloat(d_minmax+idy+0,vmin);
    device_AtomicMaxFloat(d_minmax+idy+3,vmax);
  }
}

void dfm2::cuda::cuda_Min3Max3_Points3F(
    float *h_min3,
    float *h_max3,
    const float *h_XYZ,
    unsigned int np)
{
  h_min3[0] = h_max3[0] = h_XYZ[0];
  h_min3[1] = h_max3[1] = h_XYZ[1];
  h_min3[2] = h_max3[2] = h_XYZ[2];
  // --------------------------------------
  const thrust::device_vector<float> d_XYZ(h_XYZ,h_XYZ+np*3);
  thrust::device_vector<float> d_minmax(6);
  {
    const unsigned int BLOCK = 256;
    dim3 grid((np - 1) / BLOCK + 1);
    dim3 block(BLOCK, 3);
    kernel_MinMax_TPB256 <<< grid, block >>> (
        thrust::raw_pointer_cast(d_minmax.data()),
        thrust::raw_pointer_cast(d_XYZ.data()),
        np);
  }
  thrust::copy(d_minmax.begin()+0,d_minmax.begin()+3, h_min3);
  thrust::copy(d_minmax.begin()+3,d_minmax.begin()+6, h_max3);
}

// ---------------------------------------------------------------------------------------------------------------------

__global__
void kernel_CentRad_MeshTri3D_TPB64(
    float *dXYZ_c,
    float *dMaxRad,
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
  assert( i0 < nXYZ && i1 < nXYZ && i2 < nXYZ );
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
  const float l0 = kernel_dist3(pc, p0);
  const float l1 = kernel_dist3(pc, p1);
  const float l2 = kernel_dist3(pc, p2);
  float lm = l0;
  if( l1 > lm ){ lm = l1; }
  if( l2 > lm ){ lm = l2; }
  const unsigned int TPB = 64;
  assert( blockDim.x == TPB );
  __shared__ float sRad[TPB];
  const unsigned int s_idx = threadIdx.x;
  sRad[s_idx] = lm;
  __syncthreads();
  if( s_idx == 0 ) {
    int ns = TPB;
    if( blockDim.x * (blockIdx.x+1) > nTri ) {
      ns = nTri - blockDim.x * blockIdx.x;
    }
    float blockRadMax = sRad[0];
    for(int ins=1;ins<ns;++ins){
      if( sRad[ins] > blockRadMax ){ blockRadMax = sRad[ins]; }
    }
    device_AtomicMaxFloat(dMaxRad, blockRadMax);
  }
}

void dfm2::cuda::cuda_CentsMaxRad_MeshTri3F(
    float* hXYZ_c,
    float* hMaxRad,
    const float *hXYZ,
    const unsigned int nXYZ,
    const unsigned int *hTri,
    const unsigned int nTri)
{
  const thrust::device_vector<float> dXYZ(hXYZ, hXYZ+nXYZ*3);
  const thrust::device_vector<unsigned int> dTri(hTri, hTri+nTri*3);
  thrust::device_vector<float> dMaxRad(1);
  thrust::device_vector<float> dXYZ_c(nTri*3);
  {
    const unsigned int BLOCK = 64;
    dim3 grid( (nTri-1)/BLOCK + 1 );
    dim3 block( BLOCK );
    kernel_CentRad_MeshTri3D_TPB64 <<< grid, block >>> (
        thrust::raw_pointer_cast(dXYZ_c.data()),
        thrust::raw_pointer_cast(dMaxRad.data()),
        thrust::raw_pointer_cast(dXYZ.data()),
        nXYZ,
        thrust::raw_pointer_cast(dTri.data()),
        nTri);
  }
  thrust::copy(dXYZ_c.begin(), dXYZ_c.end(), hXYZ_c);
  thrust::copy(dMaxRad.begin(), dMaxRad.end(), hMaxRad);
}

// --------------------------------------------------------------------------------

__global__
void kernel_MortonCodeId_Points3F_TPB64(
    unsigned int *dMC,
    unsigned int *dId,
    const float *dXYZ,
    const unsigned int nXYZ,
    const float* min_xyz,
    const float* max_xyz)
{
  const unsigned int idx = blockDim.x * blockIdx.x + threadIdx.x;
  if( idx >= nXYZ ) return;
  dId[idx] = idx;
  // ----------------------------
  const float x0 = (dXYZ[idx*3+0]-min_xyz[0])/(max_xyz[0]-min_xyz[0]);
  const float y0 = (dXYZ[idx*3+1]-min_xyz[1])/(max_xyz[1]-min_xyz[1]);
  const float z0 = (dXYZ[idx*3+2]-min_xyz[2])/(max_xyz[2]-min_xyz[2]);
  unsigned int mc = device_MortonCode(x0,y0,z0);
  dMC[idx] = mc;
}

void dfm2::cuda::cuda_MortonCode_Points3FSorted(
    unsigned int *hSortedId,
    std::uint32_t *hSortedMc,
    const float *hXYZ,
    const unsigned int nXYZ,
    const float* hMinXYZ,
    const float* hMaxXYZ)
{
  const thrust::device_vector<float> dXYZ(hXYZ, hXYZ+nXYZ*3);
  thrust::device_vector<unsigned int> dMC(nXYZ);
  thrust::device_vector<unsigned int> dId(nXYZ);
  thrust::device_vector<float> dMinXYZ(hMinXYZ,hMinXYZ+6);
  thrust::device_vector<float> dMaxXYZ(hMaxXYZ,hMaxXYZ+6);
  {
    const unsigned int BLOCK = 64;
    dim3 grid( (nXYZ-1)/BLOCK+1 );
    dim3 block( BLOCK );
    kernel_MortonCodeId_Points3F_TPB64 <<< grid, block >>> (
        thrust::raw_pointer_cast(dMC.data()),
        thrust::raw_pointer_cast(dId.data()),
        thrust::raw_pointer_cast(dXYZ.data()),
        nXYZ,
        thrust::raw_pointer_cast(dMinXYZ.data()),
        thrust::raw_pointer_cast(dMaxXYZ.data()));
  }
  thrust::sort_by_key(dMC.begin(),dMC.end(),dId.begin());
  thrust::copy(dMC.begin(), dMC.end(), hSortedMc);
  thrust::copy(dId.begin(), dId.end(), hSortedId);
}

// ------------------------------------------------

__device__
int device_Delta(int i, int j, const unsigned int* sortedMC, int nMC)
{
  if ( j<0 || j >= nMC ){ return -1; }
  return __clz(sortedMC[i] ^ sortedMC[j]);
}

__device__
int2 device_MortonCode_DeterminRange
(const unsigned int* sortedMC,
 int nMC,
 int imc)
{
  if( imc == 0 ){ return make_int2(0,nMC-1); }
  // ----------------------
  const std::uint32_t mc0 = sortedMC[imc-1];
  const std::uint32_t mc1 = sortedMC[imc+0];
  const std::uint32_t mc2 = sortedMC[imc+1];
  if( mc0 == mc1 && mc1 == mc2 ){ // for hash value collision
    int jmc=imc+1;
    for(;jmc<nMC;++jmc){
      if( sortedMC[jmc] != mc1 ) break;
    }
    return make_int2(imc,jmc-1);
  }
  int d = device_Delta(imc, imc + 1, sortedMC, nMC) - device_Delta(imc, imc - 1, sortedMC, nMC);
  d = d > 0 ? 1 : -1;

  //compute the upper bound for the length of the range
  const int delta_min = device_Delta(imc, imc - d, sortedMC, nMC);
  int lmax = 2;
  while (device_Delta(imc, imc + lmax*d, sortedMC, nMC)>delta_min)
  {
    lmax = lmax * 2;
  }

  //find the other end using binary search
  int l = 0;
  for (int t = lmax / 2; t >= 1; t /= 2)
  {
    if (device_Delta(imc, imc + (l + t)*d, sortedMC, nMC)>delta_min)
    {
      l = l + t;
    }
  }
  int j = imc + l*d;

  int2 range = make_int2(-1,-1);
  if (imc <= j) { range.x = imc; range.y = j; }
  else { range.x = j; range.y = imc; }
  return range;
}

__device__
int device_MortonCode_FindSplit
(const unsigned int* sortedMC,
 unsigned int iMC_start,
 unsigned int iMC_last)
{
  //return -1 if there is only
  //one primitive under this node.
  if (iMC_start == iMC_last) { return -1; }

  // ------------------------------
  const int common_prefix = __clz(sortedMC[iMC_start] ^ sortedMC[iMC_last]);

  //handle duplicated morton code
  if (common_prefix == 32 ){ return iMC_start; } // sizeof(std::uint32_t)*8

  // Use binary search to find where the next bit differs.
  // Specifically, we are looking for the highest object that
  // shares more than commonPrefix bits with the first one.
  const std::uint32_t mcStart = sortedMC[iMC_start];
  int iMC_split = iMC_start; // initial guess
  int step = iMC_last - iMC_start;
  do
  {
    step = (step + 1) >> 1; // exponential decrease
    const int newSplit = iMC_split + step; // proposed new position
    if (newSplit < iMC_last){
      const unsigned int splitCode = sortedMC[newSplit];
      int splitPrefix = __clz(mcStart ^ splitCode);
      if (splitPrefix > common_prefix){
        iMC_split = newSplit; // accept proposal
      }
    }
  }
  while (step > 1);
  return iMC_split;
}

__global__
void kernel_MortonCode_BVHTopology_TPB64(
dfm2::CNodeBVH2* dNodeBVH,
const unsigned int *dSortedMC,
const unsigned int *dSortedId,
const unsigned int nMC)
{
  const unsigned int idx = blockDim.x * blockIdx.x + threadIdx.x;
  if (idx >= nMC-1) return;
  const unsigned int ini = idx;
  const unsigned int nni = nMC-1;
  // -------------------------------
  const int2 range = device_MortonCode_DeterminRange(dSortedMC,nMC,ini);
  const int isplit = device_MortonCode_FindSplit(dSortedMC,range.x,range.y);
//  printf("%d --> %d %d  %d\n",ini, range.x, range.y, isplit);
  // -------------------------------
  if( range.x == isplit ){
    const int inlA = nni+isplit;
    dNodeBVH[ini].ichild[0] = inlA;
    dNodeBVH[inlA].iroot = ini;
    dNodeBVH[inlA].ichild[0] = dSortedId[isplit];
    dNodeBVH[inlA].ichild[1] = -1;
  }
  else{
    const int iniA = isplit;
    dNodeBVH[ini].ichild[0] = iniA;
    dNodeBVH[iniA].iroot = ini;
  }
  // ----
  if( range.y == isplit+1 ){
    const int inlB = nni+isplit+1;
    dNodeBVH[ini].ichild[1] = inlB;
    dNodeBVH[inlB].iroot = ini;
    dNodeBVH[inlB].ichild[0] = dSortedId[isplit+1];
    dNodeBVH[inlB].ichild[1] = -1;
  }
  else{
    const int iniB = isplit+1;
    dNodeBVH[ini].ichild[1] = iniB;
    dNodeBVH[iniB].iroot = ini;
  }
}




void dfm2::cuda::cuda_MortonCode_BVHTopology(
  CNodeBVH2* hNodeBVH,
  const unsigned int* aSortedId,
  const std::uint32_t* aSortedMc,
  unsigned int N)
{
  const thrust::device_vector<std::uint32_t> dMC(aSortedMc,aSortedMc+N);
  const thrust::device_vector<unsigned int> dId(aSortedId,aSortedId+N);
  thrust::device_vector<dfm2::CNodeBVH2> dNodeBVH(N*2-1);
  // ----------------------------------
  {
    const unsigned int BLOCK = 64;
    dim3 grid( (N-1)/BLOCK+1 );
    dim3 block( BLOCK );
    kernel_MortonCode_BVHTopology_TPB64 <<< grid, block >>> (
      thrust::raw_pointer_cast(dNodeBVH.data()),
      thrust::raw_pointer_cast(dMC.data()),
      thrust::raw_pointer_cast(dId.data()),
      N);
  }
  thrust::copy(dNodeBVH.begin(), dNodeBVH.end(), hNodeBVH);
}
