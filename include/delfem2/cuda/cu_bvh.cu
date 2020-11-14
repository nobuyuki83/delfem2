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

#include "cu_bvh.h"

namespace dfm2 = delfem2;

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


template <typename REAL>
__device__
float device_Distance3(
    const REAL p0[3],
    const REAL p1[3])
{
  const REAL v = (p1[0]-p0[0])*(p1[0]-p0[0]) + (p1[1]-p0[1])*(p1[1]-p0[1]) + (p1[2]-p0[2])*(p1[2]-p0[2]);
  return sqrt(v);
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
int device_Delta(int i, int j, const unsigned int* sortedMC, int nMC)
{
  if ( j<0 || j >= nMC ){ return -1; }
  return __clz(sortedMC[i] ^ sortedMC[j]);
}

__device__
int2 device_MortonCode_DeterminRange(
    const unsigned int* sortedMC,
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
int device_MortonCode_FindSplit(
    const unsigned int* sortedMC,
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

template <typename REAL>
class CudaBV_Sphere
{
public:
  __device__
  void Set_Inactive() { r = -1; }
  __device__
  void AddPoint(const REAL p[3], REAL R){
    if( R < 0 ){ return; }
    if( r < 0 ){ c[0]=p[0]; c[1]=p[1]; c[2]=p[2]; r=R; return; }
    const REAL L = device_Distance3(p,c);
    if( r>L+R ){ return; } // including
    if( R>L+r){ // included
      c[0]=p[0]; c[1]=p[1]; c[2]=p[2]; r=R;
      return;
    }
    if( fabs(L) <= 1.0e-5*fabs(r+R) ){ // almost co-centric
      r = L+R;
      return;
    }
    const REAL r0 = 0.5*(L+r-R)/L;
    const REAL r1 = 0.5*(L+R-r)/L;
    assert( r0 >= 0 && r1 >= 0 );
    c[0] = r0*c[0] + r1*p[0];
    c[1] = r0*c[1] + r1*p[1];
    c[2] = r0*c[2] + r1*p[2];
    r = 0.5*(L+r+R);
    return;
  }
  __device__
  void Add(const CudaBV_Sphere<REAL>& bb) {
    this->AddPoint(bb.c,bb.r);
  }
  __device__
  void Range_DistToPoint(REAL& min0, REAL& max0,
                         const REAL p[3]) const {
    if( r < 0 ){ return; }
    const REAL L = device_Distance3(p,c);
    if( L < r ){
      min0 = 0;
      max0 = r+L;
      return;
    }
    min0 = L-r;
    max0 = L+r;
  }
public:
  REAL r, c[3];
};

template <typename REAL>
class CudaBV_AABB3
{
public:
  __device__
  void Set_Inactive() {
    bbmin[0] = +1;
    bbmax[0] = -1;
  }
  __device__
  bool IsActive() const {
    if( bbmin[0] > bbmax[0] ){ return false; }
    return true;
  }
  __device__
  void AddPoint(const REAL p[3], REAL eps){
    if( eps < 0 ){ return; }
    if( !this->IsActive() ){ // something inside
      bbmin[0] = p[0]-eps;  bbmax[0] = p[0]+eps;
      bbmin[1] = p[1]-eps;  bbmax[1] = p[1]+eps;
      bbmin[2] = p[2]-eps;  bbmax[2] = p[2]+eps;
      return;
    }
    bbmin[0] = ( bbmin[0] < p[0]-eps ) ? bbmin[0] : p[0]-eps;
    bbmin[1] = ( bbmin[1] < p[1]-eps ) ? bbmin[1] : p[1]-eps;
    bbmin[2] = ( bbmin[2] < p[2]-eps ) ? bbmin[2] : p[2]-eps;
    bbmax[0] = ( bbmax[0] > p[0]+eps ) ? bbmax[0] : p[0]+eps;
    bbmax[1] = ( bbmax[1] > p[1]+eps ) ? bbmax[1] : p[1]+eps;
    bbmax[2] = ( bbmax[2] > p[2]+eps ) ? bbmax[2] : p[2]+eps;
  }
  __device__
  void Add(const CudaBV_AABB3<REAL>& bb){
    if( !bb.IsActive() ){ return; }
    if( !this->IsActive() ){
      bbmax[0] = bb.bbmax[0];	bbmin[0] = bb.bbmin[0];
      bbmax[1] = bb.bbmax[1];	bbmin[1] = bb.bbmin[1];
      bbmax[2] = bb.bbmax[2];	bbmin[2] = bb.bbmin[2];
      return;
    }
    bbmin[0] = ( bbmin[0] < bb.bbmin[0] ) ? bbmin[0] : bb.bbmin[0];
    bbmin[1] = ( bbmin[1] < bb.bbmin[1] ) ? bbmin[1] : bb.bbmin[1];
    bbmin[2] = ( bbmin[2] < bb.bbmin[2] ) ? bbmin[2] : bb.bbmin[2];
    bbmax[0] = ( bbmax[0] > bb.bbmax[0] ) ? bbmax[0] : bb.bbmax[0];
    bbmax[1] = ( bbmax[1] > bb.bbmax[1] ) ? bbmax[1] : bb.bbmax[1];
    bbmax[2] = ( bbmax[2] > bb.bbmax[2] ) ? bbmax[2] : bb.bbmax[2];
    return;
  }
public:
  REAL bbmin[3], bbmax[3];
};

// ---------------------------------------------------------------------------------------------------------------------
// ---------------------------------------------------------------------------------------------------------------------
// ---------------------------------------------------------------------------------------------------------------------

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
  const float l0 = device_Distance3(pc, p0);
  const float l1 = device_Distance3(pc, p1);
  const float l2 = device_Distance3(pc, p2);
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
    const unsigned int inlA = nni+isplit;
    dNodeBVH[ini].ichild[0] = inlA;
    dNodeBVH[inlA].iparent = ini;
    dNodeBVH[inlA].ichild[0] = dSortedId[isplit];
    dNodeBVH[inlA].ichild[1] = UINT_MAX;
  }
  else{
    const unsigned int iniA = isplit;
    dNodeBVH[ini].ichild[0] = iniA;
    dNodeBVH[iniA].iparent = ini;
  }
  // ----
  if( range.y == isplit+1 ){
    const unsigned int inlB = nni+isplit+1;
    dNodeBVH[ini].ichild[1] = inlB;
    dNodeBVH[inlB].iparent = ini;
    dNodeBVH[inlB].ichild[0] = dSortedId[isplit+1];
    dNodeBVH[inlB].ichild[1] = UINT_MAX;
  }
  else{
    const unsigned int iniB = isplit+1;
    dNodeBVH[ini].ichild[1] = iniB;
    dNodeBVH[iniB].iparent = ini;
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
  hNodeBVH[0].iparent = UINT_MAX;
}

// ------------------------------------------------------------------------

template <typename BBOX>
__global__
void kernel_BVHGeometry(
    BBOX* dBox,
    int* dNum,
    //
    const dfm2::CNodeBVH2* dNodeBVH,
    const float* dXYZ,
    const unsigned int* dTri,
    unsigned int nTri,
    float eps)
{
  const unsigned int ino = blockDim.x * blockIdx.x + threadIdx.x;
  if (ino >= nTri) return;

  { // make aabb for triangle
    assert( dNodeBVH[nTri - 1 + ino].ichild[1] == UINT_MAX );
    assert( dNodeBVH[nTri - 1 + ino].iparent < nTri-1 );
    const int itri = dNodeBVH[nTri - 1 + ino].ichild[0];
    assert(itri >= 0 && itri < nTri);
    const unsigned int i0 = dTri[itri * 3 + 0];
    const unsigned int i1 = dTri[itri * 3 + 1];
    const unsigned int i2 = dTri[itri * 3 + 2];
    const float *p0 = dXYZ + i0 * 3;
    const float *p1 = dXYZ + i1 * 3;
    const float *p2 = dXYZ + i2 * 3;
    dBox[nTri-1+ino].Set_Inactive();
    dBox[nTri-1+ino].AddPoint(p0,eps);
    dBox[nTri-1+ino].AddPoint(p1,eps);
    dBox[nTri-1+ino].AddPoint(p2,eps);
  }
  // ----------------------------------------------------
  unsigned int ino0 = dNodeBVH[nTri-1+ino].iparent;
  while(true){
    assert( ino0 < nTri-1 );
    //assert( dNodeBVH[ino0].ichild[0] >= 0 );
    //assert( dNodeBVH[ino0].ichild[1] >= 0 );
    const unsigned int inoc0 = dNodeBVH[ino0].ichild[0];
    const unsigned int inoc1 = dNodeBVH[ino0].ichild[1];
    assert( dNodeBVH[inoc0].iparent == ino0 );
    assert( dNodeBVH[inoc1].iparent == ino0 );
    assert( inoc0 < nTri*2-1 );
    assert( inoc1 < nTri*2-1 );
    const int iflg_old = atomicCAS(dNum+ino0,0,1);
    if( iflg_old == 0 ){ // let the another branch of the binary tree do the work
      return;
    }
    __threadfence(); // sync global memory
    // ---------------------------------------
    dBox[ino0].Set_Inactive();
    dBox[ino0].Add(dBox[inoc0]);
    dBox[ino0].Add(dBox[inoc1]);
    // ----------------------------------------
    if( dNodeBVH[ino0].iparent == UINT_MAX ){ assert(ino0==0); return; }
    ino0 = dNodeBVH[ino0].iparent;
  }
}


void dfm2::cuda::cuda_BVHGeometry_AABB3f(
    dfm2::CBV3_AABB<float>* hAABB,
    const CNodeBVH2* hNodeBVH,
    const float* hXYZ,
    unsigned int nXYZ,
    const unsigned int* hTri,
    unsigned int nTri)
{
  const thrust::device_vector<dfm2::CNodeBVH2> dNodeBVH(hNodeBVH, hNodeBVH+2*nTri-1);
  const thrust::device_vector<float> dXYZ(hXYZ, hXYZ+nXYZ*3);
  const thrust::device_vector<unsigned int> dTri(hTri, hTri+nTri*3);
  thrust::device_vector<CudaBV_AABB3<float>> dAABB( 2*nTri-1 );
  thrust::device_vector<int> dNum(nTri-1, 0);
  // -----------------------------
  {
    const unsigned int BLOCK = 512;
    dim3 grid((nTri - 1) / BLOCK + 1);
    dim3 block(BLOCK);
    kernel_BVHGeometry <<< grid, block >>> (
        thrust::raw_pointer_cast(dAABB.data()),
        thrust::raw_pointer_cast(dNum.data()),
        //
        thrust::raw_pointer_cast(dNodeBVH.data()),
        thrust::raw_pointer_cast(dXYZ.data()),
        thrust::raw_pointer_cast(dTri.data()),
        nTri,
        0.0);
  }
  // ------------------------------
  cudaMemcpy(hAABB,
      thrust::raw_pointer_cast(dAABB.data()),
      sizeof(CudaBV_AABB3<float>)*dAABB.size(),
      cudaMemcpyDeviceToHost);
}

template <typename REAL>
void dfm2::cuda::cuda_BVHGeometry_Sphere(
    dfm2::CBV3_Sphere<REAL>* hAABB,
    const CNodeBVH2* hNodeBVH,
    const REAL* hXYZ,
    unsigned int nXYZ,
    const unsigned int* hTri,
    unsigned int nTri)
{
  const thrust::device_vector<dfm2::CNodeBVH2> dNodeBVH(hNodeBVH, hNodeBVH+2*nTri-1);
  const thrust::device_vector<REAL> dXYZ(hXYZ, hXYZ+nXYZ*3);
  const thrust::device_vector<unsigned int> dTri(hTri, hTri+nTri*3);
  thrust::device_vector<CudaBV_Sphere<REAL>> dAABB( 2*nTri-1 );
  thrust::device_vector<int> dNum(nTri-1, 0);
  // -----------------------------
  {
    const unsigned int BLOCK = 512;
    dim3 grid((nTri - 1) / BLOCK + 1);
    dim3 block(BLOCK);
    kernel_BVHGeometry <<< grid, block >>> (
        thrust::raw_pointer_cast(dAABB.data()),
            thrust::raw_pointer_cast(dNum.data()),
            //
            thrust::raw_pointer_cast(dNodeBVH.data()),
            thrust::raw_pointer_cast(dXYZ.data()),
            thrust::raw_pointer_cast(dTri.data()),
            nTri,
            0.0);
  }
  // ------------------------------
  cudaMemcpy(hAABB,
             thrust::raw_pointer_cast(dAABB.data()),
             sizeof(CudaBV_Sphere<REAL>)*dAABB.size(),
             cudaMemcpyDeviceToHost);
}
template void dfm2::cuda::cuda_BVHGeometry_Sphere(
    dfm2::CBV3_Sphere<float>* hAABB,
    const CNodeBVH2* hNodeBVH,
    const float* hXYZ,
    unsigned int nXYZ,
    const unsigned int* hTri,
    unsigned int nTri);

// -------------------------------------------------------------------------

__device__
void device_BVH_IndPoint_NearestPoint(
    unsigned int* ip,
    float* cur_dist,
    //
    const float p[3],
    unsigned int ibvh,
    const delfem2::CNodeBVH2* aNodeBVH,
    const CudaBV_Sphere<float>* dBVSphere)
{
  float min0=+1.0, max0=-1.0;
  dBVSphere[ibvh].Range_DistToPoint(min0,max0,p);
  //
  if( max0 < min0 ){ return; } // ibvh is a inactive bvh the children should be inactive too
  if( *cur_dist > 0 && min0> *cur_dist ){ return; } // current range [min,max] is valid and nearer than [min0,min0].
  const unsigned int ichild0 = aNodeBVH[ibvh].ichild[0];
  const unsigned int ichild1 = aNodeBVH[ibvh].ichild[1];
  if( ichild1 == UINT_MAX ){ // leaf
    assert( min0 == max0 ); // because this is point
    if( *cur_dist < 0 || max0 < *cur_dist ){ // current range is inactive
      *cur_dist = max0;
      *ip = ichild0;
    }
    return;
  }
  // ------------------
  device_BVH_IndPoint_NearestPoint(
      ip,cur_dist,
      p, ichild0,aNodeBVH,dBVSphere);
  device_BVH_IndPoint_NearestPoint(
      ip,cur_dist,
      p, ichild1,aNodeBVH,dBVSphere);
}

template <typename REAL>
__global__
void kernel_BVHNearestPoint(
    unsigned int* dId,
    //
    const REAL* dXYZ1,
    unsigned int nXYZ1,
    const dfm2::CNodeBVH2* dNodeBVH,
    unsigned int nNodeBVH,
    const CudaBV_Sphere<REAL>* dBVSphere)
{
  const unsigned int ip1 = blockDim.x * blockIdx.x + threadIdx.x;
  if (ip1 >= nXYZ1 ) return;
  const float p1[3] = {dXYZ1[ip1*3+0], dXYZ1[ip1*3+1], dXYZ1[ip1*3+2]};
  float cur_dist = -1;
  device_BVH_IndPoint_NearestPoint(
      dId+ip1, &cur_dist,
      //
      p1,0,dNodeBVH,dBVSphere);
}

template <typename REAL>
void dfm2::cuda::cuda_BVH_NearestPoint(
    unsigned int* hInd,
    //
    const REAL* hXYZ1,
    unsigned int nXYZ1,
    const CNodeBVH2* hNodeBVH0,
    unsigned int nNodeBVH0,
    const dfm2::CBV3_Sphere<REAL>* hBVSphere0)
{
  const thrust::device_vector<REAL> dXYZ1(hXYZ1, hXYZ1+nXYZ1*3);
  const thrust::device_vector<dfm2::CNodeBVH2> dNodeBVH0(hNodeBVH0, hNodeBVH0+nNodeBVH0);
  const thrust::device_vector<CudaBV_Sphere<REAL>> dBVSphere0(nNodeBVH0);
  cudaMemcpy(
      (void*)thrust::raw_pointer_cast(dBVSphere0.data()),
      hBVSphere0,
      sizeof(CudaBV_Sphere<REAL>)*nNodeBVH0,
      cudaMemcpyHostToDevice);
  thrust::device_vector<unsigned int> dInd(nXYZ1, 0);
  // -----------------------------
  {
    cudaDeviceSetLimit(cudaLimitStackSize, 4096);
    const unsigned int BLOCK = 512;
    dim3 grid((nXYZ1 - 1) / BLOCK + 1);
    dim3 block(BLOCK);
    kernel_BVHNearestPoint << < grid, block >> > (
        thrust::raw_pointer_cast(dInd.data()),
        //
        thrust::raw_pointer_cast(dXYZ1.data()),
        nXYZ1,
        thrust::raw_pointer_cast(dNodeBVH0.data()),
        nNodeBVH0,
        thrust::raw_pointer_cast(dBVSphere0.data()));
  }
  thrust::copy(dInd.begin(), dInd.end(), hInd);
}
template void dfm2::cuda::cuda_BVH_NearestPoint(
    unsigned int* hInd,
    //
    const float* hXYZ1,
    unsigned int nXYZ1,
    const CNodeBVH2* hNodeBVH0,
    unsigned int nNodeBVH0,
    const dfm2::CBV3_Sphere<float>* hBVSphere0);
