#ifndef DFM2_CUDA_PMSNORMDIST_CUDA_H
#define DFM2_CUDA_PMDNORMDIST_CUDA_H

#include <cassert>
#include "delfem2/cuda/random.h"
#include "delfem2/cvpmsnormdist.h"

namespace delfem2 {
namespace cuda {

__device__
bool device_IsInside(int x, int y, int ubx, int uby) {
  return x >= 0 && x < ubx && y >= 0 && y < uby;
}

__device__
inline float3 device_MatVec3(const float aP[9], const float3 &q) {
  return make_float3(
      aP[0] * q.x + aP[1] * q.y + aP[2] * q.z,
      aP[3] * q.x + aP[4] * q.y + aP[5] * q.z,
      aP[6] * q.x + aP[7] * q.y + aP[8] * q.z);
}

__device__
void device_MatMat3(float X[9], const float A[9], const float B[9]) {
  X[0] = A[0] * B[0] + A[1] * B[3] + A[2] * B[6];
  X[1] = A[0] * B[1] + A[1] * B[4] + A[2] * B[7];
  X[2] = A[0] * B[2] + A[1] * B[5] + A[2] * B[8];
  X[3] = A[3] * B[0] + A[4] * B[3] + A[5] * B[6];
  X[4] = A[3] * B[1] + A[4] * B[4] + A[5] * B[7];
  X[5] = A[3] * B[2] + A[4] * B[5] + A[5] * B[8];
  X[6] = A[6] * B[0] + A[7] * B[3] + A[8] * B[6];
  X[7] = A[6] * B[1] + A[7] * B[4] + A[8] * B[7];
  X[8] = A[6] * B[2] + A[7] * B[5] + A[8] * B[8];
}

__device__
void device_MatTMat3(
    float X[9],
    const float A[9],
    const float B[9]) {
  X[0] = A[0] * B[0] + A[3] * B[3] + A[6] * B[6];
  X[1] = A[0] * B[1] + A[3] * B[4] + A[6] * B[7];
  X[2] = A[0] * B[2] + A[3] * B[5] + A[6] * B[8];
  X[3] = A[1] * B[0] + A[4] * B[3] + A[7] * B[6];
  X[4] = A[1] * B[1] + A[4] * B[4] + A[7] * B[7];
  X[5] = A[1] * B[2] + A[4] * B[5] + A[7] * B[8];
  X[6] = A[2] * B[0] + A[5] * B[3] + A[8] * B[6];
  X[7] = A[2] * B[1] + A[5] * B[4] + A[8] * B[7];
  X[8] = A[2] * B[2] + A[5] * B[5] + A[8] * B[8];
}

__device__
void HomographyPlaneInduced(
    float H[9],
    const float4 &nd, // norm dist
    const CamTransform *camTransf) {
  const float *K0inv = camTransf->K0inv;
  const float *R01 = camTransf->R01;
  const float *t01 = camTransf->t01;
  const float *K1 = camTransf->K1;
  //
  float A[9];
  const float mdinv = -1.f / nd.w;
  A[0] = R01[0] + mdinv * t01[0] * nd.x;
  A[1] = R01[1] + mdinv * t01[0] * nd.y;
  A[2] = R01[2] + mdinv * t01[0] * nd.z;
  A[3] = R01[3] + mdinv * t01[1] * nd.x;
  A[4] = R01[4] + mdinv * t01[1] * nd.y;
  A[5] = R01[5] + mdinv * t01[1] * nd.z;
  A[6] = R01[6] + mdinv * t01[2] * nd.x;
  A[7] = R01[7] + mdinv * t01[2] * nd.y;
  A[8] = R01[8] + mdinv * t01[2] * nd.z;
  float AK0inv[9];
  device_MatMat3(AK0inv, A, K0inv);
  device_MatMat3(H, K1, AK0inv);
}


__device__
void HomographyPlaneInduced2(
    float H[9],
    const float4 &nd, // norm dist
    const CamTransform *camTransf) {
  const float *K0inv = camTransf->K0inv;
  const float *R01 = camTransf->R01;
  const float *t01 = camTransf->t01;
  const float *K1 = camTransf->K1;
  //
  float A[9];
  const float mdinv = +1.f / nd.w;
  A[0] = mdinv * t01[0] * nd.x + 1.f;
  A[1] = mdinv * t01[0] * nd.y;
  A[2] = mdinv * t01[0] * nd.z;
  A[3] = mdinv * t01[1] * nd.x;
  A[4] = mdinv * t01[1] * nd.y + 1.f;
  A[5] = mdinv * t01[1] * nd.z;
  A[6] = mdinv * t01[2] * nd.x;
  A[7] = mdinv * t01[2] * nd.y;
  A[8] = mdinv * t01[2] * nd.z + 1.f;
  float AK0inv[9];
  device_MatMat3(AK0inv, A, K0inv);
  float R01AKinv[9];
  device_MatTMat3(R01AKinv, R01, AK0inv);
  device_MatMat3(H, K1, R01AKinv);
}

__device__
float device_CostPatchMatchStereoNormDist(
    const float4 &nd,
    int ixc, int iyc,
    int nx, int ny, int nrad,
    cudaTextureObject_t texi,
    const CamTransform *camTransf,
    cudaTextureObject_t texj,
    const float params[4]) {
  const float gamma = params[0];
  const float tau_color = params[1];
  const float tau_grad = params[2];
  const float alpha = params[3];

  float H[9];
  HomographyPlaneInduced2(
      H,
      nd, camTransf);

  const auto colicc = tex2D<float>(texi, ixc + 0.5, iyc + 0.5);
  float cost = 0.f;
  for (int ix1 = ixc - nrad; ix1 <= ixc + nrad; ++ix1) {
    for (int iy1 = iyc - nrad; iy1 <= iyc + nrad; ++iy1) {
      if (!device_IsInside(ix1, iy1, nx, ny)) { continue; }
      const auto coli11 = tex2D<float>(texi, ix1 + 0.5, iy1 + 0.5);
      const auto coli01 = tex2D<float>(texi, ix1 - 0.5, iy1 + 0.5);
      const auto coli21 = tex2D<float>(texi, ix1 + 1.5, iy1 + 0.5);
      const auto coli10 = tex2D<float>(texi, ix1 + 0.5, iy1 - 0.5);
      const auto coli12 = tex2D<float>(texi, ix1 + 0.5, iy1 + 1.5);
      const float gxi = (coli21 - coli01) * 0.5f;
      const float gyi = (coli12 - coli10) * 0.5f;
      const float3 pj11 = device_MatVec3(H, make_float3(float(ix1 + 0.5), float(iy1 + 0.5), 1.f));
      const float3 pj01 = device_MatVec3(H, make_float3(float(ix1 - 0.5), float(iy1 + 0.5), 1.f));
      const float3 pj21 = device_MatVec3(H, make_float3(float(ix1 + 1.5), float(iy1 + 0.5), 1.f));
      const float3 pj10 = device_MatVec3(H, make_float3(float(ix1 + 0.5), float(iy1 - 0.5), 1.f));
      const float3 pj12 = device_MatVec3(H, make_float3(float(ix1 + 0.5), float(iy1 + 1.5), 1.f));
      const auto colj11 = tex2D<float>(texj, pj11.x / pj11.z, pj11.y / pj11.z);
      const auto colj01 = tex2D<float>(texj, pj01.x / pj01.z, pj01.y / pj01.z);
      const auto colj21 = tex2D<float>(texj, pj21.x / pj21.z, pj21.y / pj21.z);
      const auto colj10 = tex2D<float>(texj, pj10.x / pj10.z, pj10.y / pj10.z);
      const auto colj12 = tex2D<float>(texj, pj12.x / pj12.z, pj12.y / pj12.z);
      const float gxj = (colj21 - colj01) * 0.5f;
      const float gyj = (colj12 - colj10) * 0.5f;
      const float diffg = fminf(fabsf(gxi - gxj) + fabsf(gyi - gyj), tau_grad);
      const float diffc = fminf(fabsf(coli11 - colj11) * 3.f, tau_color);
      const float diff = (1.f - alpha) * diffc + alpha * diffg;
      const float wi01 = expf(-fabsf(colicc - coli11) * 3.f / gamma);
      cost += wi01 * diff;
    }
  }
  return cost;
}

__global__
void kernel_CostPathMatchNormDist(
    float *aCost,
    const float4 *aNormDist,
    cudaTextureObject_t texi,
    const CamTransform *camTransf,
    cudaTextureObject_t texj,
    int ny,
    int nx,
    int nrad,
    float gamma,
    float tau_col,
    float tau_grad,
    float alpha)
//    const float* aParamCost)
{
  const int ixc = int(blockIdx.x * blockDim.x + threadIdx.x);
  const int iyc = int(blockIdx.y * blockDim.y + threadIdx.y);
  assert(blockIdx.z * blockDim.z + threadIdx.z == 0);
  assert(nrad > 0);
  if (ixc >= nx || iyc >= ny) { return; }
//  printf("%lf, %lf, %lf, %lf\n",aParamCost[0],aParamCost[1],aParamCost[2],aParamCost[3]);
  float aParamCost[4] = {gamma, tau_col, tau_grad, alpha};
  float cost = device_CostPatchMatchStereoNormDist(
      aNormDist[iyc * nx + ixc],
      ixc, iyc, nx, ny, nrad,
      texi, camTransf, texj,
      aParamCost);
  aCost[iyc * nx + ixc] = cost;
}


__global__
void kernel_DiffusionPatchMatchStereoNormDist(
    float *aCost,
    float4 *aNormDist,
    cudaTextureObject_t texi,
    const CamTransform *camTransf,
    cudaTextureObject_t texj,
    int ny,
    int nx,
    int nrad,
    int iflag,
    int ileap,
    float gamma,
    float tau_col,
    float tau_grad,
    float alpha)
{
  const int ixc = int(blockIdx.x * blockDim.x + threadIdx.x);
  int itmp = int(blockIdx.y * blockDim.y + threadIdx.y);
  assert(blockIdx.z * blockDim.z + threadIdx.z == 0);
  const int iyc = (ixc % 2 == 0) ? itmp * 2 + (1 - iflag) : itmp * 2 + iflag;
  if (ixc >= nx || iyc >= ny) { return; }
  assert(nrad > 0);
  const float aParamCost[4] = {gamma, tau_col, tau_grad, alpha};
  int aLeap[4][2] = {
      {-ileap, 0},
      {0,      -ileap},
      {+ileap, 0},
      {0,      +ileap}};
  for (auto &inb : aLeap) {
    const int ix1 = ixc + inb[0];
    const int iy1 = iyc + inb[1];
    if (!device_IsInside(ix1, iy1, nx, ny)) { continue; }
    float cost_new = device_CostPatchMatchStereoNormDist(
        aNormDist[iy1 * nx + ix1],
        ixc, iyc, nx, ny, nrad,
        texi, camTransf, texj,
        aParamCost);
//    printf("%d %d %d %d --> %lf %lf\n",ixc,iyc,ix1,iy1,aCost[iyc * nx + ixc],cost_new);
    if (cost_new < aCost[iyc * nx + ixc]) {
      aCost[iyc * nx + ixc] = cost_new;
      aNormDist[iyc * nx + ixc] = aNormDist[iy1 * nx + ix1];
    }
  }
}

__global__
void kernel_RefinePlanePatchMatchStereoNormDist(
    float *aCost,
    float4 *aNormDist,
    cudaTextureObject_t texi,
    const CamTransform *camTransf,
    cudaTextureObject_t texj,
    int ny,
    int nx,
    int nrad,
    float gamma,
    float tau_col,
    float tau_grad,
    float alpha,
    int iseedrand)
{
  const uint ix0 = blockIdx.x * blockDim.x + threadIdx.x;
  const uint iy0 = blockIdx.y * blockDim.y + threadIdx.y;
  if (ix0 >= nx || iy0 >= ny) { return; }
  assert(blockIdx.z * blockDim.z + threadIdx.z == 0);
  const uint i0 = iy0 * nx + ix0;
  uint4 seed = {
      i0 * 536523269 + 93655483 + 70001611 * iseedrand,
      i0 * 536529611 + 93660977 + 70001881 * iseedrand,
      i0 * 536527793 + 93655207 + 70003321 * iseedrand,
      i0 * 536522887 + 93652997 + 70001011 * iseedrand};
  const float aParamCost[4] = {gamma, tau_col, tau_grad, alpha};
  float mag_dn = 0.1;
  float mag_dz = 0.1;
  for (int itr = 0; itr < 10; ++itr) {
    const float4 nd0 = aNormDist[iy0 * nx + ix0];
    float4 nd1;
    {
      const float3 x0 = {float(ix0+0.5), float(iy0+0.5), 1.f};
      const float3 K0invx0 = device_MatVec3(camTransf->K0inv, x0);
      const float n0tK0invx0 = K0invx0.x * nd0.x + K0invx0.y * nd0.y + K0invx0.z * nd0.z;
      const float z0 = nd0.w / n0tK0invx0;
      seed = delfem2::cuda::device_Xorshift128(seed);
      const float nx1 = nd0.x + delfem2::cuda::device_RandomRange(seed.w, -mag_dn, +mag_dn);
      seed = delfem2::cuda::device_Xorshift128(seed);
      const float ny1 = nd0.y + delfem2::cuda::device_RandomRange(seed.w, -mag_dn, +mag_dn);
      seed = delfem2::cuda::device_Xorshift128(seed);
      const float nz1 = nd0.z + delfem2::cuda::device_RandomRange(seed.w, -mag_dn, +mag_dn);
      seed = delfem2::cuda::device_Xorshift128(seed);
      const float z1 = z0 + delfem2::cuda::device_RandomRange(seed.w, -mag_dz, +mag_dz);
      const float linv = 1.f / sqrtf(nx1 * nx1 + ny1 * ny1 + nz1 * nz1);
      nd1.x = nx1 * linv;
      nd1.y = ny1 * linv;
      nd1.z = nz1 * linv;
      const float n1tK0invx0 = K0invx0.x * nd1.x + K0invx0.y * nd1.y + K0invx0.z * nd1.z;
      nd1.w = z1 * n1tK0invx0;
//      printf("%d %d %d --> %lf,  %lf\n",iseedrand,ix0,iy0,z0,z1);
//      printf("%lf %lf %lf %lf\n",nd1.x,nd1.y,nd1.z,nd1.w);
    }
    const float cost_new = device_CostPatchMatchStereoNormDist(
        nd1,
        ix0, iy0, nx, ny, nrad,
        texi, camTransf, texj,
        aParamCost);
    if (cost_new < aCost[iy0 * nx + ix0]*0.6) {
      aNormDist[iy0 * nx + ix0] = nd1;
      aCost[iy0 * nx + ix0] = cost_new;
    }
    mag_dz *= 0.5f;
    mag_dn *= 0.5f;
  }
}

}
}

#endif