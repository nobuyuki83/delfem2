#ifndef DFM2_CUDA_PMSLINDISP_CUDA_H
#define DFM2_CUDA_PMSLINDISP_CUDA_H

#include "delfem2/cuda/viewtensor.h"
#include "delfem2/cuda/random.h"
#include "delfem2/cuda/cvimggraygrad.h"
#include <cassert>

namespace delfem2 {
namespace cuda {

__device__
bool device_IsInside(int x, int y, int ubx, int uby) {
  return x >= 0 && x < ubx && y >= 0 && y < uby;
}

__device__
float3 device_PlaneCoeff(float3 pnt, float3 nrm) {
  float3 res;
  res.x = -nrm.x / nrm.z;
  res.y = -nrm.y / nrm.z;
  res.z = (pnt.x * nrm.x + pnt.y * nrm.y + pnt.z * nrm.z) / nrm.z;
  return res;
}

__device__
float3 device_Normalize3(float3 p) {
  const float leninv = 1.f / sqrtf(p.x * p.x + p.y * p.y + p.z * p.z);
  return make_float3(
      p.x * leninv,
      p.y * leninv,
      p.z * leninv);
}

__device__
float3 device_PlaneNormal(
    const float3 coeff) {
  float3 nrm;
  nrm.z = 1.f / sqrtf(1.f + coeff.x * coeff.x + coeff.y * coeff.y);
  nrm.x = -coeff.x * nrm.z;
  nrm.y = -coeff.y * nrm.z;
  return nrm;
}


__device__
float2 device_InterpFloat2(
    float2 p,
    float2 q,
    float ratio) {
  return make_float2(
      (1.f - ratio) * p.x + ratio * q.x,
      (1.f - ratio) * p.y + ratio * q.y);
}

__device__
float3 device_InterpFloat3(
    float3 p,
    float3 q,
    float ratio) {
  return make_float3(
      (1.f - ratio) * p.x + ratio * q.x,
      (1.f - ratio) * p.y + ratio * q.y,
      (1.f - ratio) * p.z + ratio * q.z);
}

__device__
inline float device_Dissimilarity(
    const float3 pp,
    const float3 qq,
    const float2 pg,
    const float2 qg,
    float alpha,
    float tau_c,
    float tau_g) {
  float cost_c = fabs(pp.x - qq.x) + fabs(pp.y - qq.y) + fabs(pp.z - qq.z);
  float cost_g = fabs(pg.x - qg.x) + fabs(pg.y - qg.y);
  cost_c = min(cost_c, tau_c);
  cost_g = min(cost_g, tau_g);
  return (1.f - alpha) * cost_c + alpha * cost_g;
}

__device__
float device_CostPlaneMatch(
    float3 plane,
    int iv0,
    int iy0,
    int ix0,
    int ny,
    int nx,
    int ws, // window size
    const float3 *aRgb_,
    const float2 *aGrad_,
    const float *aWeight_,
    float alpha,
    float tau_c,
    float tau_g,
    float MAX_DISPARITY,
    float PLANE_PENALTY) {
  const delfem2::cuda::ViewTensor3Const<float3> aRgb(aRgb_, 2, ny, nx);
  const delfem2::cuda::ViewTensor3Const<float2> aGrad(aGrad_, 2, ny, nx);
  const delfem2::cuda::ViewTensor5Const<float> aWeight(aWeight_, 2, ny, nx, ws, ws);
  const int half = ws / 2;
  assert(ws == half * 2 + 1);
  const float disparity_orientation_sign = -1.f + 2.f * float(iv0);
  float cost = 0.f;
  for (int jx = ix0 - half; jx <= ix0 + half; ++jx) {
    for (int jy = iy0 - half; jy <= iy0 + half; ++jy) {
      if (!device_IsInside(jx, jy, nx, ny)) { continue; }
      const float dsp = plane.x * float(jx) + plane.y * float(jy) + plane.z;
      if (dsp < 0 || dsp > MAX_DISPARITY) {
        cost += PLANE_PENALTY;
        continue;
      }
      const float fx_match = float(jx) + disparity_orientation_sign * dsp;
      int kx_match = int(floor(fx_match));
      if (kx_match > nx - 2) { kx_match = nx - 2; }
      if (kx_match < 0) { kx_match = 0; }
      const float ratio = fx_match - float(kx_match);
      assert(kx_match >= 0 && kx_match < nx);
      assert(kx_match + 1 >= 0 && kx_match + 1 < nx);
      const float3 c_k0 = aRgb(1 - iv0, jy, kx_match);
      const float3 c_k1 = aRgb(1 - iv0, jy, kx_match + 1);
      const float3 c_k = device_InterpFloat3(c_k0, c_k1, ratio);
      const float2 g_k0 = aGrad(1 - iv0, jy, kx_match);
      const float2 g_k1 = aGrad(1 - iv0, jy, kx_match + 1);
      const float2 g_k = device_InterpFloat2(g_k0, g_k1, ratio);
      const float weight_ij = aWeight(iv0, iy0, ix0, jy - iy0 + half, jx - ix0 + half);
      const float dissim_ijk = device_Dissimilarity(
          aRgb(iv0, jy, jx), c_k,
          aGrad(iv0, jy, jx), g_k,
          alpha, tau_c, tau_g);
      cost += weight_ij * dissim_ijk;
    }
  }
  return cost;
}


__device__
void device_RefinePlane(
    float3 *aPlane_,
    float *aCost_,
    const float2 *aGrad_,
    const float *aWeight_,
    const float3 *aRgb_,
    int iv0,
    int iy0,
    int ix0,
    int ny,
    int nx,
    int ws,
    float max_delta_z,
    float max_delta_n,
    float end_dz,
    float alpha,
    float tau_c,
    float tau_g,
    float MAX_DISPARITY,
    float PLANE_PENALTY,
    uint4 seed) {
  delfem2::cuda::ViewTensor3<float> aCost(aCost_, 2, ny, nx);
  delfem2::cuda::ViewTensor3<float3> aPlane(aPlane_, 2, ny, nx);
  const delfem2::cuda::ViewTensor3Const<float2> aGrad(aGrad_, 2, ny, nx);
  const delfem2::cuda::ViewTensor3Const<float3> aRgb(aRgb_, 2, ny, nx);
  const delfem2::cuda::ViewTensor5Const<float> aWeight(aWeight_, 2, ny, nx, ws, ws);
  float max_dz = max_delta_z;
  float max_dn = max_delta_n;
  while (max_dz >= end_dz) {
    const float old_cost = aCost(iv0, iy0, ix0);

    float3 new_plane;
    {
      seed = delfem2::cuda::device_Xorshift128(seed);
      const float ddisp = delfem2::cuda::device_RandomRange(seed.w, -max_dz, +max_dz);
      seed = delfem2::cuda::device_Xorshift128(seed);
      const float dnx = delfem2::cuda::device_RandomRange(seed.w, -max_dn, +max_dn);
      seed = delfem2::cuda::device_Xorshift128(seed);
      const float dny = delfem2::cuda::device_RandomRange(seed.w, -max_dn, +max_dn);
      seed = delfem2::cuda::device_Xorshift128(seed);
      const float dnz = delfem2::cuda::device_RandomRange(seed.w, -max_dn, +max_dn);
      const float3 old_plane = aPlane(iv0, iy0, ix0);
      const float disp_old = old_plane.x * float(ix0) + old_plane.y * float(iy0) + old_plane.z;
      const float3 pnt_new = make_float3(float(ix0), float(iy0), disp_old + ddisp);
      const float3 nrm_old = device_PlaneNormal(old_plane);
      float3 nrm_new = {
          nrm_old.x + dnx,
          nrm_old.y + dny,
          nrm_old.z + dnz};
      nrm_new = device_Normalize3(nrm_new);
      new_plane = device_PlaneCoeff(pnt_new, nrm_new);
    }

    const float new_cost = device_CostPlaneMatch(
        new_plane,
        iv0, iy0, ix0,
        ny, nx, ws,
        aRgb_, aGrad_, aWeight_,
        alpha, tau_c, tau_g, MAX_DISPARITY, PLANE_PENALTY);

//    printf("update %f %f\n",new_cost, old_cost);

    if (new_cost < old_cost) {
      aPlane(iv0, iy0, ix0) = new_plane;
      aCost(iv0, iy0, ix0) = new_cost;
    }

    max_dz /= 2.0f;
    max_dn /= 2.0f;
  }
}

__device__
void device_ViewTransform(
    float3 &pl1,
    int &iqy,
    int &iqx,
    int iv0,
    int iy,
    int ix,
    const float3 pl0) {
  const float sign_disp = (iv0 == 0) ? -1.f : 1.f;
  const float disp = pl0.x * float(ix) + pl0.y * float(iy) + pl0.z;
  iqx = int(float(ix) + sign_disp * disp);
  iqy = iy;
  const float3 p = {float(iqx), float(iqy), disp};
  const float3 nrm = device_PlaneNormal(pl0);
  pl1 = device_PlaneCoeff(p, nrm);
}

__device__
void device_PropagationView(
    float3 *aPlane_,
    float *aCost_,
    const float2 *aGrad_,
    const float *aWeight_,
    const float3 *aRgb_,
    int iv0,
    int iy0,
    int ix0,
    const int ny,
    const int nx,
    const int ws,
    float alpha,
    float tau_c,
    float tau_g,
    float MAX_DISPARITY,
    float PLANE_PENALTY) {
  delfem2::cuda::ViewTensor3<float> aCost(aCost_, 2, ny, nx);
  delfem2::cuda::ViewTensor3<float3> aPlane(aPlane_, 2, ny, nx);
  const float3 plane_i = aPlane(iv0, iy0, ix0);

  int kx0, ky0;
  float3 plane_new;
  {
    device_ViewTransform(
        plane_new, ky0, kx0,
        iv0, iy0, ix0, plane_i);
  }
  if (!device_IsInside(kx0, ky0, nx, ny)) { return; }

  const float cost_new = device_CostPlaneMatch(
      plane_new,
      1 - iv0, ky0, kx0, ny, nx, ws,
      aRgb_, aGrad_, aWeight_,
      alpha, tau_c, tau_g, MAX_DISPARITY, PLANE_PENALTY);

  if (cost_new < aCost(1 - iv0, ky0, kx0)) {
//    aPlane(1 - iv0, my, mx).x = new_plane.x;
//    aPlane(1 - iv0, my, mx).y = new_plane.y;
//    aPlane(1 - iv0, my, mx).z = new_plane.z;
//    aCost(1 - iv0, my, mx) = new_cost;
    atomicExch(&aPlane(1 - iv0, ky0, kx0).x, plane_new.x);
    atomicExch(&aPlane(1 - iv0, ky0, kx0).y, plane_new.y);
    atomicExch(&aPlane(1 - iv0, ky0, kx0).z, plane_new.z);
    atomicExch(&aCost(1 - iv0, ky0, kx0), cost_new);
  }

}

__device__
void device_PropagationSpatial(
    float3 *aPlane_,
    float *aCost_,
    int iv0,
    int iy0,
    int ix0,
    int stride,
    const float2 *aGrad_,
    const float *aWeight_,
    const float3 *aRgb_,
    int ny,
    int nx,
    int ws,
    float alpha,
    float tau_c,
    float tau_g,
    float MAX_DISPARITY,
    float PLANE_PENALTY) {
  delfem2::cuda::ViewTensor3<float> aCost(aCost_, 2, ny, nx);
  delfem2::cuda::ViewTensor3<float3> aPlane(aPlane_, 2, ny, nx);
  const int offsets[4][2] = {
      {-stride, 0},
      {0,       -stride},
      {+stride, 0},
      {0,       +stride}};
  for (auto offset : offsets) {
    const int iy_new = iy0 + offset[0];
    const int ix_new = ix0 + offset[1];
    if (!device_IsInside(ix_new, iy_new, nx, ny)) { continue; }
    const float3 plane_new = aPlane(iv0, iy_new, ix_new);
    const float cost_new = device_CostPlaneMatch(
        plane_new,
        iv0, iy0, ix0, ny, nx, ws,
        aRgb_, aGrad_, aWeight_,
        alpha, tau_c, tau_g, MAX_DISPARITY, PLANE_PENALTY);
    if (cost_new >= aCost(iv0, iy0, ix0)) { continue; }
    aPlane(iv0, iy0, ix0).x = plane_new.x;
    aPlane(iv0, iy0, ix0).y = plane_new.y;
    aPlane(iv0, iy0, ix0).z = plane_new.z;
    aCost(iv0, iy0, ix0) = cost_new;
  }
}

// above: device functions
// ----------------------------------------------------------------------------
// below: kernel functions

__global__
void kernel_Preprocess(
    float3 *aPlane_,
    float2 *aGrad_,
    float *aWeight_,
    const float3 *aRgb_,
    const int ny,
    const int nx,
    const int ws,
    float gamma,
    float MAX_DISPARITY,
    const float3 gray_coeff) {
  const int ix1 = int(blockIdx.x * blockDim.x + threadIdx.x);
  const int iy1 = int(blockIdx.y * blockDim.y + threadIdx.y);
  const int iv1 = int(blockIdx.z * blockDim.z + threadIdx.z);
  if (ix1 >= nx || iy1 >= ny) { return; }
  assert(iv1 < 2);

  delfem2::cuda::device_ImageGrayGradient(
      aGrad_,
      aRgb_, iv1, iy1, ix1, 2, ny, nx,
      gray_coeff);

  const delfem2::cuda::ViewTensor3Const<float3> aRgb(aRgb_, 2, ny, nx);
  {
    delfem2::cuda::ViewTensor5<float> aWeight(aWeight_, 2, ny, nx, ws, ws);
    const int hws = ws / 2;
    assert(ws == hws * 2 + 1);
    for (int jx = ix1 - hws; jx <= ix1 + hws; ++jx) {
      for (int jy = iy1 - hws; jy <= iy1 + hws; ++jy) {
        if (!device_IsInside(jx, jy, nx, ny)) { continue; }
        const float3 &p = aRgb(iv1, iy1, ix1);
        const float3 &q = aRgb(iv1, jy, jx);
        const float diffcolor = fabs(q.x - p.x) + fabs(q.y - p.y) + fabs(q.z - p.z);
        const float w0 = expf(-diffcolor / gamma);
        aWeight(iv1, iy1, ix1, jy - iy1 + hws, jx - ix1 + hws) = w0;
      }
    }
  }
  {
    delfem2::cuda::ViewTensor3<float3> aPlane(aPlane_, 2, ny, nx);
    const uint ix0 = blockIdx.x * blockDim.x + threadIdx.x;
    const uint iy0 = blockIdx.y * blockDim.y + threadIdx.y;
    const uint iv0 = blockIdx.z * blockDim.z + threadIdx.z;
    const uint i0 = iv0 * nx * ny + iy0 * nx + ix0;
    uint4 seed = {
        i0 * 536523133 + 836523119,
        i0 * 336529153 + 836530489,
        i0 * 536526257 + 836525551,
        i0 * 336522877 + 836522521};
    seed = delfem2::cuda::device_Xorshift128(seed);
    const float nrm_x = delfem2::cuda::device_RandomRange(seed.w, -1.f, 1.f);
    seed = delfem2::cuda::device_Xorshift128(seed);
    const float nrm_y = delfem2::cuda::device_RandomRange(seed.w, -1.f, 1.f);
    seed = delfem2::cuda::device_Xorshift128(seed);
    const float nrm_z = delfem2::cuda::device_RandomRange(seed.w, -1.f, 1.f);
    seed = delfem2::cuda::device_Xorshift128(seed);
    const float pnt_z = delfem2::cuda::device_RandomRange(seed.w, 0.f, MAX_DISPARITY);
    float3 nrm0 = {nrm_x, nrm_y, nrm_z};
    nrm0 = device_Normalize3(nrm0);
    aPlane(iv1, iy1, ix1) = device_PlaneCoeff(
        {(float) ix1, (float) iy1, pnt_z},
        nrm0);
  }
}


__global__
void kernel_PatchMatch(
    float3 *aPlane_,
    float *aCost_,
    int iseedrand,
    const float2 *aGrad_,
    const float *aWeight_,
    const float3 *aRgb_,
    const int ny,
    const int nx,
    const int ws,
    float alpha,
    float tau_c,
    float tau_g,
    float MAX_DISPARITY,
    float PLANE_PENALTY) {
  const int ix1 = int(blockIdx.x * blockDim.x + threadIdx.x);
  const int iy1 = int(blockIdx.y * blockDim.y + threadIdx.y);
  const int iv1 = int(blockIdx.z * blockDim.z + threadIdx.z);
  if (ix1 >= nx || iy1 >= ny) { return; }
  assert(iv1 < 2);
  //
  const int aStride[5] = {8, 4, 2, 1};
  for (int stride : aStride) {
    if (ix1 % stride == 0 && iy1 % stride == 0) {
      device_PropagationSpatial(
          aPlane_, aCost_,
          iv1, iy1, ix1, stride,
          aGrad_, aWeight_, aRgb_,
          ny, nx, ws,
          alpha, tau_c, tau_g, MAX_DISPARITY, PLANE_PENALTY);
    }
    __syncthreads();
  }

  uint4 seed;
  {
    const uint ix0 = blockIdx.x * blockDim.x + threadIdx.x;
    const uint iy0 = blockIdx.y * blockDim.y + threadIdx.y;
    const uint iv0 = blockIdx.z * blockDim.z + threadIdx.z;
    const uint i0 = iv0 * nx * ny + iy0 * nx + ix0;
    seed = {
        i0 * 536523269 + 93655483 + 70001611 * iseedrand,
        i0 * 536529611 + 93660977 + 70001881 * iseedrand,
        i0 * 536527793 + 93655207 + 70003321 * iseedrand,
        i0 * 536522887 + 93652997 + 70001011 * iseedrand};
  }
  device_RefinePlane(
      aPlane_, aCost_,
      aGrad_, aWeight_, aRgb_,
      iv1, iy1, ix1,
      ny, nx, ws,
      MAX_DISPARITY * 0.5f, 1.0f, 0.1f,
      alpha, tau_c, tau_g, MAX_DISPARITY, PLANE_PENALTY,
      seed);
  device_PropagationView(
      aPlane_, aCost_,
      aGrad_, aWeight_, aRgb_,
      iv1, iy1, ix1,
      ny, nx, ws,
      alpha, tau_c, tau_g, MAX_DISPARITY, PLANE_PENALTY);
}

__global__
void kernel_Cost(
    float *aCost_,
    const float3 *aPlane_,
    const float2 *aGrad_,
    const float *aWeight_,
    const float3 *aRgb_,
    const int ny,
    const int nx,
    const int ws,
    float alpha,
    float tau_c,
    float tau_g,
    float MAX_DISPARITY,
    float PLANE_PENALTY) {
  const int ix1 = int(blockIdx.x * blockDim.x + threadIdx.x);
  const int iy1 = int(blockIdx.y * blockDim.y + threadIdx.y);
  const int iv1 = int(blockIdx.z * blockDim.z + threadIdx.z);
  if (ix1 >= nx || iy1 >= ny) { return; }
  assert(iv1 < 2);
  delfem2::cuda::ViewTensor3<float> aCost(aCost_, 2, ny, nx);
  const delfem2::cuda::ViewTensor3Const<float3> aPlane(aPlane_, 2, ny, nx);
  aCost(iv1, iy1, ix1) = device_CostPlaneMatch(
      aPlane(iv1, iy1, ix1),
      iv1, iy1, ix1, ny, nx, ws,
      aRgb_, aGrad_, aWeight_,
      alpha, tau_c, tau_g, MAX_DISPARITY, PLANE_PENALTY);
}

}
}

#endif