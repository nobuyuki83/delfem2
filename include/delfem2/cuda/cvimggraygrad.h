#ifndef DFM2_CUDA_CVIMGGRAYGRAD_H
#define DFM2_CUDA_CVIMGGRAYGRAD_H

#include <cassert>
#include <stdio.h>

namespace delfem2 {
namespace cuda {

__device__
float device_Dot(
    const float3 &p,
    const float3 &q) {
  return p.x * q.x + p.y * q.y + p.z * q.z;
}

__device__
void device_ImageGrayGradient(
    float2 *aGrad_,
    const float3 *aRgb_,
    int iv1,
    int iy1,
    int ix1,
    int nv,
    int ny,
    int nx,
    const float3 &gray_coeff) { // make grad
  const delfem2::cuda::ViewTensor3Const<float3> aRgb(aRgb_, nv, ny, nx);
  delfem2::cuda::ViewTensor3<float2> aGrad(aGrad_, nv, ny, nx);
  const int ix0 = abs(ix1 - 1); // reflectionCost(iv1, iy1, ix1)
  const int iy0 = abs(iy1 - 1); // reflection
  const int ix2 = nx - 1 - abs(nx - 2 - ix1); // reflection
  const int iy2 = ny - 1 - abs(ny - 2 - iy1); // reflection
  assert(ix0 >= 0 && ix0 < nx);
  assert(iy0 >= 0 && iy0 < ny);
  assert(ix1 >= 0 && ix1 < nx);
  assert(iy1 >= 0 && iy1 < ny);
  assert(ix2 >= 0 && ix2 < nx);
  assert(iy2 >= 0 && iy2 < ny);
  const auto v00 = device_Dot(aRgb(iv1, iy0, ix0), gray_coeff);
  const auto v01 = device_Dot(aRgb(iv1, iy0, ix1), gray_coeff);
  const auto v02 = device_Dot(aRgb(iv1, iy0, ix2), gray_coeff);
  const auto v10 = device_Dot(aRgb(iv1, iy1, ix0), gray_coeff);
  const auto v12 = device_Dot(aRgb(iv1, iy1, ix2), gray_coeff);
  const auto v20 = device_Dot(aRgb(iv1, iy2, ix0), gray_coeff);
  const auto v21 = device_Dot(aRgb(iv1, iy2, ix1), gray_coeff);
  const auto v22 = device_Dot(aRgb(iv1, iy2, ix2), gray_coeff);
  const float dx = (-v00 - 2.f * v10 - v20 + v02 + 2.f * v12 + v22) / 8.f;
  const float dy = (-v00 - 2.f * v01 - v02 + v20 + 2.f * v21 + v22) / 8.f;
  aGrad(iv1, iy1, ix1).x = dx;
  aGrad(iv1, iy1, ix1).y = dy;
}

// ----------------------------------------------------

__global__
void kernel_ImageGrayGradient(
    float2* aGrad_,
    const float3* aRgb_,
    const int nv,
    const int ny,
    const int nx,
    const float3 gray_coeff)
{
  const int ix1 = int(blockIdx.x * blockDim.x + threadIdx.x);
  const int iy1 = int(blockIdx.y * blockDim.y + threadIdx.y);
  const int iv1 = int(blockIdx.z * blockDim.z + threadIdx.z);
  if (ix1 >= nx || iy1 >= ny || iv1 >= nv ) { return; }
  device_ImageGrayGradient(
      aGrad_,
      aRgb_,
      iv1, iy1, ix1,
      nv, ny, nx,
      gray_coeff);
}


}
}

#endif
