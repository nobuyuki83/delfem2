
#include "thrust/host_vector.h"
#include "thrust/device_vector.h"

namespace delfem2 {
namespace cuda {

namespace tex {

__device__
float device_Dot(
    const float3 &p,
    const float3 &q) {
  return p.x * q.x + p.y * q.y + p.z * q.z;
}

__global__
void kernel_RGB2Gray(
    float *aGray_,
    const float3 *aRgb_,
    int ny,
    int nx,
    int pitch,
    const float3 gray_coeff) {
  const int ix1 = int(blockIdx.x * blockDim.x + threadIdx.x);
  const int iy1 = int(blockIdx.y * blockDim.y + threadIdx.y);
  assert(blockIdx.z * blockDim.z + threadIdx.z == 0);
  if (ix1 >= nx || iy1 >= ny) { return; }
  aGray_[iy1 * pitch + ix1] = ::delfem2::cuda::tex::device_Dot(aRgb_[iy1 * nx + ix1], gray_coeff);
}

inline unsigned int iDivUp(
    const unsigned int &a,
    const unsigned int &b) {
  return (a % b != 0) ? (a / b + 1) : (a / b);
}

} // tex

// ---------------------------------------------------------

__global__
void kernel_CopyFromTexture(
    float* aGray_,
    const cudaTextureObject_t texGray,
    int ny,
    int nx)
{
  const int ix1 = int(blockIdx.x * blockDim.x + threadIdx.x);
  const int iy1 = int(blockIdx.y * blockDim.y + threadIdx.y);
  const int iz1 = int(blockIdx.z * blockDim.z + threadIdx.z);
  if( ix1 >= nx || iy1 >= ny || iz1 >= 1 ){ return; }
  aGray_[iy1*nx+ix1] = tex2D<float>(texGray,ix1+0.5,iy1+0.5);
}

__global__
void kernel_GradTexture(
    float* aGrad_,
    const cudaTextureObject_t texGray,
    int ny,
    int nx)
{
  const int ix1 = int(blockIdx.x * blockDim.x + threadIdx.x);
  const int iy1 = int(blockIdx.y * blockDim.y + threadIdx.y);
  assert( blockIdx.z * blockDim.z + threadIdx.z == 0 );
  if( ix1 >= nx || iy1 >= ny ){ return; }
  const auto v21 = tex2D<float>(texGray,ix1+0.5 + 1.0,iy1+0.5);
  const auto v01 = tex2D<float>(texGray,ix1+0.5 - 1.0,iy1+0.5);
  const auto v12 = tex2D<float>(texGray,ix1+0.5,iy1+0.5 + 1.0);
  const auto v10 = tex2D<float>(texGray,ix1+0.5,iy1+0.5 - 1.0);
  aGrad_[(iy1*nx+ix1)*2+0] = (v21 - v01)*0.5f;
  aGrad_[(iy1*nx+ix1)*2+1] = (v12 - v10)*0.5f;
}

float *MakeGrayTextureFromRgb(
    cudaTextureObject_t &tex,
    int ny,
    int nx,
    const float *paRgb,
    const float gray_coeff[3]) {
  float *pdBuffer;
  size_t pitchInByte;
  cudaMallocPitch(&pdBuffer, &pitchInByte, sizeof(float) * nx, ny);

  { // from RGB to gray
    thrust::device_vector<float3> dRgb;
    {
      const auto *ptr = reinterpret_cast<const float3 *>(paRgb);
      dRgb.assign(ptr, ptr + ny * nx);
    }
    const unsigned int blockW = 32;
    const unsigned int blockH = 32;
    const dim3 grid(tex::iDivUp(nx, blockW), tex::iDivUp(ny, blockH), 1);
    const dim3 threadBlock(blockW, blockH, 1);
    tex::kernel_RGB2Gray<<< grid, threadBlock >>>(
        pdBuffer,
        thrust::raw_pointer_cast(dRgb.data()),
        ny, nx, pitchInByte / sizeof(float),
        make_float3(gray_coeff[0], gray_coeff[1], gray_coeff[2]));
    cudaDeviceSynchronize();
  }

  struct cudaResourceDesc resDesc{};
  memset(&resDesc, 0, sizeof(resDesc));
  resDesc.res.pitch2D.desc = cudaCreateChannelDesc(
      32, 0, 0, 0, cudaChannelFormatKindFloat);
  resDesc.res.pitch2D.pitchInBytes = pitchInByte;
  resDesc.res.pitch2D.width = nx;
  resDesc.res.pitch2D.height = ny;
  resDesc.res.pitch2D.devPtr = pdBuffer;// thrust::raw_pointer_cast(dGray.data());
  resDesc.resType = cudaResourceTypePitch2D;

  struct cudaTextureDesc texDesc{};
  memset(&texDesc, 0, sizeof(texDesc));
  texDesc.addressMode[0] = cudaAddressModeClamp;
  texDesc.addressMode[1] = cudaAddressModeClamp;
  // Wrap and Mirror is for Normalized coords
//  texDesc.addressMode[0] = cudaAddressModeWrap;
//  texDesc.addressMode[1] = cudaAddressModeWrap;
//  texDesc.addressMode[0] = cudaAddressModeMirror;
//  texDesc.addressMode[1] = cudaAddressModeMirror;
  texDesc.filterMode = cudaFilterModeLinear;
  texDesc.readMode = cudaReadModeElementType;
  texDesc.normalizedCoords = 0;
//  texDesc.sRGB = 0;

  {
    cudaError_t error = cudaCreateTextureObject(&tex, &resDesc, &texDesc, nullptr);
    if (error != cudaSuccess) { std::cout << cudaGetErrorString(error) << std::endl; }
  }
  return pdBuffer;
}


}
}
