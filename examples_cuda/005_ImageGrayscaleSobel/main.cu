#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"
#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"

#include "thrust/host_vector.h"
#include "thrust/device_vector.h"

template <typename REAL>
__device__
REAL GrayScale(const uchar3& ptr){
  const float f0 =
      static_cast<float>(ptr.x) * 0.299f
      + static_cast<float>(ptr.y) * 0.587f
      + static_cast<float>(ptr.z) * 0.114f;
  return static_cast<REAL>(f0);
}

__device__
bool IsInside(int ix, int iy, int nx, int ny)
{
  if( ix >=0 && ix < nx && iy >= 0 && iy < ny ){ return true; }
  return false;
}

__global__
void kernel_RGB2Gray(
    char* out,
    const uchar3* in,
    const unsigned int ny,
    const unsigned int nx)
{
  const unsigned int ix = blockIdx.x * blockDim.x + threadIdx.x;
  const unsigned int iy = blockIdx.y * blockDim.y + threadIdx.y;
  if( ix >= nx || iy >= ny ){ return; }
  out[iy*nx+ix] = GrayScale<char>(in[iy*nx+ix]);
}

__global__
void kernel_Sobel(
    float* out_dx,
    float* out_dy,
    const uchar3* in,
    const int ny,
    const int nx)
{
  const int ix1 = int(blockIdx.x * blockDim.x + threadIdx.x);
  const int iy1 = int(blockIdx.y * blockDim.y + threadIdx.y);
  if( ix1 >= nx || iy1 >= ny ){ return; }
  const int ix0 = ::abs(ix1-1); // reflection
  const int iy0 = ::abs(iy1-1); // reflection
  const int ix2 = nx - 1 - ::abs(nx - 2 - ix1); // reflection
  const int iy2 = ny - 1 - ::abs(ny - 2 - iy1); // reflection
  assert( ix0 >=0 && ix0 < nx );
  assert( iy0 >=0 && iy0 < ny );
  assert( ix1 >=0 && ix1 < nx );
  assert( iy1 >=0 && iy1 < ny );
  const auto v00 = GrayScale<float>(in[iy0*nx+ix0]);
  const auto v01 = GrayScale<float>(in[iy0*nx+ix1]);
  const auto v02 = GrayScale<float>(in[iy0*nx+ix2]);
  const auto v10 = GrayScale<float>(in[iy1*nx+ix0]);
//  const auto v11 = GrayScale<float>(in[iy1*nx+ix1]);
  const auto v12 = GrayScale<float>(in[iy1*nx+ix2]);
  const auto v20 = GrayScale<float>(in[iy2*nx+ix0]);
  const auto v21 = GrayScale<float>(in[iy2*nx+ix1]);
  const auto v22 = GrayScale<float>(in[iy2*nx+ix2]);
  const float dx = -v00 - 2*v10 - v20 + v02 + 2*v12 + v22;
  const float dy = -v00 - 2*v01 - v02 + v20 + 2*v21 + v22;
  out_dx[iy1*nx+ix1] = dx;
  out_dy[iy1*nx+ix1] = dy;
}

inline unsigned int iDivUp(
    const unsigned int &a,
    const unsigned int &b )
{
  return ( a%b != 0 ) ? (a/b+1):(a/b);
}

int main() {

  int width, height, channels;
  thrust::device_vector<uchar3> dImgRGB;
  {
    unsigned char *img = stbi_load(
        (std::string(PATH_INPUT_DIR) + "/lenna.png").c_str(),
        &width, &height, &channels, 0);
    const uchar3* ptr0 = reinterpret_cast<uchar3*>(img);
    dImgRGB.assign(ptr0,ptr0+width * height);
    delete[] img;
  }
  const unsigned int blockW = 32;
  const unsigned int blockH = 32;
  const dim3 grid( iDivUp( width, blockW ), iDivUp( height, blockH ) );
  const dim3 threadBlock( blockW, blockH );

  { // convert RGB to Grayscale
    thrust::device_vector<char> dImgGray(width * height);
    kernel_RGB2Gray<<<grid, threadBlock>>>(
        thrust::raw_pointer_cast(dImgGray.data()),
        thrust::raw_pointer_cast(dImgRGB.data()),
        height, width);
    thrust::host_vector<char> hImgGray = dImgGray;
    stbi_write_png(
        (std::string(PATH_SOURCE_DIR) + "/lenna_gray.png").c_str(),
        width, height, 1, hImgGray.data(), sizeof(char) * 0);
  }

  // ------------

  {
    thrust::device_vector<float> dImgDx(width * height);
    thrust::device_vector<float> dImgDy(width * height);
    kernel_Sobel<<<grid, threadBlock>>>(
        thrust::raw_pointer_cast(dImgDx.data()),
        thrust::raw_pointer_cast(dImgDy.data()),
        thrust::raw_pointer_cast(dImgRGB.data()),
        height, width);
    thrust::host_vector<float> hImg;
    thrust::host_vector<unsigned char> aImg(width*height);
    //
    hImg = dImgDx;
    for(int i=0;i<hImg.size();++i) { aImg[i] = (unsigned char)fabs(hImg[i]); }
    stbi_write_png(
        (std::string(PATH_SOURCE_DIR) + "/lenna_dx.png").c_str(),
        width, height, 1, aImg.data(), sizeof(char) * 0);
    //
    hImg = dImgDy;
    for(int i=0;i<hImg.size();++i) { aImg[i] = (unsigned char)fabs(hImg[i]); }
    stbi_write_png(
        (std::string(PATH_SOURCE_DIR) + "/lenna_dy.png").c_str(),
        width, height, 1, aImg.data(), sizeof(char) * 0);
  }
}