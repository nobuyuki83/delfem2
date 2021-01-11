#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"
#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"

#include "thrust/host_vector.h"
#include "thrust/device_vector.h"


__global__
void kernel(
    char* out,
    const uchar3* in,
    const unsigned int ny,
    const unsigned int nx)
{
  const unsigned int ix = blockIdx.x * blockDim.x + threadIdx.x;
  const unsigned int iy = blockIdx.y * blockDim.y + threadIdx.y;
  const float outf =
         static_cast<float>(in[iy*nx+ix].x) * 0.299f
       + static_cast<float>(in[iy*nx+ix].y) * 0.587f
       + static_cast<float>(in[iy*nx+ix].z) * 0.114f;
  out[iy*nx+ix] = static_cast<char>(outf);
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
    uchar3* ptr0 = reinterpret_cast<uchar3*>(img);
    dImgRGB.assign(ptr0,ptr0+width * height);
    delete[] img;
  }
  thrust::device_vector<char> dImgGray(width*height);

  const unsigned int blockW = 32;
  const unsigned int blockH = 32;
  const dim3 grid( iDivUp( width, blockW ), iDivUp( height, blockH ) );
  const dim3 threadBlock( blockW, blockH );

  kernel<<<grid,threadBlock>>>(
      thrust::raw_pointer_cast(dImgGray.data()),
      thrust::raw_pointer_cast(dImgRGB.data()),
      height, width);

  thrust::host_vector<char> hImgGray = dImgGray;
  stbi_write_png(
      (std::string(PATH_SOURCE_DIR)+"/lenna_gray.png").c_str(),
      width, height, 1, hImgGray.data(), sizeof(char)*0);
}