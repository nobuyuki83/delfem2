#define STB_IMAGE_IMPLEMENTATION
#include "stb/stb_image.h"
#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb/stb_image_write.h"
#include <fstream>
#include <iostream>
#include <vector>

namespace delfem2 {
namespace stb {

void ReadImagePair(
    std::vector<float> &aRGB,
    int &width, int &height,
    const std::string &path_l,
    const std::string &path_r) {
  int width_l, height_l, channels_l;
  unsigned char *img_l = stbi_load(
      path_l.c_str(),
      &width_l, &height_l, &channels_l, 0);
  int width_r, height_r, channels_r;
  unsigned char *img_r = stbi_load(
      path_r.c_str(),
      &width_r, &height_r, &channels_r, 0);
  assert(width_l == width_r && height_l == height_r && channels_l == channels_r);
  std::cout << "channel:" << channels_l << " " << channels_r << std::endl;
  width = width_l;
  height = height_l;
  const int N = width * height;
  aRGB.resize(2 * N * 3);
  for (int i = 0; i < N * 3; ++i) { aRGB[i] = float(img_l[i]); }
  for (int i = 0; i < N * 3; ++i) { aRGB[i + N * 3] = float(img_r[i]); }
  delete[] img_l;
  delete[] img_r;
}

void ReadImage(
    std::vector<float> &aRGB,
    int &width, int &height,
    const std::string &path) {
  int channels;
  unsigned char *img = stbi_load(
      path.c_str(),
      &width, &height, &channels, 0);
  const int N = width * height;
  aRGB.resize(N * 3);
  for (int i = 0; i < N * 3; ++i) { aRGB[i] = float(img[i]); }
  delete[] img;
}

void WriteImage(
    const std::vector<float> &aDisp,
    float scale,
    int ny,
    int nx,
    const std::string &path)
{
  assert( aDisp.size() == ny*nx );
  std::vector<unsigned char> aImg(nx * ny);
  for (int iy = 0; iy < ny; ++iy) {
    for (int ix = 0; ix < nx; ++ix) {
      aImg[iy * nx + ix] = char(aDisp[iy * nx + ix] * scale);
    }
  }  
  stbi_write_png(
      path.c_str(),
      nx, ny, 1, aImg.data(), sizeof(unsigned char) * 0);
}


}
}
