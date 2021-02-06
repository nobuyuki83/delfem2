

#include "delfem2/opengl/tex.h"
#define STB_IMAGE_IMPLEMENTATION
#include "stb/stb_image.h"

namespace delfem2 {
namespace openglstb {

void SetRgbToTex(
    delfem2::opengl::CTexRGB &tex0,
    const std::string &path_img0) {
  int width, height, channels;
  unsigned char *img = stbi_load(
      path_img0.c_str(),
      &width, &height, &channels, 0);
  assert(width > 0 && height > 0 && channels == 3);
  tex0.w = width;
  tex0.h = height;
  tex0.aRGB.assign(img, img + width * height * channels);
  delete[] img;
}

}
}