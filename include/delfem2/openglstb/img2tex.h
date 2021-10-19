
#ifndef DFM2_OPENGLSTB_IMG2TEX_H
#define DFM2_OPENGLSTB_IMG2TEX_H

/*
 * please define "STB_IMAGE_IMPLEMENTATION" in only one source file (probably in the main.cpp).
 */

#include "delfem2/opengl/tex.h"
#include "stb/stb_image.h"

namespace delfem2::openglstb {

void SetRgbToTex(
    delfem2::opengl::CTexRGB &tex0,
    const std::string &path_img0,
    bool is_flip_vertically) {
  int width, height, channels;
  stbi_set_flip_vertically_on_load(is_flip_vertically);
  unsigned char *img = stbi_load(
      path_img0.c_str(),
      &width, &height, &channels, 0);
  std::cout << width << " " << height << " " << channels << std::endl;
  assert(width > 0 && height > 0);
  tex0.width = width;
  tex0.height = height;
  tex0.channels = channels;
  tex0.pixel_color.assign(img, img + width * height * channels);
  stbi_image_free(img);
}

void LoadImageFileSetToTexture(const char* file_path)
{
  int width, height, channels;
  stbi_set_flip_vertically_on_load(true);
  unsigned char *img = stbi_load(
      file_path,
      &width, &height, &channels, 0);
  assert(width > 0 && height > 0);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);
  glPixelStorei(GL_UNPACK_ALIGNMENT, 1);
  if( channels == 3 ) {
    glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB,
                 static_cast<int>(width),
                 static_cast<int>(height),
                 0, GL_RGB, GL_UNSIGNED_BYTE,
                 img);
  }
  if( channels == 4 ){
    glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB,
                 static_cast<int>(width),
                 static_cast<int>(height),
                 0, GL_RGBA, GL_UNSIGNED_BYTE,
                 img);
  }
  glBindTexture(GL_TEXTURE_2D, 0);
  stbi_image_free(img);
}

}

#endif /* DFM2_OPENGLSTB_IMG2TEX_H */