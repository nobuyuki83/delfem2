#ifndef DFM2_OPENGLSTB_SCREENSHOT_H
#define DFM2_OPENGLSTB_SCREENSHOT_H

#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"
#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"

namespace delfem2 {
namespace openglstb {

void SaveScreenshot(
    const std::string &path) {
  int viewport[4];
  ::glGetIntegerv(GL_VIEWPORT, viewport);
  int width = viewport[2];
  int height = viewport[3];
  GLubyte *pixel_data = (GLubyte *) malloc((width) * (height) * 3 * (sizeof(GLubyte)));
  glPixelStorei(GL_PACK_ALIGNMENT, 1);
  glReadPixels(0, 0,
               width, height,
               GL_RGB,
               GL_UNSIGNED_BYTE,
               pixel_data);
  if (!pixel_data) std::cout << "error pixel data " << std::endl;
  stbi_flip_vertically_on_write(1);
  stbi_write_png(path.c_str(),
                 width, height, 3,
                 pixel_data,
                 0);
  free(pixel_data);
}

}
}

#endif