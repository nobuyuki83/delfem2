//
// Created by Nobuyuki Umetani on 2021-09-03.
//

#ifndef DFM2_OPENGLSTB_SCREENSHOT_H
#define DFM2_OPENGLSTB_SCREENSHOT_H

#include "stb/stb_image_write.h"

namespace delfem2::openglstb {

void TakeScreenShot(
    int width, int height,
    const std::string &fpath) {
  GLubyte *pixel_data = (GLubyte *) malloc((width) * (height) * 3 * (sizeof(GLubyte)));
  glPixelStorei(GL_PACK_ALIGNMENT, 1);
  glReadPixels(0, 0,
               width, height,
               GL_RGB,
               GL_UNSIGNED_BYTE,
               pixel_data);
  if (!pixel_data) {
    std::cout << "error pixel data " << std::endl;
  }

  stbi_flip_vertically_on_write(1);
  stbi_write_png(fpath.c_str(),
                 width, height, 3,
                 pixel_data,
                 0);
  free(pixel_data);
}

} // namespace delfem2::openstb

#endif // DFM2_OPENGLSTB_SCREENSHOT_H
