/**
 * Copyright (c) 2019 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */

/**
 * @details this file must be compile with OpenGL2 and GLEW library.
 */

#ifndef GL2EW_FUNCS_H
#define GL2EW_FUNCS_H

#include <vector>
#include <assert.h>
#include <iostream> // this must be delated in future

// -------------------------------------

class CFrameBufferManager
{
public:
  CFrameBufferManager(){
    id_framebuffer = 0;
    id_depth_render_buffer = 0;
    id_color_render_buffer = 0;
  }
  void DeleteFrameBuffer();
  void Init(int width, int height,
            std::string sFormatPixelColor,
            bool isDepth);
  void Start() const;
  void End() const;
public:
  unsigned int id_framebuffer; // id of this frame buffer
  unsigned int id_depth_render_buffer; // depth buffer
  unsigned int id_color_render_buffer; // color buffer
  int width;
  int height;
  std::string sFormatPixelColor;
};






#endif /* utility_glew_h */

