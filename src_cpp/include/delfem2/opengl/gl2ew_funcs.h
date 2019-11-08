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


class CElemBuffObj{
public:
  void SetBuffer_Elem(const std::vector<unsigned int>& aTri, unsigned int gl_elem_type);
  void DrawBuffer() const ;
public:
  unsigned int iebo;
  unsigned int gl_elem_type;
  unsigned int size_elem;
  bool is_lighting;
};

class CGLBuffer
{
public:
  CGLBuffer(){
    vbo = -1;
    vbo_nrm = -1;
  }
  void SetBuffer_Vtx(const std::vector<double>& aXYZ, int ndim);
  void SetBuffer_Nrm(const std::vector<double>& aNrm);
  void Draw_Start() const;
  void Draw_End() const ;
public:
  
public:
  unsigned int vbo;
  unsigned int vbo_nrm;
  unsigned int ndim;
};




#endif /* utility_glew_h */

