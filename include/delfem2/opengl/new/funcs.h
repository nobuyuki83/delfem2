/*
 * Copyright (c) 2019 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */


#ifndef DFM2_OPENGL_NEW_FUNCS_H
#define DFM2_OPENGL_NEW_FUNCS_H

#include "delfem2/dfm2_inline.h"

#include <vector>
#include <cassert>
#include <iostream>
// -----------
#ifdef _WIN32
  #include <windows.h>
#endif
// -----------
#ifdef EMSCRIPTEN
  #include <emscripten/emscripten.h>
  #define GLFW_INCLUDE_ES3
#else
#include <glad/glad.h>
#endif
// -----------
#if defined(__APPLE__) && defined(__MACH__) // Mac
  #include <OpenGL/gl.h>
#else
  #include <GL/gl.h>
#endif
// ----------

namespace delfem2::opengl {

/**
 * @details the OpenGL ES 2.0 only accept float array. So there is no "double" version of this file
 */
void GL4_VAO_Pos(
    unsigned int& VAO,
    unsigned int& VBO,
    const float* aP, int nP, int nDim);

void GL4_VAO_PosNrm(
    unsigned int& VAO,
    unsigned int& VBO_pos,
    unsigned int& VBO_nrm,
    const float* aP, int nP, int nDim,
    const float* aN);


class VertexArrayObject
{
public:
  class CEBO{
  public:
    int GL_MODE;
    size_t size;
    unsigned int EBO;
  };
  class CVBO{
  public:
    CVBO(){ VBO = 0; }
  public:
    unsigned int VBO;
  };
public:
  VertexArrayObject(){
    VAO = 0;
  }

  void Draw(unsigned int iel) const;

  template <typename REAL>
  void Add_VBO(
      unsigned int idx_vbo,
      const REAL *vtx_coords,
      size_t num_dof);

  template <typename REAL>
  void ADD_VBO(
      unsigned int idx_vbo,
      const std::vector<REAL>& vtx_coords){
    this->Add_VBO(idx_vbo, vtx_coords.data(), vtx_coords.size());
  }

  void Add_EBO(
      const std::vector<unsigned int>& elem_vtx,
      int gl_primitive_type);

  void Delete_EBOs();

  void DrawArray(
    int gl_primitive_type,
    unsigned int np) const;
  
public:
  unsigned int VAO;
  std::vector<CEBO> aEBO;
  std::vector<CVBO> aVBO;
};

}

#ifndef DFM2_STATIC_LIBRARY
  #include "delfem2/opengl/new/funcs.cpp"
#endif

#endif // DFM2_OPENGL_NEW_FUNCS_H

