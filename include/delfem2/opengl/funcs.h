/*
 * Copyright (c) 2019 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */

/**
 * @details the file depends on OpenGL 2.x and 4.x
 * this file need a function loader such as GLEW or GLAD.
 * compile with -USE_GLEW if you want to use GLEW, otherwise it use GLAD
 */

#ifndef DFM2_OPENGL_FUNCS_H
#define DFM2_OPENGL_FUNCS_H

#include <string>
#include <vector>
#include <array>

#include "delfem2/dfm2_inline.h"

namespace delfem2::opengl {

DFM2_INLINE int GL24_CompileShader(
    const char *vert,
    const char* frag);

/**
 *
 * @param str_glsl_vert
 * @param shaderType either GL_VERTEX_SHADER or GL_FRAGMENT_SHADER
 * @return index of shader
 */
DFM2_INLINE int compileShader(
    const std::string& str_glsl_vert,
    int shaderType);

/**
 * @function compile vertex and fragment shader
 * @return id of the shader program
 */
DFM2_INLINE int setUpGLSL(
    const std::string& str_glsl_vert,
    const std::string& str_glsl_frag);

// -------------

template <typename REAL>
constexpr int convertToGlType(){ return 0; }

template <>
constexpr int convertToGlType<float>(){ return GL_FLOAT; }

template <>
constexpr int convertToGlType<double>(){ return GL_DOUBLE; }

// -------------

template <typename T>
std::array<T,16> TransposeMat4ForOpenGL(const T* A, bool is_zflip)
{
  if( is_zflip ){
    return {
        A[0], A[4], -A[8], A[12],
        A[1], A[5], -A[9], A[13],
        A[2], A[6], -A[10], A[14],
        A[3], A[7], -A[11], A[15] };
  }
  else{
    return {
      A[0], A[4], A[8], A[12],
      A[1], A[5], A[9], A[13],
      A[2], A[6], A[10], A[14],
      A[3], A[7], A[11], A[15] };
  }
}

}

#ifndef DFM2_STATIC_LIBRARY
#  include "delfem2/opengl/funcs.cpp"
#endif


#endif  // DFM2_OPENGL_FUNCS_H
