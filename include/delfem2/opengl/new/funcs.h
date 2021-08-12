/*
 * Copyright (c) 2019 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */


#ifndef DFM2_OPENGL_NEW_FUNCS_H
#define DFM2_OPENGL_NEW_FUNCS_H

#include "delfem2/dfm2_inline.h"

// -----------
#ifdef EMSCRIPTEN
  #include <GLFW/glfw3.h>
#endif
// -----------
#ifdef _WIN32
  #include <windows.h>
#endif
// -----------
#if defined(__APPLE__) && defined(__MACH__) // Mac
  #include <OpenGL/gl.h>
#else
  #include <GL/gl.h>
#endif
// -----------
#include <vector>
#include <assert.h>
#include <iostream>

namespace delfem2 {
namespace opengl {

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


class CGL4_VAO_Mesh
{
public:
  class CEBO{
  public:
    int GL_MODE;
    unsigned int size;
    int EBO;
  };
  class CVBO{
  public:
    CVBO(){ VBO = 0; }
  public:
    unsigned int VBO;
  };
public:
  CGL4_VAO_Mesh(){
    VAO = 0;
  }

  void Draw(unsigned int iel) const;

  template <typename REAL>
  void ADD_VBO(
      unsigned int ivbo,
      const std::vector<REAL>& aF);

  /*
  void ADD_VBO(
      unsigned int ivbo,
      const std::vector<double>& aD)
  {
    std::vector<float> aF(aD.begin(),aD.end());
    this->ADD_VBO(ivbo, aF);
  }
   */
  void Add_EBO(
      const std::vector<unsigned int>& aTri,
      int GL_MODE);

  void Delete_EBOs(){
    for(size_t ie=0;ie<aEBO.size();++ie){
       unsigned int ebo = aEBO[ie].EBO;
       if( glIsBuffer(ebo) ){ glDeleteBuffers(1,&ebo); }
     }
     aEBO.clear();
  }
public:
  unsigned int VAO;
  std::vector<CEBO> aEBO;
  std::vector<CVBO> aVBO;
};


const std::string glsl33vert_projection =
"uniform mat4 matrixProjection;\n"
"uniform mat4 matrixModelView;\n"
"layout (location = 0) in vec3 posIn;\n"
"layout (location = 1) in vec3 nrmIn;\n"
"out vec3 nrmPrj;\n"
"void main()\n"
"{\n"
"  gl_Position = matrixProjection * matrixModelView * vec4(posIn.x, posIn.y, posIn.z, 1.0);\n"
"  vec4 v0 = matrixModelView * vec4(nrmIn.x, nrmIn.y, nrmIn.z, 0.0);\n"
"  nrmPrj = v0.xyz;\n"
"  if( length(nrmIn) < 1.e-30 ){ nrmPrj = vec3(0.f, 0.f, 1.f); }\n"
"}\0";

const std::string glsl33frag =
"uniform vec3 color;\n"
"in vec3 nrmPrj;\n"
"out vec4 FragColor;\n"
"void main()\n"
"{\n"
"  FragColor = abs(nrmPrj.z)*vec4(color.x, color.y, color.z, 1.0f);\n"
"}\n\0";

}
}

#ifndef DFM2_STATIC_LIBRARY
  #include "delfem2/opengl/new/funcs.cpp"
#endif

#endif /* utility_glew_h */

