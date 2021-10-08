/*
 * Copyright (c) 2019 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */

#ifdef EMSCRIPTEN
  #include <emscripten/emscripten.h>
  #define GLFW_INCLUDE_ES3
  #include <GLFW/glfw3.h>
#elif defined(USE_GLEW)
  #include <GL/glew.h>
#else
  #include <glad/glad.h>
#endif
//
#include "delfem2/opengl/funcs.h" // compile shader
#include "delfem2/opengl/new/funcs.h" // CGL4_VAO_Mesh
#include "delfem2/opengl/new/shdr_points.h"
#include "delfem2/mshmisc.h"

// ------------------------------------------

template <typename REAL>
void delfem2::opengl::Drawer_Coords::SetRawArray(
    REAL *vtx_coords,
    size_t num_vtx_,
    unsigned int ndim)
{
  this->num_vtx = num_vtx_;
  if( !::glIsVertexArray(vao.VAO) ){ ::glGenVertexArrays(1, &vao.VAO); }  // opengl ver >= 3.0
  assert(vao.aEBO.empty());
  ::glBindVertexArray(vao.VAO); // opengl ver >= 3.0
  vao.Add_VBO(0, vtx_coords, num_vtx*ndim);
  ::glEnableVertexAttribArray(0);
  ::glVertexAttribPointer(
      0, ndim, convertToGlType<REAL>(), GL_FALSE,
      ndim*sizeof(REAL), (void*)nullptr); // add this for
}
#ifdef DFM2_STATIC_LIBRARY
template void delfem2::opengl::Drawer_Coords::SetRawArray(float*, size_t, unsigned int);
template void delfem2::opengl::Drawer_Coords::SetRawArray(double*, size_t, unsigned int);
#endif

// -------------------------

void delfem2::opengl::Drawer_Coords::InitGL()
{
  const std::string glsl33vert_projection =
      "uniform mat4 matrixProjection;\n"
      "uniform mat4 matrixModelView;\n"
      "layout (location = 0) in vec3 posIn;\n"
      "void main()\n"
      "{\n"
      "  gl_Position = matrixProjection * matrixModelView * vec4(posIn.x, posIn.y, posIn.z, 1.0);\n"
      "}";

  const std::string glsl33frag =
      "uniform vec3 color;\n"
      "out vec4 FragColor;\n"
      "void main()\n"
      "{\n"
      "  FragColor = vec4(color.x, color.y, color.z, 1.0f);\n"
      "}";

#ifdef EMSCRIPTEN
  shaderProgram = GL24_CompileShader(
      (std::string("#version 300 es\n")+
      glsl33vert_projection).c_str(),
      (std::string("#version 300 es\n")+
      std::string("precision highp float;\n")+
      glsl33frag).c_str());
#else
  shaderProgram = delfem2::opengl::GL24_CompileShader(
      (std::string("#version 330 core\n")+
      glsl33vert_projection).c_str(),
      (std::string("#version 330 core\n")+
      glsl33frag).c_str());
#endif

  if( !::glIsProgram(shaderProgram) ){
    std::cerr << "shader doesnot exist" << std::endl;
  }
  ::glUseProgram(shaderProgram);
  Loc_MatrixProjection = ::glGetUniformLocation(shaderProgram,  "matrixProjection");
  Loc_MatrixModelView  = ::glGetUniformLocation(shaderProgram,  "matrixModelView");
  Loc_Color            = ::glGetUniformLocation(shaderProgram,  "color");
}


void delfem2::opengl::Drawer_Coords::Draw(
    GLenum gl_primitive_type,
    const float mP[16],
    const float mMV[16]) const
{
  ::glUseProgram(shaderProgram);
  ::glUniformMatrix4fv(
      Loc_MatrixProjection, 1, GL_FALSE,
      TransposeMat4ForOpenGL(mP,true).data());
  ::glUniformMatrix4fv(
      Loc_MatrixModelView, 1, GL_FALSE,
      TransposeMat4ForOpenGL(mMV, false).data());
  ::glUniform3f(Loc_Color, color_face.r,color_face.g, color_face.b);
  ::glBindVertexArray(this->vao.VAO);
  ::glDrawArrays(
      gl_primitive_type, 0,
      static_cast<GLsizei>(num_vtx));
}
