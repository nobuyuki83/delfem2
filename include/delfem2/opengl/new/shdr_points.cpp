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
void delfem2::opengl::CShader_Points::SetCoords(
    std::vector<REAL> &vtx_coords,
    unsigned int ndim)
{
  if( !::glIsVertexArray(vao.VAO) ){ ::glGenVertexArrays(1, &vao.VAO); }  // opengl ver >= 3.0
  assert(vao.aEBO.empty());
  ::glBindVertexArray(vao.VAO); // opengl ver >= 3.0
  vao.ADD_VBO(0, vtx_coords);
  ::glEnableVertexAttribArray(0);
  ::glVertexAttribPointer(
      0, ndim, convertToGlType<REAL>(), GL_FALSE,
      ndim*sizeof(REAL), (void*)nullptr); // add this for
  num_vtx = vtx_coords.size()/ndim;
}
#ifdef DFM2_STATIC_LIBRARY
template void delfem2::opengl::CShader_Points::SetCoords(std::vector<float>& aXYZd, unsigned int ndim);
template void delfem2::opengl::CShader_Points::SetCoords(std::vector<double>& aXYZd, unsigned int ndim);
#endif

// -------------------------

void delfem2::opengl::CShader_Points::InitGL()
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


void delfem2::opengl::CShader_Points::Draw(
    GLenum gl_primitive_type,
    float mP[16],
    float mMV[16]) const
{
  ::glUseProgram(shaderProgram);
  ::glUniformMatrix4fv(Loc_MatrixProjection, 1, GL_FALSE, mP);
  ::glUniformMatrix4fv(Loc_MatrixModelView, 1, GL_FALSE, mMV);
  ::glUniform3f(Loc_Color, color_face.r,color_face.g, color_face.b);
  ::glBindVertexArray(this->vao.VAO);
  ::glDrawArrays(
      gl_primitive_type, 0,
      static_cast<GLsizei>(num_vtx));
}
