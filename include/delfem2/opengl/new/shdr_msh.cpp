/*
 * Copyright (c) 2019 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */

#ifdef EMSCRIPTEN
#  include <emscripten/emscripten.h>
#  define GLFW_INCLUDE_ES3
#  include <GLFW/glfw3.h>
#elif defined(USE_GLEW)
#  include <GL/glew.h>
#else
#  include <glad/glad.h>
#endif

#include "delfem2/opengl/funcs.h" // compile shader
#include "delfem2/opengl/new/funcs.h" // CGL4_VAO_Mesh
#include "delfem2/opengl/new/shdr_msh.h"
#include "delfem2/mshmisc.h"

// ------------------------------------------

template<typename REAL>
void delfem2::opengl::CShader_Mesh::Initialize(
    std::vector<REAL> &aXYZd,
    unsigned int ndim,
    std::vector<unsigned int> &aLine,
    int gl_primitive_type) {
  if (!glIsVertexArray(vao.VAO)) { glGenVertexArrays(1, &vao.VAO); }
  vao.Delete_EBOs();
  vao.Add_EBO(aLine, gl_primitive_type);
  this->UpdateVertex(aXYZd, ndim, aLine);
}
#ifdef DFM2_STATIC_LIBRARY
template void delfem2::opengl::CShader_Mesh::Initialize(
    std::vector<float>& aXYZd,
    unsigned int ndim,
    std::vector<unsigned int>& aLine,
    int gl_primitive_type);
template void delfem2::opengl::CShader_Mesh::Initialize(
    std::vector<double>& aXYZd,
    unsigned int ndim,
    std::vector<unsigned int>& aLine,
    int gl_primitive_type);
#endif

template<typename REAL>
void delfem2::opengl::CShader_Mesh::UpdateVertex(
    std::vector<REAL> &aXYZd,
    unsigned int ndim,
    std::vector<unsigned int> &aLine) {
  glBindVertexArray(vao.VAO); // opengl4

  vao.ADD_VBO(0, aXYZd);
  glEnableVertexAttribArray(0);
  glVertexAttribPointer(0, ndim, convertToGlType<REAL>(), GL_FALSE, ndim * sizeof(REAL), (void *) 0); // gl24
}
#ifdef DFM2_STATIC_LIBRARY
template void delfem2::opengl::CShader_Mesh::UpdateVertex(
    std::vector<float>& aXYZd,
    unsigned int ndim,
    std::vector<unsigned int>& aLine);
template void delfem2::opengl::CShader_Mesh::UpdateVertex(
    std::vector<double>& aXYZd,
    unsigned int ndim,
    std::vector<unsigned int>& aLine);
#endif

void delfem2::opengl::CShader_Mesh::Compile() {
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

#ifdef EMSCRIPTEN
  shaderProgram = GL24_CompileShader(
      (std::string("#version 300 es\n")+
      glsl33vert_projection).c_str(),
      (std::string("#version 300 es\n")+
      std::string("precision highp float;\n")+
      glsl33frag).c_str());
#else
  shaderProgram = delfem2::opengl::GL24_CompileShader(
      (std::string("#version 330 core\n") +
      glsl33vert_projection).c_str(),
      (std::string("#version 330 core\n") +
      glsl33frag).c_str());
#endif

  if (!glIsProgram(shaderProgram)) {
    std::cout << "shader doesnot exist" << std::endl;
  }
  glUseProgram(shaderProgram);
  Loc_MatrixProjection = glGetUniformLocation(shaderProgram, "matrixProjection");
  Loc_MatrixModelView = glGetUniformLocation(shaderProgram, "matrixModelView");
  Loc_Color = glGetUniformLocation(shaderProgram, "color");
}

void delfem2::opengl::CShader_Mesh::Draw(float mP[16], float mMV[16]) const {
  glUseProgram(shaderProgram);
  glUniformMatrix4fv(Loc_MatrixProjection, 1, GL_FALSE, mP);
  glUniformMatrix4fv(Loc_MatrixModelView, 1, GL_FALSE, mMV);
  glUniform3f(Loc_Color, color[0], color[1], color[2]);
  vao.Draw(0); // draw line
}
