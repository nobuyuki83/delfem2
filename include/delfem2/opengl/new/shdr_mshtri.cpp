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
#include "delfem2/opengl/new/shdr_mshtri.h"
#include "delfem2/mshuni.h"
#include "delfem2/mshmisc.h"


namespace dfm2 = delfem2;

// ------------------------------------------

template <typename REAL>
void delfem2::opengl::CShader_TriMesh::Initialize(
    std::vector<REAL>& aXYZd,
    unsigned int ndim,
    std::vector<unsigned int>& aTri)
{
  std::vector<unsigned int> aLine;
  MeshLine_MeshElem(
      aLine,
      aTri.data(), aTri.size()/3, dfm2::MESHELEM_TRI,
      aXYZd.size()/ndim);
  // --------
  if( !glIsVertexArray(vao.VAO) ){ glGenVertexArrays(1, &vao.VAO); }
  vao.Delete_EBOs();
  vao.Add_EBO(aTri,GL_TRIANGLES);
  vao.Add_EBO(aLine,GL_LINES);
  this->UpdateVertex(aXYZd, ndim, aTri);
}
#ifdef DFM2_STATIC_LIBRARY
template void delfem2::opengl::CShader_TriMesh::Initialize(
    std::vector<float>& aXYZd,
    unsigned int ndim,
    std::vector<unsigned int>& aTri);
template void delfem2::opengl::CShader_TriMesh::Initialize(
    std::vector<double>& aXYZd,
    unsigned int ndim,
    std::vector<unsigned int>& aTri);
#endif

template <typename REAL>
void delfem2::opengl::CShader_TriMesh::UpdateVertex(
    std::vector<REAL>& aXYZd,
    unsigned int ndim,
    std::vector<unsigned int>& aTri)
{
  std::vector<REAL> aNrmd(aXYZd.size());
  delfem2::Normal_MeshTri3D(
      aNrmd.data(),
      aXYZd.data(), aXYZd.size()/ndim,
      aTri.data(), aTri.size()/3);

  glBindVertexArray(vao.VAO); // opengl4
  
  vao.ADD_VBO(0,aXYZd);
  glEnableVertexAttribArray(0);
  glVertexAttribPointer(0, ndim, convertToGlType<REAL>(), GL_FALSE, ndim*sizeof(REAL), (void*)0); // gl24
  
  vao.ADD_VBO(1,aNrmd);
  glEnableVertexAttribArray(1);
  glVertexAttribPointer(1, ndim, convertToGlType<REAL>(), GL_FALSE, ndim*sizeof(REAL), (void*)0); // gl24
}
#ifdef DFM2_STATIC_LIBRARY
template void delfem2::opengl::CShader_TriMesh::UpdateVertex(
    std::vector<float>& aXYZd,
    unsigned int ndim,
    std::vector<unsigned int>& aTri);
template void delfem2::opengl::CShader_TriMesh::UpdateVertex(
    std::vector<double>& aXYZd,
    unsigned int ndim,
    std::vector<unsigned int>& aTri);
#endif


void delfem2::opengl::CShader_TriMesh::Compile()
{
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
  "  nrmPrj = normalize(nrmPrj);\n"
  "}\0";
  
  const std::string glsl33frag =
  "uniform vec3 color;\n"
  "in vec3 nrmPrj;\n"
  "out vec4 FragColor;\n"
  "void main()\n"
  "{\n"
  "  FragColor = abs(nrmPrj.z)*vec4(color.x, color.y, color.z, 1.0f);\n"
//  "  FragColor = vec4(color.x, color.y, color.z, 1.0f);\n"
//  "  FragColor = vec4(1.f, 0.f, 0.f, 1.0f);\n"
  "}\n\0";
  
#ifdef EMSCRIPTEN
  shaderProgram = GL24_CompileShader((std::string("#version 300 es\n")+
                                      glsl33vert_projection).c_str(),
                                     (std::string("#version 300 es\n")+
                                      std::string("precision highp float;\n")+
                                      glsl33frag).c_str());
#else
  shaderProgram = dfm2::opengl::GL24_CompileShader((std::string("#version 330 core\n")+
                                                    glsl33vert_projection).c_str(),
                                                   (std::string("#version 330 core\n")+
                                                    glsl33frag).c_str());
#endif
  
  if( !glIsProgram(shaderProgram) ){
    std::cout << "shader doesnot exist" << std::endl;
  }
  glUseProgram(shaderProgram);
  Loc_MatrixProjection = glGetUniformLocation(shaderProgram,  "matrixProjection");
  Loc_MatrixModelView  = glGetUniformLocation(shaderProgram,  "matrixModelView");
  Loc_Color            = glGetUniformLocation(shaderProgram,  "color");
}


void delfem2::opengl::CShader_TriMesh::Draw(float mP[16], float mMV[16]) const
{
  glUseProgram(shaderProgram);
  glUniformMatrix4fv(Loc_MatrixProjection, 1, GL_FALSE, mP);
  glUniformMatrix4fv(Loc_MatrixModelView, 1, GL_FALSE, mMV);
  glUniform3f(Loc_Color, this->color_face.r , this->color_face.g, this->color_face.b);
  glLineWidth(this->line_width);
  vao.Draw(0); // draw face
  glUniform3f(Loc_Color, 0,0,0);
  vao.Draw(1); // draw line
}
