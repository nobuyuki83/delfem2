/*
 * Copyright (c) 2019 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */

#include "delfem2/opengl/new/drawer_mshdisp.h"

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
#include "delfem2/mshuni.h"
#include "delfem2/mshmisc.h"

template <typename REAL>
void delfem2::opengl::CShader_TriMesh_Disp::Initialize(
    std::vector<REAL>& aPosD,
    unsigned int ndim,
    std::vector<unsigned int>& aTri,
    std::vector<REAL>& aDispD)
{
  assert( ndim == 2 || ndim == 3 );
  assert( aPosD.size()%ndim == 0 );
  
  std::vector<unsigned int> aLine;
  MeshLine_MeshElem(
      aLine,
      aTri.data(), aTri.size()/3, delfem2::MESHELEM_TRI,
      aPosD.size()/ndim);
  // ------
  if( !glIsVertexArray(vao.idx_vao) ){ glGenVertexArrays(1, &vao.idx_vao); }
  vao.Delete_EBOs();
  vao.Add_EBO(aTri,GL_TRIANGLES);
  vao.Add_EBO(aLine,GL_LINES);
  this->UpdateVertex(aPosD, ndim, aDispD);
}
#ifdef DFM2_STATIC_LIBRARY
template void delfem2::opengl::CShader_TriMesh_Disp::Initialize(
    std::vector<float>& aPosD,
    unsigned int ndim,
    std::vector<unsigned int>& aTri,
    std::vector<float>& aDispD);
template void delfem2::opengl::CShader_TriMesh_Disp::Initialize(
    std::vector<double>& aPosD,
    unsigned int ndim,
    std::vector<unsigned int>& aTri,
    std::vector<double>& aDispD);
#endif

template <typename REAL>
void delfem2::opengl::CShader_TriMesh_Disp::UpdateVertex(
    std::vector<REAL>& aPosD,
    [[maybe_unused]] unsigned int ndim,
    std::vector<REAL>& aDispD)
{
  assert( aPosD.size() == aDispD.size() );
  vao.ADD_VBO(0,aPosD);
  glEnableVertexAttribArray(0);
  glVertexAttribPointer(0, 2, convertToGlType<REAL>(), GL_FALSE, 2*sizeof(REAL), (void*)0); // gl24
  
  vao.ADD_VBO(1,aDispD);
  glEnableVertexAttribArray(1);
  glVertexAttribPointer(1, 2, convertToGlType<REAL>(), GL_FALSE, 2*sizeof(REAL), (void*)0); // gl24
}
#ifdef DFM2_STATIC_LIBRARY
template void delfem2::opengl::CShader_TriMesh_Disp::UpdateVertex(
    std::vector<float>& aPosD,
    unsigned int ndim,
    std::vector<float>& aDispD);
template void delfem2::opengl::CShader_TriMesh_Disp::UpdateVertex(
    std::vector<double>& aPosD,
    unsigned int ndim,
    std::vector<double>& aDispD);
#endif

void delfem2::opengl::CShader_TriMesh_Disp::Compile()
{
  const std::string glsl33vert_projection =
  "uniform mat4 matrixProjection;\n"
  "uniform mat4 matrixModelView;\n"
  "layout (location = 0) in vec3 posIn;\n"
  "layout (location = 1) in vec3 dispIn;\n"
  "void main()\n"
  "{\n"
  "  vec4 pos = vec4(posIn.x, posIn.y, posIn.z, 1.0) + vec4(dispIn.x, dispIn.y, dispIn.z, 1.0);"
  "  gl_Position = matrixProjection * matrixModelView * pos;\n"
  "}\0";
  
  const std::string glsl33frag =
  "uniform vec3 color0;\n"
  "out vec4 FragColor;\n"
  "void main()\n"
  "{\n"
  "  FragColor = vec4(color0.x, color0.y, color0.z, 1.0f);\n"
  "}\n\0";
  
#ifdef EMSCRIPTEN
  shaderProgram = GL24_CompileShader((std::string("#version 300 es\n")+
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
  
  if( !glIsProgram(shaderProgram) ){
    std::cout << "shader doesnot exist" << std::endl;
  }
  glUseProgram(shaderProgram);
  Loc_MatrixProjection = glGetUniformLocation(shaderProgram, "matrixProjection");
  Loc_MatrixModelView  = glGetUniformLocation(shaderProgram, "matrixModelView");
  Loc_Color0           = glGetUniformLocation(shaderProgram, "color0");
  Loc_Color1           = glGetUniformLocation(shaderProgram, "color1");
  Loc_ValMin           = glGetUniformLocation(shaderProgram, "val_min");
  Loc_ValMax           = glGetUniformLocation(shaderProgram, "val_max");
}

void delfem2::opengl::CShader_TriMesh_Disp::Draw(
    float mP[16],
    float mMV[16])
{
  glUseProgram(shaderProgram);
  glUniformMatrix4fv(Loc_MatrixProjection, 1, GL_FALSE, mP);
  glUniformMatrix4fv(Loc_MatrixModelView, 1, GL_FALSE, mMV);
  glUniform3f(Loc_Color0, 0.8, 0.8, 0.8);
  vao.Draw(0); // draw face
  glUniform3f(Loc_Color0, 0,0,0);
  vao.Draw(1); // draw line
}

// ------------------------------------------------

