/*
 * Copyright (c) 2019 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */

#ifdef USE_GLEW
  #include <GL/glew.h>
#else
  #include <glad/glad.h>
  #ifdef EMSCRIPTEN
    #include <emscripten/emscripten.h>
    #define GLFW_INCLUDE_ES3
  #endif
#endif

#include "delfem2/mshtopo.h"
#include "delfem2/msh.h"
#include "delfem2/color.h"

#include "delfem2/gl24_funcs.h" // compile shader
#include "delfem2/gl4_funcs.h" // CGL4_VAO_Mesh
#include "delfem2/gl4_mshcolor.h"


void CShader_TriMesh::Initialize
 (std::vector<double>& aXYZd,
  std::vector<unsigned int>& aTri)
{
  std::vector<unsigned int> aLine;
  MeshLine_MeshElem(aLine,
                    aTri.data(), aTri.size()/3, MESHELEM_TRI,
                    aXYZd.size()/3);
  // --------
  if( !glIsVertexArray(vao.VAO) ){ glGenVertexArrays(1, &vao.VAO); }
  vao.Delete_EBOs();
  vao.Add_EBO(aTri,GL_TRIANGLES);
  vao.Add_EBO(aLine,GL_LINES);
  this->UpdateVertex(aXYZd, aTri);
}

void CShader_TriMesh::UpdateVertex
 (std::vector<double>& aXYZd,
  std::vector<unsigned int>& aTri)
{
  std::vector<double> aNrmd(aXYZd.size());
  Normal_MeshTri3D(aNrmd.data(),
                   aXYZd.data(), aXYZd.size()/3,
                   aTri.data(), aTri.size()/3);

  glBindVertexArray(vao.VAO); // opengl4
  
  vao.ADD_VBO(0,aXYZd);
  glEnableVertexAttribArray(0);
  glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 3*sizeof(float), (void*)0); // gl24
  
  vao.ADD_VBO(1,aNrmd);
  glEnableVertexAttribArray(1);
  glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, 3*sizeof(float), (void*)0); // gl24
}


void CShader_TriMesh::Compile()
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
  shaderProgram = GL24_CompileShader((std::string("#version 300 es\n")+
                                      glsl33vert_projection).c_str(),
                                     (std::string("#version 300 es\n")+
                                      std::string("precision highp float;\n")+
                                      glsl33frag).c_str());
#else
  shaderProgram = GL24_CompileShader((std::string("#version 330 core\n")+
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


void CShader_TriMesh::Draw(float mP[16], float mMV[16])
{
  glUseProgram(shaderProgram);
  glUniformMatrix4fv(Loc_MatrixProjection, 1, GL_FALSE, mP);
  glUniformMatrix4fv(Loc_MatrixModelView, 1, GL_FALSE, mMV);
  glUniform3f(Loc_Color, 1,0,0);
  vao.Draw(0);
  glUniform3f(Loc_Color, 0,0,0);
  vao.Draw(1);
}


// ----------------------------------------------------------------

void CShader_TriMesh_Scalar::Initialize
(std::vector<double>& aPosD,
 unsigned int ndim,
 std::vector<unsigned int>& aTri,
 std::vector<double>& aValD)
{
  assert( ndim == 2 || ndim == 3 );
  assert( aPosD.size()%ndim == 0 );
  
  std::vector<unsigned int> aLine;
  MeshLine_MeshElem(aLine,
                    aTri.data(), aTri.size()/3, MESHELEM_TRI,
                    aPosD.size()/ndim);
  // ------
  if( !glIsVertexArray(vao.VAO) ){ glGenVertexArrays(1, &vao.VAO); }
  vao.Delete_EBOs();
  vao.Add_EBO(aTri,GL_TRIANGLES);
  vao.Add_EBO(aLine,GL_LINES);
  this->UpdateVertex(aPosD, ndim, aValD);
}

void CShader_TriMesh_Scalar::UpdateVertex
(std::vector<double>& aPosD,
 unsigned int ndim,
 std::vector<double>& aValD)
{
  vao.ADD_VBO(0,aPosD);
  glEnableVertexAttribArray(0);
  glVertexAttribPointer(0, 2, GL_FLOAT, GL_FALSE, 2*sizeof(float), (void*)0); // gl24
  
  vao.ADD_VBO(1,aValD);
  glEnableVertexAttribArray(1);
  glVertexAttribPointer(1, 1, GL_FLOAT, GL_FALSE, 1*sizeof(float), (void*)0); // gl24
}


void CShader_TriMesh_Scalar::Compile()
{
  const std::string glsl33vert_projection =
  "uniform mat4 matrixProjection;\n"
  "uniform mat4 matrixModelView;\n"
  "layout (location = 0) in vec3 posIn;\n"
  "layout (location = 1) in float valIn;\n"
  "out float val;\n"
  "void main()\n"
  "{\n"
  "  gl_Position = matrixProjection * matrixModelView * vec4(posIn.x, posIn.y, posIn.z, 1.0);\n"
  "  val = valIn;\n"
  "}\0";
  
  const std::string glsl33frag =
  "uniform vec3 color0;\n"
  "uniform vec3 color1;\n"
  "uniform float val_min;\n"
  "uniform float val_max;\n"
  "in float val;\n"
  "out vec4 FragColor;\n"
  "void main()\n"
  "{\n"
  "  float r01 = (val-val_min)/(val_max-val_min);\n"
  "  vec3 clr01 = (1.f-r01)*color0 + r01*color1;\n"
  "  FragColor = vec4(clr01.x, clr01.y, clr01.z, 1.0f);\n"
  "}\n\0";
  
#ifdef EMSCRIPTEN
  shaderProgram = GL24_CompileShader((std::string("#version 300 es\n")+
                                      glsl33vert_projection).c_str(),
                                     (std::string("#version 300 es\n")+
                                      std::string("precision highp float;\n")+
                                      glsl33frag).c_str());
#else
  shaderProgram = GL24_CompileShader((std::string("#version 330 core\n")+
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


void CShader_TriMesh_Scalar::Draw(float mP[16], float mMV[16])
{
  glUseProgram(shaderProgram);
  glUniformMatrix4fv(Loc_MatrixProjection, 1, GL_FALSE, mP);
  glUniformMatrix4fv(Loc_MatrixModelView, 1, GL_FALSE, mMV);
  glUniform1f(Loc_ValMin, val_min);
  glUniform1f(Loc_ValMax, val_max);
  glUniform3f(Loc_Color0, color_min.r,color_min.g,color_min.b);
  glUniform3f(Loc_Color1, color_max.r,color_max.g,color_max.b);
  vao.Draw(0);
  glUniform3f(Loc_Color0, 0,0,0);
  glUniform3f(Loc_Color1, 0,0,0);
  vao.Draw(1);
}

