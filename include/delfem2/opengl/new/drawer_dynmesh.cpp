/*
 * Copyright (c) 2019 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */

#include "delfem2/opengl/new/drawer_dynmesh.h"

#include <cstdio>
#ifdef EMSCRIPTEN
  #include <emscripten/emscripten.h>
  #define GLFW_INCLUDE_ES3
#else
  #include <glad/glad.h>
#endif
#include <GLFW/glfw3.h>

#include "delfem2/opengl/new/funcs.h"
#include "delfem2/opengl/funcs.h"
#include "delfem2/cad2_dtri2.h"

namespace dfm2 = delfem2;

// ------------------------------------------------------------------------------


void dfm2::opengl::CShader_MeshDTri2D::MakeBuffer (
  const std::vector<CVec2d>& aVec2,
  const std::vector<CDynTri>& aETri)
{
  std::vector<float> aXYf;
  std::vector<unsigned int> aTri;
  std::vector<unsigned int> aLine;
  {
    aXYf.resize(aVec2.size()*2);
    for(size_t iv=0;iv<aVec2.size();++iv){
      aXYf[iv*2+0] = aVec2[iv].x;
      aXYf[iv*2+1] = aVec2[iv].y;
    }
    aTri.resize(aETri.size()*3);
    for(size_t it=0;it<aETri.size();++it){
      aTri[it*3+0] = aETri[it].v[0];
      aTri[it*3+1] = aETri[it].v[1];
      aTri[it*3+2] = aETri[it].v[2];
    }
    aLine.reserve(aTri.size()*1.5*1.1);
    for(const auto & it : aETri){
      int i0 = it.v[0];
      int i1 = it.v[1];
      int i2 = it.v[2];
      if( it.s2[0] == UINT_MAX || i2 > i1 ){ aLine.push_back(i1); aLine.push_back(i2); }
      if( it.s2[1] == UINT_MAX || i0 > i2 ){ aLine.push_back(i2); aLine.push_back(i0); }
      if( it.s2[2] == UINT_MAX || i1 > i0 ){ aLine.push_back(i0); aLine.push_back(i1); }
    }
  }
  // ----------
  if( !glIsVertexArray(vao.VAO) ){ glGenVertexArrays(1, &vao.VAO); }
  glBindVertexArray(vao.VAO);

  vao.ADD_VBO(0,aXYf);
  glEnableVertexAttribArray(0);
  glVertexAttribPointer(0, 2, GL_FLOAT, GL_FALSE, 2*sizeof(float), (void*)0); // gl24

  vao.Delete_EBOs();
  {
    unsigned int EBO_Tri;
    glGenBuffers(1, &EBO_Tri);
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, EBO_Tri);
    glBufferData(GL_ELEMENT_ARRAY_BUFFER, sizeof(unsigned int)*aTri.size(), aTri.data(), GL_STATIC_DRAW);
    VertexArrayObject::CEBO e0;
    e0.size = aTri.size();
    e0.GL_MODE = GL_TRIANGLES;
    e0.EBO = EBO_Tri;
    vao.aEBO.push_back(e0);
  }
  {
    unsigned int EBO_Line;
    glGenBuffers(1, &EBO_Line);
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, EBO_Line);
    glBufferData(GL_ELEMENT_ARRAY_BUFFER, sizeof(unsigned int)*aLine.size(), aLine.data(), GL_STATIC_DRAW);
    VertexArrayObject::CEBO e0;
    e0.size = aLine.size();
    e0.GL_MODE = GL_LINES;
    e0.EBO = EBO_Line;
    vao.aEBO.push_back(e0);
  }
}

void dfm2::opengl::CShader_MeshDTri2D::Compile()
{
  const std::string glsl33vert_projection =
  "uniform mat4 matrixProjection;\n"
  "uniform mat4 matrixModelView;\n"
  "layout (location = 0) in vec2 posIn;\n"
  "void main()\n"
  "{\n"
  "  gl_Position = matrixProjection * matrixModelView * vec4(posIn.x, posIn.y, 0.0, 1.0);\n"
  "}";

  const std::string glsl33frag =
  "uniform vec3 color;\n"
  "out vec4 FragColor;\n"
  "void main()\n"
  "{\n"
  "  FragColor = vec4(color.x, color.y, color.z, 1.0f);\n"
  "}";
  
  #ifdef EMSCRIPTEN
    shdr0_program = GL24_CompileShader((std::string("#version 300 es\n")+
                                        glsl33vert_projection).c_str(),
                                       (std::string("#version 300 es\n")+
                                        std::string("precision highp float;\n")+
                                        glsl33frag).c_str());
  #else
    shdr0_program = GL24_CompileShader((std::string("#version 330 core\n")+
                                        glsl33vert_projection).c_str(),
                                       (std::string("#version 330 core\n")+
                                        glsl33frag).c_str());
  #endif
    
    
  assert( glIsProgram(shdr0_program) );
  glUseProgram(shdr0_program);
  shdr0_Loc_MatrixProjection = glGetUniformLocation(shdr0_program,  "matrixProjection");
  shdr0_Loc_MatrixModelView  = glGetUniformLocation(shdr0_program,  "matrixModelView");
  shdr0_Loc_Color            = glGetUniformLocation(shdr0_program,  "color");
  std::cout << "  shaderProgram: " << shdr0_program;
  std::cout << "  projectionMatrixLoc: " << shdr0_Loc_MatrixProjection;
  std::cout << "  LocColor: " << shdr0_Loc_Color << std::endl;
}


void dfm2::opengl::CShader_MeshDTri2D::Draw
 (const float mP[16],
  const float mMV[16]) const
{
  glUseProgram(shdr0_program);
  glUniformMatrix4fv(shdr0_Loc_MatrixProjection, 1, GL_FALSE, mP);
  glUniformMatrix4fv(shdr0_Loc_MatrixModelView, 1, GL_FALSE, mMV);
  glUniform3f(shdr0_Loc_Color, 1,1,1);
  vao.Draw(0);
  glUniform3f(shdr0_Loc_Color, 0,0,0);
  vao.Draw(1);
}
