/*
 * Copyright (c) 2019 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */

#ifdef EMSCRIPTEN
  #include <emscripten/emscripten.h>
  #define GLFW_INCLUDE_ES3
#else
  #include <glad/glad.h>
#endif
#include "delfem2/opengl/new/funcs.h"
#include "delfem2/opengl/new/v23dtricad.h"
#include "delfem2/opengl/funcs.h"
#include "delfem2/cad2_dtri2.h"
#include <GLFW/glfw3.h>
#include <stdio.h>

namespace dfm2 = delfem2;

// -----------------------------------

void AddLine
(std::vector<float>& aPxyNxyf,
 std::vector< std::vector<unsigned int>>& aaTri,
 const std::vector<float>& aXY0f,
 const std::vector<unsigned int>& aLine)
{
  const size_t npl = aLine.size();
  const size_t icnt_p0 = aPxyNxyf.size()/8;
  aPxyNxyf.resize(aPxyNxyf.size()+npl*8,0.0);
  for(unsigned int ipl=0;ipl<npl;++ipl){
    aPxyNxyf[(icnt_p0+ipl)*8+0] = aXY0f[aLine[ipl]*2+0];
    aPxyNxyf[(icnt_p0+ipl)*8+1] = aXY0f[aLine[ipl]*2+1];
    aPxyNxyf[(icnt_p0+ipl)*8+2] = 0.0;
    aPxyNxyf[(icnt_p0+ipl)*8+3] = 0.0;
    aPxyNxyf[(icnt_p0+ipl)*8+4] = aXY0f[aLine[ipl]*2+0];
    aPxyNxyf[(icnt_p0+ipl)*8+5] = aXY0f[aLine[ipl]*2+1];
    aPxyNxyf[(icnt_p0+ipl)*8+6] = 0.0;
    aPxyNxyf[(icnt_p0+ipl)*8+7] = 0.0;
  }
  const size_t nseg = aLine.size()-1;
  for(size_t iseg=0;iseg<nseg;++iseg){
    const int i0 = aLine[iseg+0];
    const int i1 = aLine[iseg+1];
    double v01x = aXY0f[i1*2+0] - aXY0f[i0*2+0];
    double v01y = aXY0f[i1*2+1] - aXY0f[i0*2+1];
    double len = sqrt(v01x*v01x+v01y*v01y);
    double n01x = +v01y/len;
    double n01y = -v01x/len;
    aPxyNxyf[(icnt_p0+iseg+0)*8+2] += +n01x;
    aPxyNxyf[(icnt_p0+iseg+0)*8+3] += +n01y;
    aPxyNxyf[(icnt_p0+iseg+0)*8+6] += -n01x;
    aPxyNxyf[(icnt_p0+iseg+0)*8+7] += -n01y;
    aPxyNxyf[(icnt_p0+iseg+1)*8+2] += +n01x;
    aPxyNxyf[(icnt_p0+iseg+1)*8+3] += +n01y;
    aPxyNxyf[(icnt_p0+iseg+1)*8+6] += -n01x;
    aPxyNxyf[(icnt_p0+iseg+1)*8+7] += -n01y;
  }
  for(unsigned ipl=0;ipl<npl*2;++ipl){
    double nx0 = aPxyNxyf[icnt_p0*8+ipl*4+2];
    double ny0 = aPxyNxyf[icnt_p0*8+ipl*4+3];
    double len = sqrt(nx0*nx0+ny0*ny0);
    aPxyNxyf[icnt_p0*8+ipl*4+2] /= len;
    aPxyNxyf[icnt_p0*8+ipl*4+3] /= len;
  }
  const unsigned int iatri = aaTri.size();
  aaTri.resize(aaTri.size()+1);
  aaTri[iatri].resize(nseg*2*3);
  for(unsigned int iseg=0;iseg<nseg;++iseg){
    const size_t i0 = icnt_p0+iseg+0;
    const size_t i1 = icnt_p0+iseg+1;
    aaTri[iatri][iseg*6+0] = i0*2+0;
    aaTri[iatri][iseg*6+1] = i0*2+1;
    aaTri[iatri][iseg*6+2] = i1*2+0;
    aaTri[iatri][iseg*6+3] = i0*2+1;
    aaTri[iatri][iseg*6+4] = i1*2+0;
    aaTri[iatri][iseg*6+5] = i1*2+1;
  }
}

void AddPoint
(std::vector<float>& aPxyNxyf,
 std::vector< std::vector<unsigned int> >& aaTri,
 double x0,
 double y0,
 unsigned int ndiv)
{
  const unsigned int np0 = aPxyNxyf.size()/4;
  const unsigned int npa = ndiv+1;
  aPxyNxyf.resize(np0*4+npa*4);
  const double rdiv = 2.0*3.14156/ndiv;
  for(unsigned int ipa=0;ipa<npa;++ipa){
    if( ipa == 0 ){
      aPxyNxyf[np0*4+ipa*4+0] = x0;
      aPxyNxyf[np0*4+ipa*4+1] = y0;
      aPxyNxyf[np0*4+ipa*4+2] = 0.0;
      aPxyNxyf[np0*4+ipa*4+3] = 0.0;
    }
    else{
      aPxyNxyf[np0*4+ipa*4+0] = x0;
      aPxyNxyf[np0*4+ipa*4+1] = y0;
      aPxyNxyf[np0*4+ipa*4+2] = cos((ipa-1)*rdiv);
      aPxyNxyf[np0*4+ipa*4+3] = sin((ipa-1)*rdiv);
    }
  }
  unsigned int itria0 = aaTri.size();
  aaTri.resize(aaTri.size()+1);
  for(size_t idiv=0;idiv<ndiv;++idiv){
    aaTri[itria0].push_back(np0);
    aaTri[itria0].push_back(np0+1+(idiv+0)%ndiv);
    aaTri[itria0].push_back(np0+1+(idiv+1)%ndiv);
  }
}

void dfm2::opengl::CShader_Cad2D::MakeBuffer(const CCad2D& cad)
{
  {
    vao_face.aEBO.clear();
    vao_edge.aEBO.clear();
  }
  
  std::vector<float> aXY0f;
  std::vector<std::vector<unsigned int> > aaLine;
  std::vector<unsigned int> aTri;
  cad.MeshForVisualization(aXY0f,aaLine,aTri);
  { // triangles for faces
    if( !glIsVertexArray(vao_face.VAO) ){ glGenVertexArrays(1, &vao_face.VAO); }
    glBindVertexArray(vao_face.VAO);
    
    vao_face.ADD_VBO(0,aXY0f);
    glEnableVertexAttribArray(0);
    glVertexAttribPointer(0, 2, GL_FLOAT, GL_FALSE, 2*sizeof(float), (void*)0); // gl24
    /*
    {
      if( !glIsBuffer(vao_face.VBO_pos) ){ glGenBuffers(1, &vao_face.VBO_pos); }
      glBindBuffer(GL_ARRAY_BUFFER, vao_face.VBO_pos); // gl24
      glBufferData(GL_ARRAY_BUFFER, sizeof(float)*aXY0f.size(), aXY0f.data(), GL_STATIC_DRAW); // gl24
      glEnableVertexAttribArray(0);
      glVertexAttribPointer(0, 2, GL_FLOAT, GL_FALSE, 2*sizeof(float), (void*)0); // gl24
    }
     */

    vao_face.Delete_EBOs();
    {
      unsigned int EBO_Tri;
      glGenBuffers(1, &EBO_Tri);
      glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, EBO_Tri);
      glBufferData(GL_ELEMENT_ARRAY_BUFFER, sizeof(unsigned int)*aTri.size(), aTri.data(), GL_STATIC_DRAW);
      CGL4_VAO_Mesh::CEBO e0;
      e0.size = aTri.size();
      e0.GL_MODE = GL_TRIANGLES;
      e0.EBO = EBO_Tri;
      vao_face.aEBO.push_back(e0);
    }
  }
  { // point and edges
    std::vector<float> aPxyNxyf;
    std::vector< std::vector<unsigned int>> aaTri;
    for(const auto & iv : cad.aVtx){
      AddPoint(aPxyNxyf, aaTri,
               iv.pos.x,
               iv.pos.y,
               16);
    }
    for(const auto & il : aaLine){
      AddLine(aPxyNxyf,aaTri,
              aXY0f,il);
    }
    if( !glIsVertexArray(vao_edge.VAO) ){ glGenVertexArrays(1, &vao_edge.VAO); }
    glBindVertexArray(vao_edge.VAO);
    
    vao_edge.ADD_VBO(0,aPxyNxyf);
    glVertexAttribPointer(0, 2, GL_FLOAT, GL_FALSE, 4*sizeof(float), (void*)0); // gl24
    glEnableVertexAttribArray(0);
    glVertexAttribPointer(1, 2, GL_FLOAT, GL_FALSE, 4*sizeof(float), (void*)(2*sizeof(float))); // gl24
    glEnableVertexAttribArray(1);
    /*
    {
      if( !glIsBuffer(vao_edge.VBO_pos) ){ glGenBuffers(1, &vao_edge.VBO_pos); }
      glBindBuffer(GL_ARRAY_BUFFER, vao_edge.VBO_pos); // gl24
      glBufferData(GL_ARRAY_BUFFER, sizeof(float)*aPxyNxyf.size(), aPxyNxyf.data(), GL_STATIC_DRAW); // gl24
    }
     */
    ///
    vao_edge.Delete_EBOs();
    for(auto & il : aaTri){
      unsigned int EBO_Tri;
      glGenBuffers(1, &EBO_Tri);
      glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, EBO_Tri);
      glBufferData(GL_ELEMENT_ARRAY_BUFFER, sizeof(unsigned int)*il.size(), il.data(), GL_STATIC_DRAW);
      CGL4_VAO_Mesh::CEBO e0;
      e0.size = il.size();
      e0.GL_MODE = GL_TRIANGLES;
      e0.EBO = EBO_Tri;
      vao_edge.aEBO.push_back(e0);
    }
  }
}


void dfm2::opengl::CShader_Cad2D::Compile_Face()
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


void dfm2::opengl::CShader_Cad2D::Compile_Edge()
{

  const std::string glsl33vert_projection =
  "uniform mat4 matrixProjection;\n"
  "uniform mat4 matrixModelView;\n"
  "uniform float line_width;\n"
  "layout (location = 0) in vec2 posIn;\n"
  "layout (location = 1) in vec2 nrmIn;\n"
  "void main()\n"
  "{\n"
  "  vec2 v0 = posIn + nrmIn*0.5*line_width;\n"
  "  gl_Position = matrixProjection * matrixModelView * vec4(v0.x, v0.y, line_width*0.01, 1.0);\n"
  "}\0";

  const std::string glsl33frag =
  "uniform vec3 color;\n"
  "out vec4 FragColor;\n"
  "void main()\n"
  "{\n"
  "  FragColor = vec4(color.x, color.y, color.z, 1.0f);\n"
  "}\n\0";
  
  #ifdef EMSCRIPTEN
    shdr1_program = GL24_CompileShader((std::string("#version 300 es\n")+
                                        glsl33vert_projection).c_str(),
                                       (std::string("#version 300 es\n")+
                                        std::string("precision highp float;\n")+
                                        glsl33frag).c_str());
  #else
    shdr1_program = GL24_CompileShader((std::string("#version 330 core\n")+
                                        glsl33vert_projection).c_str(),
                                       (std::string("#version 330 core\n")+
                                        glsl33frag).c_str());
  #endif
    
    
  assert( glIsProgram(shdr1_program) );
  glUseProgram(shdr1_program);
  shdr1_Loc_MatrixProjection = glGetUniformLocation(shdr1_program,  "matrixProjection");
  shdr1_Loc_MatrixModelView  = glGetUniformLocation(shdr1_program,  "matrixModelView");
  shdr1_Loc_Color            = glGetUniformLocation(shdr1_program,  "color");
  shdr1_Loc_LineWidth        = glGetUniformLocation(shdr1_program,  "line_width");
  std::cout << "  shaderProgram: " << shdr1_program;
  std::cout << "  projectionMatrixLoc: " << shdr1_Loc_MatrixProjection;
  std::cout << "  LocColor: " << shdr1_Loc_Color << std::endl;
}

void dfm2::opengl::CShader_Cad2D::Draw(
    const float mP[16],
    const float mMV[16],
    const CCad2D& cad) const
{
//  assert( vao_face.aElem.size() >= 2 );
  if( is_show_face ){
    glUseProgram(shdr0_program);
    glUniformMatrix4fv(shdr0_Loc_MatrixProjection, 1, GL_FALSE, mP);
    glUniformMatrix4fv(shdr0_Loc_MatrixModelView, 1, GL_FALSE, mMV);
    glUniform3f(shdr0_Loc_Color, 1,1,1);
    vao_face.Draw(0);
  }
  
  glUseProgram(shdr1_program);
  glUniformMatrix4fv(shdr1_Loc_MatrixProjection, 1, GL_FALSE, mP);
  glUniformMatrix4fv(shdr1_Loc_MatrixModelView, 1, GL_FALSE, mMV);
  
  assert( vao_edge.aEBO.size() == cad.aVtx.size()+cad.aEdge.size() );
  
  const int ipicked_iv = cad.ivtx_picked;
  const int ipicked_ie = cad.iedge_picked;
//  const int ipicked_elem = cad.ipicked_elem;
    
  const double view_height = 1.0/mP[5];
  const int nv = cad.aVtx.size();
  /////
  glUniform1f(shdr1_Loc_LineWidth, view_height*0.04);
  for(size_t iv=0;iv<cad.aVtx.size();++iv){
    if( (int)iv == ipicked_iv ){ glUniform3f(shdr1_Loc_Color, 1.0f,0.9f,0.f); }
    else{                   glUniform3f(shdr1_Loc_Color, 1.f,0.f,0.f); }
    vao_edge.Draw(iv);
  }
  ///
  glUniform1f(shdr1_Loc_LineWidth, view_height*0.02);
  for(size_t ie=0;ie<cad.aEdge.size();++ie){
    if( (int)ie == ipicked_ie ){ glUniform3f(shdr1_Loc_Color, 1.0f,0.9f,0.f); }
    else{                   glUniform3f(shdr1_Loc_Color, 0.f,0.f,0.f); }
    vao_edge.Draw(nv+ie);
  }
}






// ------------------------------------------------------------------------------


void dfm2::opengl::CShader_MeshDTri2D::MakeBuffer
 (const std::vector<CVec2d>& aVec2,
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
    CGL4_VAO_Mesh::CEBO e0;
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
    CGL4_VAO_Mesh::CEBO e0;
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
  "}\0";

  const std::string glsl33frag =
  "uniform vec3 color;\n"
  "out vec4 FragColor;\n"
  "void main()\n"
  "{\n"
  "  FragColor = vec4(color.x, color.y, color.z, 1.0f);\n"
  "}\n\0";
  
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
