/**
* Copyright (c) 2019 Nobuyuki Umetani
*
* This source code is licensed under the MIT license found in the
* LICENSE file in the root directory of this source tree.
*/

#include <stdio.h>

#include <glad/glad.h>
#include <GLFW/glfw3.h>

#ifdef EMSCRIPTEN
  #include <emscripten/emscripten.h>
  #define GLFW_INCLUDE_ES3
#endif

#include "delfem2/cad2d.h"
#include "delfem2/gl24_funcs.h"
#include "delfem2/gl4_funcs.h"
#include "delfem2/gl4_v23dtricad.h"


void CShader_CCad2D::MakeBuffer(const CCad2D& cad)
{
  std::vector<float> aXY0f;
  std::vector<std::vector<unsigned int> > aaLine;
  std::vector<unsigned int> aTri;
  cad.MeshForVisualization(aXY0f,aaLine,aTri);
  {
    unsigned int VAO;
    glGenVertexArrays(1, &VAO); // opengl4
    glBindVertexArray(VAO); // opengl4

    unsigned int VBO_pos;
    glGenBuffers(1, &VBO_pos);
    
    glBindBuffer(GL_ARRAY_BUFFER, VBO_pos); // gl24
    glBufferData(GL_ARRAY_BUFFER, sizeof(float)*aXY0f.size(), aXY0f.data(), GL_STATIC_DRAW); // gl24
    glEnableVertexAttribArray(0);
    glVertexAttribPointer(0, 2, GL_FLOAT, GL_FALSE, 2*sizeof(float), (void*)0); // gl24
    ///
    unsigned int EBO_Tri;
    glGenBuffers(1, &EBO_Tri);
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, EBO_Tri);
    glBufferData(GL_ELEMENT_ARRAY_BUFFER, sizeof(unsigned int)*aTri.size(), aTri.data(), GL_STATIC_DRAW);
    
    vao_face.VAO = VAO;
    vao_face.VBO_pos = VBO_pos;
    
    CGL4_VAO_Mesh::CElem e0;
    e0.size = aTri.size();
    e0.GL_MODE = GL_TRIANGLES;
    e0.EBO = EBO_Tri;
    vao_face.aElem.push_back(e0);
  }
  {
    std::vector<float> aPxyNxyf;
    {
      int icnt = 0;
      for(int il=0;il<aaLine.size();++il){ icnt += aaLine[il].size(); }
      aPxyNxyf.assign(icnt*8,0.0);
    }
    {
      int icnt = 0;
      for(int il=0;il<aaLine.size();++il){
        const int npl = aaLine[il].size();
        for(int ipl=0;ipl<npl;++ipl){
          aPxyNxyf[(icnt+ipl)*8+0] = aXY0f[aaLine[il][ipl]*2+0];
          aPxyNxyf[(icnt+ipl)*8+1] = aXY0f[aaLine[il][ipl]*2+1];
          aPxyNxyf[(icnt+ipl)*8+2] = 0.0;
          aPxyNxyf[(icnt+ipl)*8+3] = 0.0;
          aPxyNxyf[(icnt+ipl)*8+4] = aXY0f[aaLine[il][ipl]*2+0];
          aPxyNxyf[(icnt+ipl)*8+5] = aXY0f[aaLine[il][ipl]*2+1];
          aPxyNxyf[(icnt+ipl)*8+6] = 0.0;
          aPxyNxyf[(icnt+ipl)*8+7] = 0.0;
        }
        const int nseg = aaLine[il].size()-1;
        for(int iseg=0;iseg<nseg;++iseg){
          const int i0 = aaLine[il][iseg+0];
          const int i1 = aaLine[il][iseg+1];
          double v01x = aXY0f[i1*2+0] - aXY0f[i0*2+0];
          double v01y = aXY0f[i1*2+1] - aXY0f[i0*2+1];
          double len = sqrt(v01x*v01x+v01y*v01y);
          double n01x = +v01y/len;
          double n01y = -v01x/len;
          aPxyNxyf[(icnt+iseg+0)*8+2] += +n01x;
          aPxyNxyf[(icnt+iseg+0)*8+3] += +n01y;
          aPxyNxyf[(icnt+iseg+0)*8+6] += -n01x;
          aPxyNxyf[(icnt+iseg+0)*8+7] += -n01y;
          aPxyNxyf[(icnt+iseg+1)*8+2] += +n01x;
          aPxyNxyf[(icnt+iseg+1)*8+3] += +n01y;
          aPxyNxyf[(icnt+iseg+1)*8+6] += -n01x;
          aPxyNxyf[(icnt+iseg+1)*8+7] += -n01y;
        }
        icnt += aaLine[il].size();
      }
    }
    for(int i=0;i<aPxyNxyf.size()/4;++i){
      double nx0 = aPxyNxyf[i*4+2];
      double ny0 = aPxyNxyf[i*4+3];
      double len = sqrt(nx0*nx0+ny0*ny0);
      aPxyNxyf[i*4+2] /= len;
      aPxyNxyf[i*4+3] /= len;
    }
    std::vector< std::vector<unsigned int>> aaTri(aaLine.size());
    {
      int icnt = 0;
      for(int il=0;il<aaLine.size();++il){
        const int nseg = aaLine[il].size()-1;
        aaTri[il].resize(nseg*2*3);
        for(int iseg=0;iseg<nseg;++iseg){
          int i0 = icnt+iseg+0;
          int i1 = icnt+iseg+1;
          aaTri[il][iseg*6+0] = i0*2+0;
          aaTri[il][iseg*6+1] = i0*2+1;
          aaTri[il][iseg*6+2] = i1*2+0;
          aaTri[il][iseg*6+3] = i0*2+1;
          aaTri[il][iseg*6+4] = i1*2+0;
          aaTri[il][iseg*6+5] = i1*2+1;
        }
        icnt += aaLine[il].size();
      }
    }
    /*
    for(int il=0;il<aaTri.size();++il){
      for(int it=0;it<aaTri[il].size()/3;++it){
        std::cout << il << " " << it << " " << aaTri[il][it*3+0] << " " << aaTri[il][it*3+1] << " " << aaTri[il][it*3+2] << std::endl;
      }
    }
     */
    glGenVertexArrays(1, &vao_edge.VAO); // opengl4
    glBindVertexArray(vao_edge.VAO); // opengl4

    glGenBuffers(1, &vao_edge.VBO_pos);
    
    glBindBuffer(GL_ARRAY_BUFFER, vao_edge.VBO_pos); // gl24
    glBufferData(GL_ARRAY_BUFFER, sizeof(float)*aPxyNxyf.size(), aPxyNxyf.data(), GL_STATIC_DRAW); // gl24
    glVertexAttribPointer(0, 2, GL_FLOAT, GL_FALSE, 4*sizeof(float), (void*)0); // gl24
    glEnableVertexAttribArray(0);
    glVertexAttribPointer(1, 2, GL_FLOAT, GL_FALSE, 4*sizeof(float), (void*)(2*sizeof(float))); // gl24
    glEnableVertexAttribArray(1);
    ///
    for(int il=0;il<aaTri.size();++il){
      unsigned int EBO_Tri;
      glGenBuffers(1, &EBO_Tri);
      glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, EBO_Tri);
      glBufferData(GL_ELEMENT_ARRAY_BUFFER, sizeof(unsigned int)*aaTri[il].size(), aaTri[il].data(), GL_STATIC_DRAW);
      CGL4_VAO_Mesh::CElem e0;
      e0.size = aaTri[il].size();
      e0.GL_MODE = GL_TRIANGLES;
      e0.EBO = EBO_Tri;
      vao_edge.aElem.push_back(e0);
    }
    
  }
}


void CShader_CCad2D::Compile_Face()
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
  std::cout << "projectionMatrixLoc: " << shdr0_Loc_MatrixProjection << "   shaderProgram: " << shdr0_program << "  LocColor: " << shdr0_Loc_Color << std::endl;
}


void CShader_CCad2D::Compile_Edge()
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
  "  gl_Position = matrixProjection * matrixModelView * vec4(v0.x, v0.y, 0.0, 1.0);\n"
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
  std::cout << "projectionMatrixLoc: " << shdr1_Loc_MatrixProjection << "   shaderProgram: " << shdr1_program << "  LocColor: " << shdr1_Loc_Color << std::endl;
}


void CShader_CCad2D::Draw(float mP[16], float mMV[16]) const
{
//  assert( vao_face.aElem.size() >= 2 );
  glUseProgram(shdr0_program);
  glUniformMatrix4fv(shdr0_Loc_MatrixProjection, 1, GL_FALSE, mP);
  glUniformMatrix4fv(shdr0_Loc_MatrixModelView, 1, GL_FALSE, mMV);
  glUniform3f(shdr0_Loc_Color, 1,1,1);
  vao_face.Draw(0);
  
  glUseProgram(shdr1_program);
  glUniformMatrix4fv(shdr1_Loc_MatrixProjection, 1, GL_FALSE, mP);
  glUniformMatrix4fv(shdr1_Loc_MatrixModelView, 1, GL_FALSE, mMV);
  glUniform3f(shdr1_Loc_Color, 0,0,0);
  double view_height = 1.0/mP[5];
  glUniform1f(shdr1_Loc_LineWidth, view_height*0.02);
  for(int ie=0;ie<vao_edge.aElem.size();++ie){
    vao_edge.Draw(ie);
  }
}
