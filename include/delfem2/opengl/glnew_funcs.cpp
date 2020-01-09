/*
 * Copyright (c) 2019 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */

#include <cstdio>

#ifdef USE_GLEW
  #include <GL/glew.h>
#else
  #include <glad/glad.h>
  #ifdef EMSCRIPTEN
    #include <emscripten/emscripten.h>
    #define GLFW_INCLUDE_ES3
  #endif
#endif

#include "delfem2/opengl/glnew_funcs.h"

// ---------------------------------------------------------

void GL4_VAO_Pos
(unsigned int& VAO,
 unsigned int& VBO,
 const float* aP, int nP, int nDim)
{
  glGenVertexArrays(1, &VAO); // opengl4
  glBindVertexArray(VAO); // opengl4
  
  glGenBuffers(1, &VBO);
  glBindBuffer(GL_ARRAY_BUFFER, VBO); // gl24
  glBufferData(GL_ARRAY_BUFFER, sizeof(float)*nP*nDim, aP, GL_STATIC_DRAW); // gl24
  glEnableVertexAttribArray(0);
  glVertexAttribPointer(0, nDim, GL_FLOAT, GL_FALSE, nDim*sizeof(float), (void*)0); // gl24
  
  glBindBuffer(GL_ARRAY_BUFFER, 0); // gl24
  
  glBindVertexArray(0); // gl4
}


void GL4_VAO_PosNrm
(unsigned int& VAO,
 unsigned int& VBO_pos,
 unsigned int& VBO_nrm,
 const float* aP, int nP, int nDim,
 const float* aN)
{
  // set up vertex data (and buffer(s)) and configure vertex attributes
  glGenVertexArrays(1, &VAO); // opengl4
  // bind the Vertex Array Object first, then bind and set vertex buffer(s), and then configure vertex attributes(s).
  glBindVertexArray(VAO); // opengl4
  
  // vertex position array
  glGenBuffers(1, &VBO_pos);
  glBindBuffer(GL_ARRAY_BUFFER, VBO_pos); // gl24
  glBufferData(GL_ARRAY_BUFFER, sizeof(float)*nP*nDim, aP, GL_STATIC_DRAW); // gl24
  glEnableVertexAttribArray(0);
  glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 3*sizeof(float), (void*)0); // gl24
  
  glGenBuffers(1, &VBO_nrm);
  glBindBuffer(GL_ARRAY_BUFFER, VBO_nrm); // gl24
  glBufferData(GL_ARRAY_BUFFER, sizeof(float)*nP*nDim, aN, GL_STATIC_DRAW); // gl24
  glEnableVertexAttribArray(1);
  glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, 3*sizeof(float), (void*)0); // gl24
  
  glBindBuffer(GL_ARRAY_BUFFER, 0); // gl24
  
//  glEnableVertexAttribArray(0); // gl24
  
  // note that this is allowed, the call to glVertexAttribPointer registered VBO as the vertex attribute's bound vertex buffer object so afterwards we can safely unbind
  
  // remember: do NOT unbind the EBO while a VAO is active as the bound element buffer object IS stored in the VAO; keep the EBO bound.
  //glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0);
  
  // You can unbind the VAO afterwards so other VAO calls won't accidentally modify this VAO, but this rarely happens. Modifying other
  // VAOs requires a call to glBindVertexArray anyways so we generally don't unbind VAOs (nor VBOs) when it's not directly necessary.
  glBindVertexArray(0);
}


void CGL4_VAO_Mesh::Draw(unsigned int iel) const {
  if( iel >= aEBO.size() ){ assert(0); return; }
  glBindVertexArray(VAO);
  glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, aEBO[iel].EBO);
  glDrawElements(aEBO[iel].GL_MODE,
                 aEBO[iel].size,
                 GL_UNSIGNED_INT,
                 0);
}


void CGL4_VAO_Mesh::ADD_VBO
(unsigned int ivbo,
 const std::vector<float>& aXY0f)
{
  glBindVertexArray(this->VAO); // opengl4
  assert( glIsVertexArray(this->VAO) );
  // ----
  if( ivbo >= aVBO.size() ){ aVBO.resize(ivbo+1); }
  if( !glIsBuffer(aVBO[ivbo].VBO) ){ glGenBuffers(1, &aVBO[ivbo].VBO); }
  glBindBuffer(GL_ARRAY_BUFFER, aVBO[ivbo].VBO); // gl24
  glBufferData(GL_ARRAY_BUFFER, sizeof(float)*aXY0f.size(), aXY0f.data(), GL_STATIC_DRAW); // gl24
}


void CGL4_VAO_Mesh::Add_EBO
 (const std::vector<unsigned int>& aElem,
  int GL_MODE)
{
  glBindVertexArray(this->VAO); // opengl4
  assert( glIsVertexArray(this->VAO) );
  // -----
  unsigned int EBO_Tri;
  glGenBuffers(1, &EBO_Tri);
  glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, EBO_Tri); // gl24
  glBufferData(GL_ELEMENT_ARRAY_BUFFER, sizeof(unsigned int)*aElem.size(), aElem.data(), GL_STATIC_DRAW); // gl24
  CGL4_VAO_Mesh::CEBO e0;
  e0.EBO = EBO_Tri;
  e0.GL_MODE = GL_MODE;
  e0.size = aElem.size();
  this->aEBO.push_back(e0);
}
