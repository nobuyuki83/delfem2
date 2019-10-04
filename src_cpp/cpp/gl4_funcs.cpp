/*
 * Copyright (c) 2019 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */

#include <stdio.h>

#ifdef USE_GLEW
  #include <GL/glew.h>
#else
  #include <glad/glad.h>
  #ifdef EMSCRIPTEN
    #include <emscripten/emscripten.h>
    #define GLFW_INCLUDE_ES3
  #endif
#endif

#include "delfem2/gl4_funcs.h"

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
  return;
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
  glBindVertexArray(0); // gl4
  return;
}


void CGL4_VAO_Mesh::Draw(unsigned int iel) const {
  if( iel >= aElem.size() ){ assert(0); return; }
  glBindVertexArray(VAO);
  glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, aElem[iel].EBO);
  glDrawElements(aElem[iel].GL_MODE,
                 aElem[iel].size,
                 GL_UNSIGNED_INT,
                 0);
}

