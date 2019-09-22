/*
 * Copyright (c) 2019 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */

#include <stdio.h>

/*
#if defined(__APPLE__) && defined(__MACH__)
#include <GL/glew.h>
//#include <OpenGL/glext.h>
#include <OpenGL/gl.h>
#else
#include <GL/glew.h>
#endif
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

#include "delfem2/gl4_funcs.h"



/////////////////////////////////////////////


/**
 * @details the OpenGL ES 2.0 only accept float array. So there is no "double" version of this file
 */
int GL4_VAO_MeshTri3D
(const float* aP, int nP, int nDim,
 const unsigned int* aTri, int nTri)
{
  unsigned int VBO, VAO, EBO;
  // set up vertex data (and buffer(s)) and configure vertex attributes
  glGenVertexArrays(1, &VAO); // opengl4
  glGenBuffers(1, &VBO);
  glGenBuffers(1, &EBO);
  // bind the Vertex Array Object first, then bind and set vertex buffer(s), and then configure vertex attributes(s).
  glBindVertexArray(VAO); // opengl4
  
  glBindBuffer(GL_ARRAY_BUFFER, VBO); // gl24
  glBufferData(GL_ARRAY_BUFFER, sizeof(float)*nP*nDim, aP, GL_STATIC_DRAW); // gl24
  glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, EBO); // gl24
  glBufferData(GL_ELEMENT_ARRAY_BUFFER, sizeof(unsigned int)*nTri*3, aTri, GL_STATIC_DRAW); // gl24
  glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 3*sizeof(float), (void*)0); // gl24
  glEnableVertexAttribArray(0); // gl24
  
  // note that this is allowed, the call to glVertexAttribPointer registered VBO as the vertex attribute's bound vertex buffer object so afterwards we can safely unbind
  glBindBuffer(GL_ARRAY_BUFFER, 0); // gl24
  
  // remember: do NOT unbind the EBO while a VAO is active as the bound element buffer object IS stored in the VAO; keep the EBO bound.
  //glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0);
  
  // You can unbind the VAO afterwards so other VAO calls won't accidentally modify this VAO, but this rarely happens. Modifying other
  // VAOs requires a call to glBindVertexArray anyways so we generally don't unbind VAOs (nor VBOs) when it's not directly necessary.
  glBindVertexArray(0); // gl4
  return VAO;
}
