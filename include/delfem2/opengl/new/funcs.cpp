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
#include "delfem2/opengl/new/funcs.h"

// ---------------------------------------------------------

DFM2_INLINE void delfem2::opengl::GL4_VAO_Pos(
    unsigned int& VAO,
    unsigned int& VBO,
    const float* aP,
    int nP,
    int nDim)
{
  glGenVertexArrays(1, &VAO); // opengl4
  glBindVertexArray(VAO); // opengl4
  
  glGenBuffers(1, &VBO);
  glBindBuffer(GL_ARRAY_BUFFER, VBO); // gl24
  glBufferData(
      GL_ARRAY_BUFFER,
      static_cast<GLsizeiptr>(sizeof(float)*nP*nDim),
      aP, GL_STATIC_DRAW); // gl24
  glEnableVertexAttribArray(0);
  glVertexAttribPointer(
      0, nDim, GL_FLOAT, GL_FALSE,
      static_cast<GLsizei>(nDim*sizeof(float)),
      (void*)nullptr); // gl24
  
  glBindBuffer(GL_ARRAY_BUFFER, 0); // gl24
  
  glBindVertexArray(0); // gl4
}


DFM2_INLINE void delfem2::opengl::GL4_VAO_PosNrm(
    unsigned int& VAO,
    unsigned int& VBO_pos,
    unsigned int& VBO_nrm,
    const float* aP,
    int nP,
    int nDim,
    const float* aN)
{
  // set up vertex data (and buffer(s)) and configure vertex attributes
  glGenVertexArrays(1, &VAO); // opengl4
  // bind the Vertex Array Object first,
  // then bind and set vertex buffer(s),
  // and then configure vertex attributes(s).
  glBindVertexArray(VAO); // opengl4
  
  // vertex position array
  glGenBuffers(1, &VBO_pos);
  glBindBuffer(GL_ARRAY_BUFFER, VBO_pos); // gl24
  glBufferData(
      GL_ARRAY_BUFFER,
      static_cast<GLsizeiptr>(sizeof(float)*nP*nDim),
      aP, GL_STATIC_DRAW); // gl24
  glEnableVertexAttribArray(0);
  glVertexAttribPointer(
      0, 3, GL_FLOAT, GL_FALSE,
      3*sizeof(float), (void*)nullptr); // gl24
  
  glGenBuffers(1, &VBO_nrm);
  glBindBuffer(GL_ARRAY_BUFFER, VBO_nrm); // gl24
  glBufferData(
      GL_ARRAY_BUFFER,
      static_cast<GLsizeiptr>(sizeof(float)*nP*nDim),
      aN, GL_STATIC_DRAW); // gl24
  glEnableVertexAttribArray(1);
  glVertexAttribPointer(
      1, 3, GL_FLOAT, GL_FALSE,
      3*sizeof(float), (void*)nullptr); // gl24
  
  glBindBuffer(GL_ARRAY_BUFFER, 0); // gl24
  
//  glEnableVertexAttribArray(0); // gl24
  
  // note that this is allowed, the call to glVertexAttribPointer registered VBO as the vertex attribute's bound vertex buffer object so afterwards we can safely unbind
  
  // remember: do NOT unbind the EBO while a VAO is active as the bound element buffer object IS stored in the VAO; keep the EBO bound.
  //glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0);
  
  // You can unbind the VAO afterwards so other VAO calls won't accidentally modify this VAO, but this rarely happens. Modifying other
  // VAOs requires a call to glBindVertexArray anyways so we generally don't unbind VAOs (nor VBOs) when it's not directly necessary.
  glBindVertexArray(0);
}


DFM2_INLINE void delfem2::opengl::VertexArrayObject::Draw(
    unsigned int iel) const
{
  if( iel >= aEBO.size() ){ return; }
  glBindVertexArray(VAO);
  glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, aEBO[iel].EBO);
  glDrawElements(aEBO[iel].GL_MODE,
                 static_cast<GLsizei>(aEBO[iel].size),
                 GL_UNSIGNED_INT,
                 nullptr);
}


// ------------------------------------

template <typename REAL>
void delfem2::opengl::VertexArrayObject::Add_VBO(
    unsigned int idx_vbo,
    const REAL *vtx_coords,
    size_t num_dof){
  glBindVertexArray(this->VAO); // opengl4
  assert( glIsVertexArray(this->VAO) );
  if( idx_vbo >= aVBO.size() ){ aVBO.resize(idx_vbo+1); }
  if( !glIsBuffer(aVBO[idx_vbo].VBO) ){ glGenBuffers(1, &aVBO[idx_vbo].VBO); }
  glBindBuffer(GL_ARRAY_BUFFER, aVBO[idx_vbo].VBO); // gl24
  glBufferData(
      GL_ARRAY_BUFFER,
      sizeof(REAL)*num_dof,
      vtx_coords,
      GL_STATIC_DRAW); // gl24
}
#ifdef DFM2_STATIC_LIBRARY
template void delfem2::opengl::VertexArrayObject::Add_VBO(
    unsigned int idx_vbo,
    const float *vtx_coords,
    size_t num_dof);
template void delfem2::opengl::VertexArrayObject::Add_VBO(
    unsigned int idx_vbo,
    const double *vtx_coords,
    size_t num_dof);
#endif

// ------------------------------------

DFM2_INLINE void delfem2::opengl::VertexArrayObject::Add_EBO(
    const std::vector<unsigned int>& elem_vtx,
    int gl_primitive_type)
{
  glBindVertexArray(this->VAO); // opengl4
  assert( glIsVertexArray(this->VAO) );
  assert( gl_primitive_type != GL_QUADS ); // quad is for legacy OpenGL
  // -----
  unsigned int idEBO;
  glGenBuffers(1, &idEBO);
  glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, idEBO); // gl24
  glBufferData(
      GL_ELEMENT_ARRAY_BUFFER,
      static_cast<GLsizeiptr>(sizeof(unsigned int)*elem_vtx.size()),
      elem_vtx.data(),
      GL_STATIC_DRAW); // gl24
  VertexArrayObject::CEBO e0{
    gl_primitive_type,
    elem_vtx.size(),
    idEBO};
  this->aEBO.push_back(e0);
}


void delfem2::opengl::VertexArrayObject::Delete_EBOs(){
  for(auto & ie : aEBO){
    unsigned int ebo = ie.EBO;
    if( glIsBuffer(ebo) ){ glDeleteBuffers(1,&ebo); }
  }
  aEBO.clear();
}
