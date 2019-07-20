//
//  glew_funcs.cpp
//  c_gl
//
//  Created by Nobuyuki Umetani on 2019-07-19.
//

#include <stdio.h>

#include "delfem2/glew_funcs.h"


void CElemBuffObj::SetBuffer_Elem(const std::vector<unsigned int>& aTri, unsigned int gl_elem_type)
{
  this->gl_elem_type = gl_elem_type;
  size_elem = aTri.size();
  glGenBuffers(1,&iebo);
  glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, iebo);
  glBufferData(GL_ELEMENT_ARRAY_BUFFER,
               aTri.size()*(sizeof(int)),
               aTri.data(),
               GL_STATIC_DRAW);
  glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0);
}

void CElemBuffObj::DrawBuffer() const
{
  glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, this->iebo);
  glDrawElements(gl_elem_type, size_elem, GL_UNSIGNED_INT, 0);
}


void CGLBuffer::SetBuffer_Vtx(const std::vector<double>& aXYZ, int ndim)
{
  this->ndim = ndim;
  glGenBuffers(1,&vbo);
  glBindBuffer(GL_ARRAY_BUFFER, vbo);
  glBufferData(GL_ARRAY_BUFFER,
               aXYZ.size() * (sizeof(double)),
               aXYZ.data(),
               GL_STATIC_DRAW);
  glBindBuffer(GL_ARRAY_BUFFER, 0);
}

void CGLBuffer::SetBuffer_Nrm(const std::vector<double>& aNrm)
{
  glGenBuffers(1,&vbo_nrm);
  glBindBuffer(GL_ARRAY_BUFFER, vbo_nrm);
  glBufferData(GL_ARRAY_BUFFER,
               aNrm.size() * (sizeof(double)),
               aNrm.data(),
               GL_STATIC_DRAW);
  glBindBuffer(GL_ARRAY_BUFFER, 0);
}


void CGLBuffer::Draw_Start() const
{
  assert( glIsBuffer(this->vbo) );
  glEnableClientState(GL_VERTEX_ARRAY);
  glBindBuffer(GL_ARRAY_BUFFER, this->vbo);
  glVertexPointer(ndim, GL_DOUBLE, 0, 0);
  if( glIsBuffer(this->vbo_nrm) ){
    glEnableClientState(GL_NORMAL_ARRAY);
    glBindBuffer(GL_ARRAY_BUFFER, this->vbo_nrm);
    glNormalPointer(GL_DOUBLE, 0, 0);
  }
}
void CGLBuffer::Draw_End() const
{
  glDisableClientState(GL_NORMAL_ARRAY);
  glDisableClientState(GL_VERTEX_ARRAY);
}
