#include <iostream>
#include <fstream>
#include <sstream>
#include <cmath>
#include <stack>
#include <cstring>
#include <cstdlib>

#if defined(__APPLE__) && defined(__MACH__) // Mac
  #include <OpenGL/gl.h>
  #include <OpenGL/glu.h>
#elif defined(__MINGW32__) // probably I'm using Qt and don't want to use GLUT
  #include <GL/glu.h>
#elif defined(WIN32) // windows
  #include <windows.h>
  #include <GL/gl.h>
#else
  #include <GL/glew.h>
  #include <GL/gl.h>
#endif
////
#include "delfem2/mat3.h"
#include "delfem2/v23m3q.h"

#include "delfem2/v23_gl.h"  // vec3, mat3
#include "delfem2/funcs_glew.h"           // shader

#include "delfem2/depth_v3_gl.h"

//////////////////////////////////////

void CDepthContext::DeleteFrameBuffer()
{
  if( id_framebuffer != -1 ){
    glDeleteFramebuffers(1, &id_framebuffer);
    id_framebuffer = -1;
  }
  // TODO: delete depth_texture here
  if( id_depth_render_buffer != -1  ){
    glDeleteRenderbuffersEXT(1, &id_depth_render_buffer);
    id_depth_render_buffer = -1;
  }
}

void CDepthContext::SetFrameBufferSize
(int width, int height)
{
  DeleteFrameBuffer();
  glGenFramebuffers(1, &id_framebuffer);
  glBindFramebuffer(GL_FRAMEBUFFER, id_framebuffer);
  ////
  glDrawBuffer(GL_NONE);
  glReadBuffer(GL_NONE);
  ////
  glGenRenderbuffers(1, &id_depth_render_buffer);
  glBindRenderbuffer(GL_RENDERBUFFER, id_depth_render_buffer);
  glRenderbufferStorage(GL_RENDERBUFFER, GL_DEPTH_COMPONENT32F, width, height);
  glFramebufferRenderbuffer(GL_FRAMEBUFFER, GL_DEPTH_ATTACHMENT, GL_RENDERBUFFER, id_depth_render_buffer);
  ////
  GLenum status = glCheckFramebufferStatus(GL_FRAMEBUFFER) ;
  glBindFramebuffer(GL_FRAMEBUFFER, 0);
  glBindRenderbuffer(GL_RENDERBUFFER, 0);
  if(status != GL_FRAMEBUFFER_COMPLETE){
    std::cout << "error!: " << status << std::endl;
    std::cout << GL_FRAMEBUFFER_INCOMPLETE_ATTACHMENT << std::endl;
    std::cout << GL_FRAMEBUFFER_INCOMPLETE_MISSING_ATTACHMENT << std::endl;
    std::cout << GL_FRAMEBUFFER_UNSUPPORTED << std::endl;
    std::cout << GL_FRAMEBUFFER_INCOMPLETE_READ_BUFFER << std::endl;
    std::cout << GL_FRAMEBUFFER_ATTACHMENT_TEXTURE_CUBE_MAP_FACE << std::endl;
    std::cout << GL_FRAMEBUFFER_INCOMPLETE_ATTACHMENT << std::endl;
    std::cout << GL_FRAMEBUFFER_INCOMPLETE_READ_BUFFER << std::endl;
    std::cout << GL_FRAMEBUFFER_INCOMPLETE_DRAW_BUFFER << std::endl;
    return;
  }
}

void CDepth::SetColor(double r, double g, double b){
  color[0] = r;
  color[1] = g;
  color[2] = b;
}
/////////////////////////////////


void CDepth::SetView(){
  ::glMatrixMode(GL_MODELVIEW);
  ::glLoadIdentity();
  {
    const CVector3& dx = dirWidth;
    const CVector3& dz = dirPrj;
    const CVector3& dy = Cross(dz,dx);
    CMatrix3 R = Mat3(dx,dy,dz);
    CVector3 o = R.Trans()*orgPrj;
    double A[16];
    A[ 0] = dx.x;  A[ 1] = dy.x;  A[ 2] = dz.x;  A[ 3] = 0;
    A[ 4] = dx.y;  A[ 5] = dy.y;  A[ 6] = dz.y;  A[ 7] = 0;
    A[ 8] = dx.z;  A[ 9] = dy.z;  A[10] = dz.z;  A[11] = 0;
    A[12] = -o.x;  A[13] = -o.y;  A[14] = -o.z;  A[15] = 1;
    ::glMultMatrixd(A);
  }
  ::glMatrixMode(GL_PROJECTION);
  ::glLoadIdentity();
  ::glOrtho(-lengrid*nResW*0.5, +lengrid*nResW*0.5,
            -lengrid*nResH*0.5, +lengrid*nResH*0.5,
            0, depth_max);
  ::glMatrixMode(GL_MODELVIEW);
}


void CDepth::TakeDepthShot
(const CInputDepth& obj)
{
  GLint view[4]; glGetIntegerv(GL_VIEWPORT, view); // current viewport
  ::glViewport(0, 0, nResW, nResH);
  
  ::glClearColor(0.0f, 0.0f, 0.0f ,1.0f);
  ::glClear(GL_COLOR_BUFFER_BIT|GL_DEPTH_BUFFER_BIT|GL_STENCIL_BUFFER_BIT);
  ::glDisable(GL_BLEND);
  ::glEnable(GL_DEPTH_TEST);
  this->SetView();
  obj.Draw();
//  glFinish();
  aDepth.resize(nResW*nResH);
  
  glReadBuffer(GL_DEPTH_ATTACHMENT);
  glReadPixels(0, 0, nResW, nResH, GL_DEPTH_COMPONENT, GL_FLOAT, aDepth.data());
  
  ::glViewport(view[0], view[1], view[2], view[3]);
}

void CDepth::Draw_Point(bool is_draw_miss) const
{
  ::glDisable(GL_LIGHTING);
  if( aDepth.size() != nResW*nResH ) return;
  ::glColor3dv(color);
  ::glBegin(GL_POINTS);
  CVector3 o;
  for(int i=0;i<nResW*nResH;++i){
    this->getGPos(o,i,aDepth[i]);
    ::glVertex3d(o.x,o.y,o.z);
  }
  ::glEnd();
}

void CDepth::getGPos(CVector3& p, int i, double depth) const{
  const CVector3& dx = dirWidth;
  const CVector3& dz = dirPrj;
  const CVector3 dy = Cross(dz,dx);
  CMatrix3 R = Mat3(dx,dy,dz);
  int ih = i/nResW;
  int iw = i-ih*nResW;
  double lz = -depth*depth_max;
  double lx = (iw+0.5)*lengrid-0.5*nResW*lengrid;
  double ly = (ih+0.5)*lengrid-0.5*nResH*lengrid;
  p = R*CVector3(lx,ly,lz) + orgPrj;
}
