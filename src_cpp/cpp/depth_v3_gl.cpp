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
#include "delfem2/funcs_gl.h"
#include "delfem2/v23_gl.h"  // vec3, mat3

#include "delfem2/depth_v3_gl.h"

/////////////////////////////////

void CGPUSampler::SetColor(double r, double g, double b){
  color[0] = r;
  color[1] = g;
  color[2] = b;
}

void CGPUSampler::Init(int nw, int nh, bool isColor, bool isDepth)
{
  this->nResX = nw;
  this->nResY = nh;
  this->isColor = isColor;
  this->isDepth = isDepth;
  const int npix = nw*nh;
  /////
  if( isDepth ){ aZ.resize(npix,0); }
  else{ aZ.clear(); }
  ////
  if( isColor ){ aRGBA.resize(npix*4,128); }
  else{ aRGBA.clear(); }
  ////////
  if( id_tex_color > 0 ){ glDeleteTextures(1, &id_tex_color); }
  id_tex_color = 0;
}

void CGPUSampler::SetCoord
(double elen, double depth_max,
 const std::vector<double>& org_prj,
 const std::vector<double>& dir_prj,
 const std::vector<double>& dir_width)
{
  this->lengrid = elen;
  this->z_range = depth_max;
  z_axis[0] = dir_prj[0];  z_axis[1] = dir_prj[1];  z_axis[2] = dir_prj[2];
  origin[0] = org_prj[0];  origin[1] = org_prj[1];  origin[2] = org_prj[2];
  x_axis[0] = dir_width[0];  x_axis[1] = dir_width[1];  x_axis[2] = dir_width[2];
}

void CGPUSampler::SetView(){
  ::glMatrixMode(GL_MODELVIEW);
  ::glLoadIdentity();
  ViewTransformation(x_axis,z_axis,origin);
  ::glMatrixMode(GL_PROJECTION);
  ::glLoadIdentity();
  ::glOrtho(0.0, +lengrid*nResX,
            0.0, +lengrid*nResY,
            0, z_range);
  ::glMatrixMode(GL_MODELVIEW);
}

void CGPUSampler::Start()
{
  glGetIntegerv(GL_VIEWPORT, view); // current viewport
  ::glViewport(0, 0, nResX, nResY);

  if(      bgcolor.size() == 4 ){ ::glClearColor(bgcolor[0], bgcolor[1], bgcolor[2], bgcolor[3]); }
  else if( bgcolor.size() == 3 ){ ::glClearColor(bgcolor[0], bgcolor[1], bgcolor[2], 1.0); }
  else if( bgcolor.size() > 0  ){ ::glClearColor(bgcolor[0], bgcolor[0], bgcolor[0], 1.0 ); }
  else{                           ::glClearColor(1.0, 1.0, 1.0, 1.0 ); }
  ::glClear(GL_COLOR_BUFFER_BIT|GL_DEPTH_BUFFER_BIT|GL_STENCIL_BUFFER_BIT);
  ::glDisable(GL_BLEND);
  ::glEnable(GL_DEPTH_TEST);
  this->SetView();
}

void CGPUSampler::End()
{
  const int npix = nResX*nResY;
  
  if( isDepth ){
    assert( (int)aZ.size() == npix );
    glReadBuffer(GL_DEPTH_ATTACHMENT);
    glReadPixels(0, 0, nResX, nResY, GL_DEPTH_COMPONENT, GL_FLOAT, aZ.data());
    for(int i=0;i<npix;++i){ aZ[i] *= (-1.0*z_range); }
  }
  else{ aZ.clear(); }
  
  ///////
  if( isColor ){
    assert( (int)aRGBA.size() == npix*4 );
    glReadBuffer(GL_COLOR_ATTACHMENT0);
    glReadPixels(0, 0, nResX, nResY, GL_RGBA, GL_UNSIGNED_BYTE, aRGBA.data());
  }
  else{ aRGBA.clear(); }
  
  ::glViewport(view[0], view[1], view[2], view[3]);
}

void CGPUSampler::LoadTex()
{
  if( id_tex_color == 0 ){
    glGenTextures(1, &id_tex_color);
  }
  glBindTexture(GL_TEXTURE_2D, id_tex_color);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);
  glTexEnvi(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_MODULATE);
  glTexEnvf(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_DECAL);
  glPixelStorei(GL_UNPACK_ALIGNMENT,1);
  if( (int)aRGBA.size() == nResX*nResY*4 ){
    glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA,
                 nResX, nResY, 0, GL_RGBA, GL_UNSIGNED_BYTE,
                 aRGBA.data());
  }
  glBindTexture(GL_TEXTURE_2D, 0);
}

void CGPUSampler::Draw() const {
  
  ::glPointSize(1);
  this->Draw_Point();
  /////
  ::glLineWidth(3);
  this->Draw_Axis();
  ////
  ::glLineWidth(1);
  ::glColor3d(0,0,0);
  this->Draw_BoundingBox();
  
  if( id_tex_color > 0 ){
    const CVector3& dx = x_axis;
    const CVector3& dy = Cross(z_axis,dx);
    const double lx = lengrid*nResX;
    const double ly = lengrid*nResY;
    CVector3 p0 = origin;
    CVector3 p1 = origin + lx*dx;
    CVector3 p2 = origin + lx*dx + ly*dy;
    CVector3 p3 = origin + ly*dy;
    ::glEnable(GL_TEXTURE_2D);
    ::glDisable(GL_LIGHTING);
    ::glBindTexture(GL_TEXTURE_2D, id_tex_color);
    ::glColor3d(1,1,1);
    ::glBegin(GL_QUADS);
    ::glTexCoord2d(0.0, 0.0); myGlVertex(p0);
    ::glTexCoord2d(1.0, 0.0); myGlVertex(p1);
    ::glTexCoord2d(1.0, 1.0); myGlVertex(p2);
    ::glTexCoord2d(0.0, 1.0); myGlVertex(p3);
    ::glEnd();
    ::glBindTexture(GL_TEXTURE_2D, 0);
    ::glDisable(GL_TEXTURE_2D);
  }
}

void CGPUSampler::Draw_Axis() const
{
  ::glMatrixMode(GL_MODELVIEW);
  ::glPushMatrix();
  ModelTransformation(x_axis, z_axis, origin);
  DrawAxis(draw_len_axis);
    ::glPopMatrix();
}

void CGPUSampler::Draw_BoundingBox() const
{
  ::glMatrixMode(GL_MODELVIEW);
  ::glPushMatrix();
  ModelTransformation(x_axis, z_axis, origin);
  ::glLineWidth(3);
  Draw_AABB3D_MinMaxXYZ_Edge(0.0, lengrid*nResX, 0.0, lengrid*nResY, 0.0, -z_range);
  ::glPopMatrix();
}

void CGPUSampler::Draw_Point() const
{
  ::glDisable(GL_LIGHTING);
  if( (int)aZ.size() != nResX*nResY ) return;
  if( color.size() == 3 ){ ::glColor3dv(color.data()); }
  if( color.size() == 4 ){ ::glColor4dv(color.data()); }
  ::glBegin(GL_POINTS);
  double o[3];
  for(int iy=0;iy<nResY;++iy){
    for(int ix=0;ix<nResX;++ix){
      this->getGPos(o,ix,iy,aZ[iy*nResX+ix]);
      ::glVertex3dv(o);
    }
  }
  ::glEnd();
}

void CGPUSampler::getGPos(double p[3], int ix, int iy, double depth) const
{
  const CVector3& dx = x_axis;
  const CVector3& dz = z_axis;
  const CVector3& dy = Cross(dz,dx);
  double lz = depth;
  double lx = (ix+0.5)*lengrid;
  double ly = (iy+0.5)*lengrid;
  CVector3 vp = lx*dx+ly*dy+lz*dz + origin;
  p[0] = vp.x;
  p[1] = vp.y;
  p[2] = vp.z;
}
