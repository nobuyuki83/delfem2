/*
 * Copyright (c) 2019 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */


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
#elif defined(_WIN32) // windows
  #include <windows.h>
  #include <GL/glew.h>
  #include <GL/gl.h>
#else
  #include <GL/glew.h>
  #include <GL/gl.h>
#endif
////
#include "delfem2/gl_funcs.h"
#include "delfem2/gl_v23q.h"  // vec3, mat3
#include "delfem2/gl_gpusampler.h"

/*
double Dot(const std::vector<double>& p0, const std::vector<double>& p1){
  const int n = p0.size();
  double v=0;
  for(int i=0;i<n;++i){ v += p0[i]*p1[i]; }
  return v;
}

void Scale(std::vector<double>& p0, double s)
{
  const int n = p0.size();
  for(int i=0;i<n;++i){ p0[i] *= s; }
}

void Normalize(std::vector<double>& p0)
{
  const double ss = Dot(p0,p0);
  Scale(p0,1.0/sqrt(ss));
}
 */

double Dot(const double* p0, const double* p1, int ndof){
  double v=0;
  for(int i=0;i<ndof;++i){ v += p0[i]*p1[i]; }
  return v;
}

void Scale(double* p0, int n, double s)
{
  for(int i=0;i<n;++i){ p0[i] *= s; }
}

void Normalize(double* p0, int n)
{
  const double ss = Dot(p0,p0,n);
  Scale(p0,n,1.0/sqrt(ss));
}



/////////////////////////////////

void CGPUSampler::SetColor(double r, double g, double b){
  color[0] = r;
  color[1] = g;
  color[2] = b;
}

void CGPUSampler::Init(int nw, int nh,
                       std::string sFormatPixelColor, bool isDepth)
{
  this->nResX = nw;
  this->nResY = nh;
  this->sFormatPixelColor = sFormatPixelColor;
  this->isDepth = isDepth;
  const int npix = nw*nh;
  /////
  if( isDepth ){ aZ.resize(npix,0); }
  else{ aZ.clear(); }
  ////
  aF_RGBA.clear();
  aUC_RGBA.clear();
  if( sFormatPixelColor == "4byte"  ){ aUC_RGBA.resize(npix*4,128); }
  else if( sFormatPixelColor == "4float" ){ aF_RGBA.resize(npix*4,128); }
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
  z_axis[0] = dir_prj[0];   z_axis[1] = dir_prj[1];   z_axis[2] = dir_prj[2];
  origin[0] = org_prj[0];   origin[1] = org_prj[1];   origin[2] = org_prj[2];
  x_axis[0] = dir_width[0]; x_axis[1] = dir_width[1]; x_axis[2] = dir_width[2];
  Normalize(z_axis,3);
  Normalize(x_axis,3);
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
  ::glClear(GL_COLOR_BUFFER_BIT|GL_DEPTH_BUFFER_BIT);
  ::glDisable(GL_BLEND);
  ::glEnable(GL_DEPTH_TEST);
  this->SetView();
}

void CGPUSampler::End()
{
  const int npix = nResX*nResY;
  if( isDepth ){
    assert( (int)aZ.size() == npix );
//    glReadBuffer(GL_DEPTH_ATTACHMENT); // this caused crash in PyOpenGL. I didn't understand but I comment this out.
    glReadPixels(0, 0, nResX, nResY, GL_DEPTH_COMPONENT, GL_FLOAT, aZ.data());
    for(int i=0;i<npix;++i){ aZ[i] *= (-1.0*z_range); }
  }
  else{
    aZ.clear();
  }
  ///////
  
  if( sFormatPixelColor == "4byte" || sFormatPixelColor == "4float" ){
    glReadBuffer(GL_COLOR_ATTACHMENT0);
    if( sFormatPixelColor == "4byte" ){
      assert( (int)aUC_RGBA.size() == npix*4 );
      glReadPixels(0, 0, nResX, nResY, GL_RGBA, GL_UNSIGNED_BYTE, aUC_RGBA.data());
      aF_RGBA.clear();
    }
    if( sFormatPixelColor == "4float" ){
      assert( (int)aF_RGBA.size() == npix*4 );
      glReadPixels(0, 0, nResX, nResY, GL_RGBA, GL_FLOAT, aF_RGBA.data());
      aUC_RGBA.clear();
    }
  }
  else{
    aUC_RGBA.clear();
    aF_RGBA.clear();
  }
  
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
  if( sFormatPixelColor == "4byte" ){
    assert( (int)aUC_RGBA.size() == nResX*nResY*4 );
    glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA,
                 nResX, nResY, 0, GL_RGBA, GL_UNSIGNED_BYTE,
                 aUC_RGBA.data());
  }
  else if( sFormatPixelColor == "4float" ){
    assert( (int)aF_RGBA.size() == nResX*nResY*4 );
    glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA,
                 nResX, nResY, 0, GL_RGBA, GL_FLOAT,
                 aF_RGBA.data());
  }
  glBindTexture(GL_TEXTURE_2D, 0);
}

void CGPUSampler::Draw() const {
  
  ::glPointSize(this->pointSize);
  this->Draw_Point();
  /////
  ::glLineWidth(3);
  this->Draw_Axis();
  ////
  ::glLineWidth(1);
  ::glColor3d(0,0,0);
  this->Draw_BoundingBox();
  
  if( id_tex_color > 0 && this->isDrawTex ){
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
  /////
  const CVector3& dx = x_axis;
  const CVector3& dz = z_axis;
  const CVector3& dy = Cross(dz,dx);
  ::glBegin(GL_POINTS);
  for(int iy=0;iy<nResY;++iy){
    for(int ix=0;ix<nResX;++ix){
      double lz = aZ[iy*nResX+ix];
      double lx = (ix+0.5)*lengrid;
      double ly = (iy+0.5)*lengrid;
      CVector3 vp = lx*dx+ly*dy+lz*dz + origin;
      myGlVertex(vp);
    }
  }
  ::glEnd();
}

std::vector<double> CGPUSampler::getGPos(int ix, int iy) const
{
  const CVector3& dx = x_axis;
  const CVector3& dz = z_axis;
  const CVector3& dy = Cross(dz,dx);
  double lz = aZ[iy*nResX+ix];
  double lx = (ix+0.5)*lengrid;
  double ly = (iy+0.5)*lengrid;
  CVector3 vp = lx*dx+ly*dy+lz*dz + origin;
  std::vector<double> res;
  res.push_back(vp.x);
  res.push_back(vp.y);
  res.push_back(vp.z);
  return res;
}
