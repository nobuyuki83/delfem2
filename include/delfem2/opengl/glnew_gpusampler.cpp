/*
 * Copyright (c) 2019 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */


#include <iostream>
#include <cmath>
#include <stack>
#include <fstream>
#include <sstream>
#include <cstring>
#include <cstdlib>

#include "glad/glad.h" // gl3.0+
#if defined(__APPLE__) && defined(__MACH__) // Mac
  #include <OpenGL/gl.h>
#elif defined(_WIN32) // windows
  #include <windows.h>
  #include <GL/gl.h>
#else
  #include <GL/gl.h>
#endif

#include "delfem2/vec3.h"
#include "delfem2/opengl/glnew_gpusampler.h"

namespace dfm2 = delfem2;

// --------------------------------------------------------

double DotX(const double* p0, const double* p1, int ndof){
  double v=0;
  for(int i=0;i<ndof;++i){ v += p0[i]*p1[i]; }
  return v;
}

void ScaleX(double* p0, int n, double s)
{
  for(int i=0;i<n;++i){ p0[i] *= s; }
}

void Normalize(double* p0, int n)
{
  const double ss = DotX(p0,p0,n);
  ScaleX(p0,n,1.0/sqrt(ss));
}

//PURPOSE:      For square matrices. This is column major for OpenGL
inline void MultiplyMatrices4by4OpenGL_FLOAT(
    float *result,
    const float *m1,
    const float *m2)
{
  result[ 0]=m1[0]*m2[0]+  m1[4]*m2[1]+  m1[8]*m2[2]+  m1[12]*m2[3];
  result[ 4]=m1[0]*m2[4]+  m1[4]*m2[5]+  m1[8]*m2[6]+  m1[12]*m2[7];
  result[ 8]=m1[0]*m2[8]+  m1[4]*m2[9]+  m1[8]*m2[10]+  m1[12]*m2[11];
  result[12]=m1[0]*m2[12]+   m1[4]*m2[13]+  m1[8]*m2[14]+  m1[12]*m2[15];

  result[ 1]=m1[1]*m2[0]+  m1[5]*m2[1]+  m1[9]*m2[2]+  m1[13]*m2[3];
  result[ 5]=m1[1]*m2[4]+  m1[5]*m2[5]+  m1[9]*m2[6]+  m1[13]*m2[7];
  result[ 9]=m1[1]*m2[8]+  m1[5]*m2[9]+  m1[9]*m2[10]+  m1[13]*m2[11];
  result[13]=m1[1]*m2[12]+  m1[5]*m2[13]+  m1[9]*m2[14]+  m1[13]*m2[15];

  result[ 2]=m1[2]*m2[0]+  m1[6]*m2[1]+  m1[10]*m2[2]+  m1[14]*m2[3];
  result[ 6]=m1[2]*m2[4]+  m1[6]*m2[5]+  m1[10]*m2[6]+  m1[14]*m2[7];
  result[10]=m1[2]*m2[8]+  m1[6]*m2[9]+  m1[10]*m2[10]+  m1[14]*m2[11];
  result[14]=m1[2]*m2[12]+   m1[6]*m2[13]+  m1[10]*m2[14]+  m1[14]*m2[15];

  result[ 3]=m1[3]*m2[0]+  m1[7]*m2[1]+  m1[11]*m2[2]+  m1[15]*m2[3];
  result[ 7]=m1[3]*m2[4]+  m1[7]*m2[5]+  m1[11]*m2[6]+  m1[15]*m2[7];
  result[11]=m1[3]*m2[8]+   m1[7]*m2[9]+  m1[11]*m2[10]+  m1[15]*m2[11];
  result[15]=m1[3]*m2[12]+   m1[7]*m2[13]+  m1[11]*m2[14]+  m1[15]*m2[15];
}

// --------------------------------------------

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
  const int npix = nw*nh;
  /////
  if( isDepth ){ aZ.resize(npix,0); }
  else{ aZ.clear(); }
  // -------------
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

void CGPUSampler::Matrix_MVP
 (float mMV[16],
  float mP[16]) const
{
  {
    const double* ax = this->x_axis;
    const double* az = this->z_axis;
    double ay[3]; Cross3D(ay, az, ax);
    const double o[3] = { Dot3D(ax,origin), Dot3D(ay,origin), Dot3D(az,origin) };
    mMV[ 0] = ax[0];  mMV[ 1] = ay[0];  mMV[ 2] = az[0];  mMV[ 3] = 0;
    mMV[ 4] = ax[1];  mMV[ 5] = ay[1];  mMV[ 6] = az[1];  mMV[ 7] = 0;
    mMV[ 8] = ax[2];  mMV[ 9] = ay[2];  mMV[10] = az[2];  mMV[11] = 0;
    mMV[12] = -o[0];  mMV[13] = -o[1];  mMV[14] = -o[2];  mMV[15] = 1;
  }
  {
    double l = 0.0;
    double r = +lengrid*nResX;
    double b = 0.0;
    double t = +lengrid*nResY;
    double n = -z_range;
    double f = 0;
    /*
    mP[0*4+0] = 2.0/(r-l);
    mP[0*4+1] = 0.0;
    mP[0*4+2] = 0.0;
    mP[0*4+3] = -(l+r)/(r-l);
    mP[1*4+0] = 0.0;
    mP[1*4+1] = 2.0/(t-b);
    mP[1*4+2] = 0.0;
    mP[1*4+3] = -(t+b)/(t-b);
    mP[2*4+0] = 0.0;
    mP[2*4+1] = 0.0;
    mP[2*4+2] = 2.0/(n-f);
    mP[2*4+3] = -(n+f)/(n-f);
    mP[3*4+0] = 0.0;
    mP[3*4+1] = 0.0;
    mP[3*4+2] = 0.0;
    mP[3*4+3] = 1.0;
     */
    mP[0*4+0] = 2.0/(r-l);
    mP[1*4+0] = 0.0;
    mP[2*4+0] = 0.0;
    mP[3*4+0] = -(l+r)/(r-l);
    mP[0*4+1] = 0.0;
    mP[1*4+1] = 2.0/(t-b);
    mP[2*4+1] = 0.0;
    mP[3*4+1] = -(t+b)/(t-b);
    mP[0*4+2] = 0.0;
    mP[1*4+2] = 0.0;
    mP[2*4+2] = 2.0/(n-f);
    mP[3*4+2] = -(n+f)/(n-f);
    mP[0*4+3] = 0.0;
    mP[1*4+3] = 0.0;
    mP[2*4+3] = 0.0;
    mP[3*4+3] = 1.0;
  }
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
}

void CGPUSampler::End()
{
  ::glBindFramebuffer(GL_FRAMEBUFFER, 0);

  { // get depth from texture
#ifdef EMSCRIPTEN
    std::cout << "the function \"glGetTexImage\" is not supported in emscripten" << std::endl;
      ::glViewport(view[0], view[1], view[2], view[3]);
    return;
#endif
    std::vector<float> aDepth;
    aDepth.resize(nResX*nResY);
    ::glBindTexture(GL_TEXTURE_2D, id_tex_depth);
    ::glGetTexImage(GL_TEXTURE_2D, 0,
        GL_DEPTH_COMPONENT, GL_FLOAT,
        (void*)aDepth.data());
    std::vector<double> aXYZ(nResX*nResY*3);
    const double* ax = this->x_axis;
    const double* az = this->z_axis;
    double ay[3]; Cross3D(ay, az, ax);
    for(int iy=0;iy<nResY;++iy){
      for(int ix=0;ix<nResX;++ix){
        int ip = iy*nResX+ix;
        double lz = -aDepth[ip]*this->z_range;
        double lx = (ix+0.5)*lengrid;
        double ly = (iy+0.5)*lengrid;
        aXYZ[ip*3+0] = origin[0] + lx*ax[0] + ly*ay[0] + lz*az[0];
        aXYZ[ip*3+1] = origin[1] + lx*ax[1] + ly*ay[1] + lz*az[1];
        aXYZ[ip*3+2] = origin[2] + lx*ax[2] + ly*ay[2] + lz*az[2];
      }
    }
    shdr2.Initialize(aXYZ);
  }
  
  ::glViewport(view[0], view[1], view[2], view[3]);
}

void CGPUSampler::InitGL() {
  ::glEnable(GL_TEXTURE_2D);
  ::glActiveTexture(GL_TEXTURE0);

  { // initialize color texture
    // create to render to
    if( id_tex_color > 0 ){ glDeleteTextures(1, &id_tex_color); }
    ::glGenTextures(1, &id_tex_color);
    ::glBindTexture(GL_TEXTURE_2D, id_tex_color);
    // define size and format of level 0
    ::glTexImage2D(GL_TEXTURE_2D, 0,
                   GL_RGBA, nResX, nResY, 0,
                   GL_RGBA, GL_UNSIGNED_BYTE, nullptr);
//                   GL_RGBA, GL_FLOAT, nullptr);
    // set the filtering so we don't need mips
    ::glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
    ::glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
    ::glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
    ::glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);
  }

  { // depth texture
    // create to render to
    if( id_tex_depth > 0 ){ glDeleteTextures(1, &id_tex_depth); }
    ::glGenTextures(1, &id_tex_depth);
    ::glBindTexture(GL_TEXTURE_2D, id_tex_depth);
    // define size and format of level 0
    ::glTexImage2D(GL_TEXTURE_2D, 0,
                   GL_DEPTH_COMPONENT32F, nResX, nResY, 0,
                   GL_DEPTH_COMPONENT, GL_FLOAT, nullptr);
    // set the filtering so we don't need mips
    ::glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
    ::glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
    ::glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
    ::glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);
  }
  ::glBindTexture(GL_TEXTURE_2D,0);
  {
    ::glGenFramebuffers(1,&id_framebuffer);
    ::glBindFramebuffer(GL_FRAMEBUFFER, id_framebuffer);
    // attach the texture as the first color attachment
    ::glFramebufferTexture2D(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_TEXTURE_2D,
                             id_tex_color, 0);
    ::glFramebufferTexture2D(GL_FRAMEBUFFER, GL_DEPTH_ATTACHMENT, GL_TEXTURE_2D,
                             id_tex_depth, 0);
    // Always check that our framebuffer is ok
    GLenum status = glCheckFramebufferStatus(GL_FRAMEBUFFER) ;
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
    }
  }
  //
  { // draw grid
    this->shdr0.Compile();
    double xmin = 0.0;
    double xmax = lengrid*nResX;
    double ymin = 0.0;
    double ymax = lengrid*nResY;
    double zmin = 0.0;
    double zmax = -z_range;
    std::vector<double> aPos3d = {
        xmin, ymin, zmin,
        xmin, ymin, zmax,
        xmin, ymax, zmin,
        xmin, ymax, zmax,
        xmax, ymin, zmin,
        xmax, ymin, zmax,
        xmax, ymax, zmin,
        xmax, ymax, zmax,
    };
    std::vector<unsigned int> aTri = {
        0,  1,
        1,  3,
        2,  3,
        0,  2,
        4,  5,
        5,  7,
        6,  7,
        4,  6,
        0,  4,
        1,  5,
        2,  6,
        3,  7,
    };
    shdr0.Initialize(aPos3d, aTri);
  }
  // -----
  { // draw texture
    shdr1.Compile();
    // --------------
    const CVector3& dx = x_axis;
    const CVector3& dy = Cross(z_axis,dx);
    const double lx = lengrid*nResX;
    const double ly = lengrid*nResY;
    CVector3 p0 = origin;
    CVector3 p1 = origin + lx*dx;
    CVector3 p2 = origin + lx*dx + ly*dy;
    CVector3 p3 = origin + ly*dy;
    std::vector<double> aPos3d = {
        p0.x, p0.y, p0.z,
        p1.x, p1.y, p1.z,
        p2.x, p2.y, p2.z,
        p3.x, p3.y, p3.z,
    };
    std::vector<unsigned int> aTri = {
        0, 1, 2,
        0, 2, 3,
    };
    std::vector<double> aTex2d = {
        0.0, 0.0,
        1.0, 0.0,
        1.0, 1.0,
        0.0, 1.0
    };
    shdr1.Initialize(aPos3d, aTri, aTex2d);
  }

  {
    shdr2.Compile();
  }
}

void CGPUSampler::Draw(float mP[16], float mV[16]) const
{
 ::glLineWidth(5);
  float mM[16];
  {
    const double* ax = this->x_axis;
    const double* az = this->z_axis;
    double ay[3]; Cross3D(ay, az, ax);
    const double* o = this->origin;
    mM[ 0] = ax[0];  mM[ 1] = ax[1];  mM[ 2] = ax[2];  mM[ 3] = 0;
    mM[ 4] = ay[0];  mM[ 5] = ay[1];  mM[ 6] = ay[2];  mM[ 7] = 0;
    mM[ 8] = az[0];  mM[ 9] = az[1];  mM[10] = az[2];  mM[11] = 0;
    mM[12] = +o[0];  mM[13] = +o[1];  mM[14] = +o[2];  mM[15] = 1;
  }
  float mMV[16];
  for(int i=0;i<4;++i){
    for(int j=0;j<4;++j) {
      mMV[i*4+j] = 0;
      for(int k=0;k<4;++k){
        mMV[i*4+j] += mM[i*4+k] *  mV[k*4+j];
      }
    }
  }
  shdr0.Draw(mP,mMV);
  shdr2.Draw(mP,mV);
  glEnable(GL_TEXTURE_2D);
  glActiveTexture(0);
  glBindTexture(GL_TEXTURE_2D, this->id_tex_color);
  shdr1.Draw(mP,mV);
  glBindTexture(GL_TEXTURE_2D, 0);
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
