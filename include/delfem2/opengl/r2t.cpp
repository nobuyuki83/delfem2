/*
 * Copyright (c) 2019 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */

#include <iostream>

#ifdef _WIN32
#  include <windows.h>  // should be included before glfw3.h
#endif

#ifdef EMSCRIPTEN
#  include <GLFW/glfw3.h>
#else
#  include "glad/glad.h" // gl3.0+
#endif

#if defined(__APPLE__) && defined(__MACH__) // Mac
#  include <OpenGL/gl.h>
#else
#  include <GL/gl.h>
#endif

#include "delfem2/opengl/r2t.h"

// ---------------------------------------------

namespace delfem2 {
namespace opengl {
namespace r2t {

DFM2_INLINE double DotX(const double *p0, const double *p1, int ndof) {
  double v = 0;
  for (int i = 0; i < ndof; ++i) { v += p0[i] * p1[i]; }
  return v;
}

DFM2_INLINE void ScaleX(double *p0, int n, double s) {
  for (int i = 0; i < n; ++i) { p0[i] *= s; }
}

DFM2_INLINE void NormalizeX(double *p0, int n) {
  const double ss = DotX(p0, p0, n);
  ScaleX(p0, n, 1.0 / sqrt(ss));
}

template<typename T>
static T MyDot3(const T a[3], const T b[3]) {
  return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}
#ifdef DFM2_STATIC_LIBRARY
template float MyDot3(const float a[3], const float b[3]);
template double MyDot3(const double a[3], const double b[3]);
#endif

template<typename T>
void MyCross3(T r[3], const T v1[3], const T v2[3]) {
  r[0] = v1[1] * v2[2] - v2[1] * v1[2];
  r[1] = v1[2] * v2[0] - v2[2] * v1[0];
  r[2] = v1[0] * v2[1] - v2[0] * v1[1];
}
#ifdef DFM2_STATIC_LIBRARY
template void MyCross3(float r[3], const float v1[3], const float v2[3]);
template void MyCross3(double r[3], const double v1[3], const double v2[3]);
#endif

}  // namespace r2t
}  // namespace opengl
}  // namespace delfem2

DFM2_INLINE void delfem2::Mat4_OrthongoalProjection_AffineTrans(
    double mMV[16],
    double mP[16],
    const double origin[3],
    const double az[3],
    const double ax[3],
    unsigned int nResX,
    unsigned int nResY,
    double lengrid,
    double z_range) {
  namespace lcl = delfem2::opengl::r2t;
  { // global to local
    double ay[3];
    lcl::MyCross3(ay, az, ax);
    const double o[3] = {lcl::MyDot3(ax, origin), lcl::MyDot3(ay, origin), lcl::MyDot3(az, origin)};
    mMV[0] = ax[0];
    mMV[4] = ax[1];
    mMV[8] = ax[2];
    mMV[12] = -o[0];
    mMV[1] = ay[0];
    mMV[5] = ay[1];
    mMV[9] = ay[2];
    mMV[13] = -o[1];
    mMV[2] = az[0];
    mMV[6] = az[1];
    mMV[10] = az[2];
    mMV[14] = -o[2];
    mMV[3] = 0;
    mMV[7] = 0;
    mMV[11] = 0;
    mMV[15] = 1;
  }
  { //
    double l = 0.0;
    double r = +lengrid * nResX;
    double b = 0.0;
    double t = +lengrid * nResY;
    double n = -z_range;
    double f = 0;
    mP[0 * 4 + 0] = 2.0 / (r - l);
    mP[1 * 4 + 0] = 0.0;
    mP[2 * 4 + 0] = 0.0;
    mP[3 * 4 + 0] = -(l + r) / (r - l);
    mP[0 * 4 + 1] = 0.0;
    mP[1 * 4 + 1] = 2.0 / (t - b);
    mP[2 * 4 + 1] = 0.0;
    mP[3 * 4 + 1] = -(t + b) / (t - b);
    mP[0 * 4 + 2] = 0.0;
    mP[1 * 4 + 2] = 0.0;
    mP[2 * 4 + 2] = 2.0 / (n - f);
    mP[3 * 4 + 2] = -(n + f) / (n - f);
    mP[0 * 4 + 3] = 0.0;
    mP[1 * 4 + 3] = 0.0;
    mP[2 * 4 + 3] = 0.0;
    mP[3 * 4 + 3] = 1.0;
  }
}

// --------------------------------------------

DFM2_INLINE void delfem2::opengl::CRender2Tex::Start() {
  ::glGetIntegerv(GL_VIEWPORT, view); // current viewport
  ::glViewport(0, 0,
               static_cast<int>(width),
               static_cast<int>(height));
  ::glBindFramebuffer(GL_FRAMEBUFFER, id_framebuffer);
  ::glBindTexture(GL_TEXTURE_2D, 0);  // unbind texture
}

DFM2_INLINE void delfem2::opengl::CRender2Tex::End() {
  ::glBindFramebuffer(GL_FRAMEBUFFER, 0);
  ::glViewport(view[0], view[1], view[2], view[3]);
  this->CopyToCPU_Depth();
  this->CopyToCPU_RGBA8UI();
  this->CopyToCPU_RGBA32F();
}

DFM2_INLINE void delfem2::opengl::CRender2Tex::CopyToCPU_Depth() {

#ifdef EMSCRIPTEN
  std::cout << "In delfem2::opengl::CRender2Tex::CopyToCPU_Depth()" << std::endl;
  std::cout << "the function \"glGetTexImage\" is not supported in emscripten" << std::endl;
  return;
#endif
  //
  aZ.resize(width * height);
  ::glBindFramebuffer(GL_FRAMEBUFFER, id_framebuffer);
  ::glBindTexture(GL_TEXTURE_2D, id_tex_depth);
  ::glGetTexImage(GL_TEXTURE_2D, 0,
                  GL_DEPTH_COMPONENT, GL_FLOAT,
                  (void *) aZ.data());
  ::glBindTexture(GL_TEXTURE_2D, 0);
  ::glBindFramebuffer(GL_FRAMEBUFFER, 0);
}

DFM2_INLINE void delfem2::opengl::CRender2Tex::CopyToCPU_RGBA8UI() {
#ifdef EMSCRIPTEN
  std::cout << "the function \"glGetTexImage\" is not supported in emscripten" << std::endl;
  return;
#endif
  aRGBA_8ui.resize(width * height * 4);
  ::glBindTexture(GL_TEXTURE_2D, id_tex_color);
  ::glGetTexImage(GL_TEXTURE_2D, 0,
                  GL_RGBA, GL_UNSIGNED_BYTE,
                  (void *) aRGBA_8ui.data());
  ::glBindTexture(GL_TEXTURE_2D, 0);
  ::glBindFramebuffer(GL_FRAMEBUFFER, 0);
}

DFM2_INLINE void delfem2::opengl::CRender2Tex::CopyToCPU_RGBA32F() {

#ifdef EMSCRIPTEN
  std::cout << "the function \"glGetTexImage\" is not supported in emscripten" << std::endl;
  return;
#endif

  aRGBA_32f.resize(width * height * 4);
  std::vector<float> aF_RGBA(width * height * 4);
  ::glBindTexture(GL_TEXTURE_2D, id_tex_color);
  ::glGetTexImage(GL_TEXTURE_2D, 0,
                  GL_RGBA, GL_FLOAT,
                  (void *) aRGBA_32f.data());
  ::glBindTexture(GL_TEXTURE_2D, 0);
  ::glBindFramebuffer(GL_FRAMEBUFFER, 0);
}

DFM2_INLINE void delfem2::opengl::CRender2Tex::InitGL() {
  ::glEnable(GL_TEXTURE_2D);
  ::glActiveTexture(GL_TEXTURE0);

  { // initialize color texture
    // create and bind texture
    if (id_tex_color > 0) { glDeleteTextures(1, &id_tex_color); }
    ::glGenTextures(1, &id_tex_color);
    ::glBindTexture(GL_TEXTURE_2D, id_tex_color);
    // define size and format of level 0
    if (is_rgba_8ui) {
      ::glTexImage2D(GL_TEXTURE_2D, 0,
                     GL_RGBA, static_cast<int>(width), static_cast<int>(height), 0,
                     GL_RGBA, GL_UNSIGNED_BYTE, nullptr);
    } else {
      ::glTexImage2D(GL_TEXTURE_2D, 0,
                     GL_RGBA, static_cast<int>(width), static_cast<int>(height), 0,
                     GL_RGBA, GL_FLOAT, nullptr);
    }
    // set the filtering so we don't need mips
    ::glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
    ::glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
    ::glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
    ::glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);
  }

  { // depth texture
    // create to render to
    if (id_tex_depth > 0) { glDeleteTextures(1, &id_tex_depth); }
    ::glGenTextures(1, &id_tex_depth);
    ::glBindTexture(GL_TEXTURE_2D, id_tex_depth);
    // define size and format of level 0
    ::glTexImage2D(GL_TEXTURE_2D, 0,
                   GL_DEPTH_COMPONENT32F,
                   static_cast<int>(width),
                   static_cast<int>(height), 0,
                   GL_DEPTH_COMPONENT, GL_FLOAT, nullptr);
    // set the filtering so we don't need mips
    ::glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
    ::glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
    ::glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
    ::glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);
  }
  ::glBindTexture(GL_TEXTURE_2D, 0);
  {
    ::glGenFramebuffers(1, &id_framebuffer);
    ::glBindFramebuffer(GL_FRAMEBUFFER, id_framebuffer);
    // attach the texture as the first color attachment
    ::glFramebufferTexture2D(
        GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_TEXTURE_2D,
        id_tex_color, 0);
    ::glFramebufferTexture2D(
        GL_FRAMEBUFFER, GL_DEPTH_ATTACHMENT, GL_TEXTURE_2D,
        id_tex_depth, 0);
    // Always check that our framebuffer is ok
    GLenum status = glCheckFramebufferStatus(GL_FRAMEBUFFER);
    if (status != GL_FRAMEBUFFER_COMPLETE) {
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
    ::glBindFramebuffer(GL_FRAMEBUFFER, 0);
  }

  if (aRGBA_8ui.size() == width * height * 4) {
    ::glBindTexture(GL_TEXTURE_2D, id_tex_color);
    // define size and format of level 0
    ::glTexImage2D(GL_TEXTURE_2D, 0,
                   GL_RGBA, static_cast<int>(width), static_cast<int>(height), 0,
                   GL_RGBA, GL_UNSIGNED_BYTE, aRGBA_8ui.data());
  }
}

DFM2_INLINE void delfem2::opengl::CRender2Tex::BoundingBox3(
    double *pmin,
    double *pmax) const {
  /*
  if( aZ.size() != nResX*nResY ){ return; }
  for(unsigned int ix=0;ix<nResX;++ix){
    for(unsigned int iy=0;iy<nResY;++iy){
      CVec3d vp;
      {
        const CVec3d& dx = x_axis;
        const CVec3d& dz = z_axis;
        const CVec3d& dy = Cross(dz,dx);
        double lz = aZ[iy*nResX+ix];
        double lx = (ix+0.5)*lengrid;
        double ly = (iy+0.5)*lengrid;
        vp = lx*dx+ly*dy+lz*dz + CVec3d(origin);
        if( -lz > z_range*0.99 ) continue;
      }
      // ----------
      if( pmin[0] > pmax[0] ){
        pmin[0] = pmax[0] = vp.x();
        pmin[1] = pmax[1] = vp.y();
        pmin[2] = pmax[2] = vp.z();
        continue;
      }
      const double x0 = vp.x();
      if( x0 < pmin[0] ){ pmin[0] = x0; }
      if( x0 > pmax[0] ){ pmax[0] = x0; }
      const double y0 = vp.y();
      if( y0 < pmin[1] ){ pmin[1] = y0; }
      if( y0 > pmax[1] ){ pmax[1] = y0; }
      const double z0 = vp.z();
      if( z0 < pmin[2] ){ pmin[2] = z0; }
      if( z0 > pmax[2] ){ pmax[2] = z0; }
    }
  }
   */
}

DFM2_INLINE bool delfem2::opengl::GetProjectedPoint(
    CVec3d &p0,
    CVec3d &n0,
    const CVec3d &ps,
    const CRender2Tex &smplr) {
  double mMVPG[16];
  smplr.GetMVPG(mMVPG);
  double mMVPGinv[16];
  Inverse_Mat4(mMVPGinv, mMVPG);
  double pg[3];
  Vec3_Vec3Mat4_AffineProjection(pg, ps.p, mMVPG);
  const unsigned int nx = smplr.width;
  const unsigned int ny = smplr.height;
  const int ix0 = (int) floor(pg[0]);
  const int iy0 = (int) floor(pg[1]);
  const int ix1 = ix0 + 1;
  const int iy1 = iy0 + 1;
  if (ix0 < 0 || ix0 >= (int) nx) { return false; }
  if (ix1 < 0 || ix1 >= (int) nx) { return false; }
  if (iy0 < 0 || iy0 >= (int) ny) { return false; }
  if (iy1 < 0 || iy1 >= (int) ny) { return false; }
  if (smplr.aZ[iy0 * nx + ix0] > 0.99) return false;
  if (smplr.aZ[iy0 * nx + ix1] > 0.99) return false;
  if (smplr.aZ[iy1 * nx + ix0] > 0.99) return false;
  if (smplr.aZ[iy1 * nx + ix1] > 0.99) return false;
  const CVec3d p00(ix0, iy0, smplr.aZ[iy0 * nx + ix0]);
  const CVec3d p01(ix0, iy1, smplr.aZ[iy1 * nx + ix0]);
  const CVec3d p10(ix1, iy0, smplr.aZ[iy0 * nx + ix1]);
  const CVec3d p11(ix1, iy1, smplr.aZ[iy1 * nx + ix1]);
  const double rx = pg[0] - ix0;
  const double ry = pg[1] - iy0;
  CVec3d p3 = (1 - rx) * (1 - ry) * p00 + rx * (1 - ry) * p10 + (1 - rx) * ry * p01 + rx * ry * p11;
  CVec3d dpx = (ry - 1) * p00 + (1 - ry) * p10 - ry * p01 + ry * p11;
  CVec3d dpy = (rx - 1) * p00 - rx * p10 + (1 - rx) * p01 + rx * p11;
  CVec3d n3 = Cross(dpx, dpy);
  //
  Vec3_Vec3Mat4_AffineProjection(p0.p, p3.p, mMVPGinv);
  CVec3d p4;
  Vec3_Vec3Mat4_AffineProjection(p4.p, (p3 + n3).p, mMVPGinv);
  n0 = (p4 - p0).normalized();
  return true;
}
