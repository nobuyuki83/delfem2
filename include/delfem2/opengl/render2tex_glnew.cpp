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

// ----------------
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
#include "delfem2/opengl/render2tex_glnew.h"

namespace dfm2 = delfem2;

// --------------------------------------------

void dfm2::opengl::CRender2Tex_DrawNewGL::SetDepth()
{
  std::vector<float> aZ;
  CRender2Tex::ExtractFromTexture_Depth(aZ);
  assert( aZ.size() == nResX*nResY );
  std::vector<double> aXYZ(nResX*nResY*3);
  const double* ax = this->x_axis;
  const double* az = this->z_axis;
  double ay[3]; dfm2::Cross3(ay, az, ax);
  for(int iy=0;iy<nResY;++iy){
    for(int ix=0;ix<nResX;++ix){
      int ip = iy*nResX+ix;
      double lz = -aZ[ip]*this->z_range;
      double lx = (ix+0.5)*lengrid;
      double ly = (iy+0.5)*lengrid;
      aXYZ[ip*3+0] = origin[0] + lx*ax[0] + ly*ay[0] + lz*az[0];
      aXYZ[ip*3+1] = origin[1] + lx*ax[1] + ly*ay[1] + lz*az[1];
      aXYZ[ip*3+2] = origin[2] + lx*ax[2] + ly*ay[2] + lz*az[2];
    }
  }
  shdr2.Initialize(aXYZ);
}

void dfm2::opengl::CRender2Tex_DrawNewGL::InitGL() {
  CRender2Tex::InitGL();
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
    const dfm2::CVec3d& dx = x_axis;
    const dfm2::CVec3d& dy = dfm2::CVec3d(z_axis)^dx;
    const double lx = lengrid*nResX;
    const double ly = lengrid*nResY;
    dfm2::CVec3d p0 = origin;
    dfm2::CVec3d p1 = p0 + lx*dx;
    dfm2::CVec3d p2 = p0 + lx*dx + ly*dy;
    dfm2::CVec3d p3 = p0 + ly*dy;
    std::vector<double> aPos3d = {
        p0.x(), p0.y(), p0.z(),
        p1.x(), p1.y(), p1.z(),
        p2.x(), p2.y(), p2.z(),
        p3.x(), p3.y(), p3.z(),
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

void dfm2::opengl::CRender2Tex_DrawNewGL::Draw(float mP[16], float mV[16]) const
{
 ::glLineWidth(5);
  float mM[16];
  {
    const double* ax = this->x_axis;
    const double* az = this->z_axis;
    double ay[3]; dfm2::Cross3(ay, az, ax);
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
