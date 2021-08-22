/*
 * Copyright (c) 2019 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */

// ----------------
#ifdef EMSCRIPTEN
  #include <GLFW/glfw3.h>
#else
  #include "glad/glad.h" // gl3.0+
#endif
// ----------------
#ifdef _WIN32
  #include <windows.h>
#endif
// ----------------
#if defined(__APPLE__) && defined(__MACH__) // Mac
  #include <OpenGL/gl.h>
#else
  #include <GL/gl.h>
#endif
// ----------------
#include "delfem2/mat4.h"
#include "delfem2/opengl/new/r2tgln.h"
#include <stack>

// --------------------------------------------

DFM2_INLINE void delfem2::opengl::CRender2Tex_DrawNewGL::SetDepth(
    const delfem2::opengl::CRender2Tex& r2t)
{
  if( r2t.aZ.size() != r2t.width*r2t.height ){ return; }
  //
  unsigned int nx = r2t.width;
  unsigned int ny = r2t.height;
  std::vector<double> aXYZ(nx*ny*3);
  for(unsigned int iy=0;iy<ny;++iy){
    for(unsigned int ix=0;ix<nx;++ix){
      const int ip = iy*nx+ix;
      double q0[3] = {
          (ix+0.5)/nx*2.0-1.0,
          (iy+0.5)/ny*2.0-1.0,
          r2t.aZ[ip]*2.0-1.0 };
      aXYZ[ip*3+0] = q0[0];
      aXYZ[ip*3+1] = q0[1];
      aXYZ[ip*3+2] = q0[2];
    }
  }
  shdr2.Initialize(aXYZ);
}

DFM2_INLINE void delfem2::opengl::CRender2Tex_DrawNewGL::InitGL()
{
  //
  { // draw grid
    this->shdr0.Compile();
    std::vector<double> aPos3d = {
        -1, -1, -1,
        -1, -1, +1,
        -1, +1, -1,
        -1, +1, +1,
        +1, -1, -1,
        +1, -1, +1,
        +1, +1, -1,
        +1, +1, +1,
    };
    std::vector<unsigned int> aLine = {
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
    shdr0.Initialize(aPos3d, 3, aLine, GL_LINES);
  }
  // -----
  { // draw texture
    shdr1.Compile();
    // --------------
    std::vector<double> aPos3d = {
        -1, -1, -1,
        +1, -1, -1,
        +1, +1, -1,
        -1, +1, -1
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
    shdr1.setCoords(aPos3d,3);
    shdr1.setTexCoords(aTex2d);
    shdr1.setElement( aTri, GL_TRIANGLES);
  }
  {
    shdr2.Compile();
  }
}

DFM2_INLINE void delfem2::opengl::CRender2Tex_DrawNewGL::Draw(
    const delfem2::opengl::CRender2Tex& r2t,
    float mP0[16],
    float mMV0[16]) const
{
  double mMVP[16]; delfem2::MatMat4(mMVP,r2t.mat_modelview_colmajor,r2t.mat_projection_colmajor);
  double mMVPinv[16]; delfem2::Inverse_Mat4(mMVPinv,mMVP);
  float mMVP1[16]; delfem2::MatMat4(mMVP1,mMVPinv,mMV0);
  shdr0.Draw(mP0,mMVP1);
  shdr2.Draw(mP0,mMVP1);
  glEnable(GL_TEXTURE_2D);
  glActiveTexture(0);
  glBindTexture(GL_TEXTURE_2D, r2t.id_tex_color);
  shdr1.Draw(mP0,mMVP1);
  glBindTexture(GL_TEXTURE_2D, 0);
}
