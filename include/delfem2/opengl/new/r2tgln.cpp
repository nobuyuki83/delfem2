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
  if( r2t.aDepth.size() != r2t.width*r2t.height ){ return; }
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
          1.-r2t.aDepth[ip]*2.0 };
      aXYZ[ip*3+0] = q0[0];
      aXYZ[ip*3+1] = q0[1];
      aXYZ[ip*3+2] = q0[2];
    }
  }
  drawer_projected_points.SetRawArray(aXYZ.data(),aXYZ.size()/3, 3);
}

DFM2_INLINE void delfem2::opengl::CRender2Tex_DrawNewGL::InitGL()
{
  //
  { // draw grid
    this->drawer_view_frustrum.Compile();
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
    drawer_view_frustrum.Initialize(aPos3d, 3, aLine, GL_LINES);
  }
  // -----
  // draw texture
  drawer_projected_image.InitGL();
  drawer_projected_points.InitGL();
}

DFM2_INLINE void delfem2::opengl::CRender2Tex_DrawNewGL::Draw(
    const delfem2::opengl::CRender2Tex& r2t,
    const float mat4_projection[16],
    const float mat4_modelview[16]) const
{
  namespace dfm2 = delfem2;
  const dfm2::CMat4d mvp0 = dfm2::CMat4d(r2t.mat_modelview) * dfm2::CMat4d(r2t.mat_projection);
  const dfm2::CMat4f mv1 = (dfm2::CMat4d(mat4_modelview) * mvp0.Inverse()).cast<float>();
  drawer_view_frustrum.Draw(mat4_projection, mv1.data());
#ifndef EMSCRIPTEN
  glPointSize(this->pointSize);
#endif
  drawer_projected_points.Draw(GL_POINTS, mat4_projection, mv1.data());
  glEnable(GL_TEXTURE_2D);
  glActiveTexture(0);
  glBindTexture(GL_TEXTURE_2D, r2t.id_tex_color);
  drawer_projected_image.Draw(mat4_projection, mv1.data());
  glBindTexture(GL_TEXTURE_2D, 0);
}
