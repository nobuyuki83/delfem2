/*
 * Copyright (c) 2019 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */

#include "delfem2/opengl/old/color.h"
#include "delfem2/opengl/old/v3q.h"
#include "delfem2/opengl/old/colorv3.h"
#include "delfem2/vec3.h"

// -------

#if defined(_WIN32) // windows
  #define NOMINMAX   // to remove min,max macro
  #include <windows.h>
#endif

#if defined(__APPLE__) && defined(__MACH__)
  #define GL_SILENCE_DEPRECATION
  #include <OpenGL/gl.h>
#else
  #include <GL/gl.h>
#endif

// -----------------------------

DFM2_INLINE void DrawQuad_ScalarQ1
(const delfem2::CVec3d& p0,
 const delfem2::CVec3d& p1,
 const delfem2::CVec3d& p2,
 const delfem2::CVec3d& p3,
 double v0, double v1, double v2, double v3,
 const std::vector<std::pair<double, delfem2::CColor> >& colorMap)
{
  ::glBegin(GL_QUADS);
  {
    delfem2::opengl::heatmap(v0, colorMap);
    delfem2::CVec3d n0; UnitNormal(n0,  p0, p1, p3);
    delfem2::opengl::myGlNormal(n0);
    delfem2::opengl::myGlVertex(p0);
  }
  {
    delfem2::opengl::heatmap(v1, colorMap);
    delfem2::CVec3d n1; UnitNormal(n1,  p0, p1, p2);
    delfem2::opengl::myGlNormal(n1);
    delfem2::opengl::myGlVertex(p1);
  }
  {
    delfem2::opengl::heatmap(v2, colorMap);
    delfem2::CVec3d n2; UnitNormal(n2,  p1, p2, p3);
    delfem2::opengl::myGlNormal(n2);
    delfem2::opengl::myGlVertex(p2);
  }
  {
    delfem2::opengl::heatmap(v3, colorMap);
    delfem2::CVec3d n3; UnitNormal(n3,  p2, p3, p0);
    delfem2::opengl::myGlNormal(n3);
    delfem2::opengl::myGlVertex(p3);
  }
  ::glEnd();
}

