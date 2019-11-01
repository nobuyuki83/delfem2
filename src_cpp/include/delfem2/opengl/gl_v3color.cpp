/*
 * Copyright (c) 2019 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */

#include "delfem2/vec3.h"

// -------

#if defined(__APPLE__)
#include <OpenGL/gl.h>
#elif defined(_MSC_VER) // windows
#include <windows.h>
#include <GL/gl.h>
#else
#include <GL/gl.h>
#endif

#include "delfem2/opengl/gl_v3color.h"
#include "delfem2/opengl/gl2_color.h"
#include "delfem2/opengl/gl2_v23.h"

// -----------------------------

void DrawQuad_ScalarQ1
(const CVector3& p0, const CVector3& p1, const CVector3& p2, const CVector3& p3,
 double v0, double v1, double v2, double v3,
 const std::vector<std::pair<double, CColor> >& colorMap)
{
  ::glBegin(GL_QUADS);
  {
    heatmap(v0, colorMap);
    CVector3 n0; UnitNormal(n0,  p0, p1, p3);
    myGlNormal(n0);
    myGlVertex(p0);
  }
  {
    heatmap(v1, colorMap);
    CVector3 n1; UnitNormal(n1,  p0, p1, p2);
    myGlNormal(n1);
    myGlVertex(p1);
  }
  {
    heatmap(v2, colorMap);
    CVector3 n2; UnitNormal(n2,  p1, p2, p3);
    myGlNormal(n2);
    myGlVertex(p2);
  }
  {
    heatmap(v3, colorMap);
    CVector3 n3; UnitNormal(n3,  p2, p3, p0);
    myGlNormal(n3);
    myGlVertex(p3);
  }
  ::glEnd();
}
