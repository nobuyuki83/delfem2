/*
 * Copyright (c) 2019 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */

#include "delfem2/opengl/old/v2.h"

#include "delfem2/vec2.h"

#if defined(_WIN32) // windows
#define NOMINMAX   // to remove min,max macro
#include <windows.h>
#endif

#if defined(__APPLE__) && defined(__MACH__) // mac
#define GL_SILENCE_DEPRECATION
#include <OpenGL/gl.h>
#else
#include <GL/gl.h>
#endif


//----------------------------------------------------

void delfem2::opengl::myGlVertex2(
    int i,
    const std::vector<double> &vec) {
  ::glVertex3d(vec[i * 2], vec[i * 2 + 1], +0.0);
}

void delfem2::opengl::myGlVertex(
    unsigned int i,
    const std::vector<CVec2d> &aP) {
  ::glVertex3d(aP[i].x, aP[i].y, +0.0);
}

void delfem2::opengl::myGlVertex(const CVec2d &v) {
  ::glVertex2d(v.x, v.y);
}

//--------------------------------------------------------

void delfem2::opengl::drawPolyLine(
    const std::vector<CVec2d> &aP) {
  ::glBegin(GL_LINES);
  for (unsigned int ip = 0; ip < aP.size() - 1; ip++) {
    unsigned int jp = ip + 1;
    opengl::myGlVertex(ip, aP);
    opengl::myGlVertex(jp, aP);
  }
  ::glEnd();
  //
  ::glBegin(GL_POINTS);
  for (unsigned int ip = 0; ip < aP.size(); ip++) {
    opengl::myGlVertex(ip, aP);
  }
  ::glEnd();
}

void delfem2::opengl::drawPolyLine2D(
    const std::vector<CVec2d> &aP) {
  ::glBegin(GL_LINES);
  for (unsigned int ip = 0; ip < aP.size() - 1; ip++) {
    unsigned int jp = ip + 1;
    myGlVertex(ip, aP);
    myGlVertex(jp, aP);
  }
  ::glEnd();

  // ----------
  ::glBegin(GL_POINTS);
  for (unsigned int ip = 0; ip < aP.size(); ip++) {
    myGlVertex(ip, aP);
  }
  ::glEnd();
}

void delfem2::opengl::Draw_MeshTri(
    const std::vector<CVec2d> &aP,
    const std::vector<unsigned int> &aTri) {
  const size_t nTri = aTri.size() / 3;
  ::glBegin(GL_TRIANGLES);
  for (unsigned int itri = 0; itri < nTri; itri++) {
    const unsigned int i0 = aTri[itri * 3 + 0];
    const unsigned int i1 = aTri[itri * 3 + 1];
    const unsigned int i2 = aTri[itri * 3 + 2];
    const CVec2d &v0 = aP[i0];
    const CVec2d &v1 = aP[i1];
    const CVec2d &v2 = aP[i2];
    myGlVertex(v0);
    myGlVertex(v1);
    myGlVertex(v2);
  }
  ::glEnd();
}

void delfem2::opengl::Draw_MeshTri_Edge(
    const std::vector<CVec2d> &aP,
    const std::vector<unsigned int> &aTri) {
  //  const unsigned int nxys = (int)aXY.size()/2;
  ::glColor3d(0, 0, 0);
  ::glBegin(GL_LINES);
  const size_t nTri = aTri.size() / 3;
  for (unsigned int itri = 0; itri < nTri; itri++) {
    const unsigned int i0 = aTri[itri * 3 + 0];
    const unsigned int i1 = aTri[itri * 3 + 1];
    const unsigned int i2 = aTri[itri * 3 + 2];
    const CVec2d &v0 = aP[i0];
    const CVec2d &v1 = aP[i1];
    const CVec2d &v2 = aP[i2];
    myGlVertex(v0);
    myGlVertex(v1);
    myGlVertex(v1);
    myGlVertex(v2);
    myGlVertex(v2);
    myGlVertex(v0);
  }
  ::glEnd();
}

// -------------------------------------------------------