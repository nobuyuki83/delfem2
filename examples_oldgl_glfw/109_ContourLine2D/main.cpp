/*
 * Copyright (c) 2019 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */

#include <vector>
#include <cstdlib>
#if defined(_WIN32) // windows
#  define NOMINMAX   // to remove min,max macro
#  include <windows.h>  // this should come before glfw3.h
#endif
#define GL_SILENCE_DEPRECATION
#include <GLFW/glfw3.h>

#include "delfem2/points.h"
#include "delfem2/mshuni.h"
#include "delfem2/slice.h"
#include "delfem2/mshprimitive.h"
#include "delfem2/color.h"
#include "delfem2/glfw/viewer3.h"
#include "delfem2/glfw/util.h"
#include "delfem2/opengl/old/mshuni.h"
#include "delfem2/opengl/old/color.h"

// ---------------------------

void Draw(
    const std::vector<double> &vtx_xy,
    const std::vector<unsigned int> &tri_vtx,
    const std::vector<double> &vtx_val,
    const std::vector<delfem2::CSegInfo> &aSeg) {
  ::glDisable(GL_LIGHTING);
  ::glLineWidth(1);
  delfem2::opengl::DrawMeshTri2D_Edge(
      vtx_xy.data(), vtx_xy.size() / 2,
      tri_vtx.data(), tri_vtx.size() / 3);

  std::vector<std::pair<double, delfem2::CColor> > colorMap;
  delfem2::ColorMap_BlueCyanGreenYellowRed(
      colorMap,
      -1, 1);
  delfem2::opengl::DrawMeshTri2D_ScalarP1(
      vtx_xy.data(), vtx_xy.size() / 2,
      tri_vtx.data(), tri_vtx.size() / 3,
      vtx_val.data(), 1,
      colorMap);
  ::glLineWidth(5);
  ::glBegin(GL_LINES);
  ::glColor3d(0, 0, 0);
  for (auto &iseg : aSeg) {
    double pA[2], pB[2];
    iseg.Pos2D(pA, pB,
               vtx_xy.data(), tri_vtx.data());
    ::glVertex2dv(pA);
    ::glVertex2dv(pB);
  }
  ::glEnd();

}

void Initialize(
    std::vector<double> &vtx_xy,
    std::vector<unsigned int> &tri_vtx,
    std::vector<double> &vtx_val) {
  { // make mesh
    const int ndiv = 16;
    std::vector<unsigned int> aQuad;
    delfem2::MeshQuad2D_Grid(
        vtx_xy, aQuad,
        ndiv, ndiv);
    delfem2::convert2Tri_Quad(
        tri_vtx,
        aQuad);
    delfem2::Translate_Points2(
        vtx_xy,
        -ndiv * 0.5, -ndiv * 0.5);
    delfem2::Scale_PointsX(
        vtx_xy,
        1.0 / ndiv);
  }

  { // make value
    const size_t np = vtx_xy.size() / 2;
    vtx_val.resize(np);
    for (size_t ip = 0; ip < np; ++ip) {
      double x0 = vtx_xy[ip * 2 + 0];
      double y0 = vtx_xy[ip * 2 + 1];
      vtx_val[ip] = sqrt(x0 * x0 + y0 * y0) * 4.0 - 1.5;
    }
  }

}

int main() {
  std::vector<double> vtx_xy;
  std::vector<unsigned int> tri_vtx;
  std::vector<double> vtx_val;
  Initialize(vtx_xy, tri_vtx, vtx_val);
  std::vector<delfem2::CSegInfo> segments;
  delfem2::glfw::CViewer3 viewer;
  //
  delfem2::glfw::InitGLOld();
  viewer.InitGL();
  while (!glfwWindowShouldClose(viewer.window)) {
    {
      static int iframe = 0;
      double thresthold = 0.9 * sin(iframe * 0.001);
      segments.clear();
      delfem2::AddContour(
          segments,
          thresthold,
          tri_vtx.data(), tri_vtx.size() / 3,
          vtx_val.data());
      iframe += 1;
    }
    // ------------------
    viewer.DrawBegin_oldGL();
    Draw(vtx_xy, tri_vtx, vtx_val, segments);
    viewer.SwapBuffers();
    glfwPollEvents();
  }
  glfwDestroyWindow(viewer.window);
  glfwTerminate();
  exit(EXIT_SUCCESS);
}
