/*
 * Copyright (c) 2019 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */

#include <vector>
#include <cstdlib>
#include <list>
#if defined(_WIN32) // windows
#  define NOMINMAX   // to remove min,max macro
#  include <windows.h>
#endif
#define GL_SILENCE_DEPRECATION
#include <GLFW/glfw3.h>

#include "delfem2/vec2.h"
#include "delfem2/geo_convhull2.h"
#include "delfem2/geo_gjk2.h"
#include "delfem2/glfw/viewer3.h"
#include "delfem2/glfw/util.h"

namespace dfm2 = delfem2;

void DrawPointsAndTheirConvexhull(
    const std::vector<dfm2::CVec2d> &vtxA_xy,
    const std::vector<unsigned int> &edgeA_vtx) {
  ::glPointSize(3);
  ::glBegin(GL_POINTS);
  for (const auto &xy: vtxA_xy) {
    ::glVertex2dv(xy.data());
  }
  ::glEnd();
  //
  ::glBegin(GL_LINE_LOOP);
  for (unsigned int iedge = 0; iedge < edgeA_vtx.size() / 2; iedge++) {
    const unsigned int i1 = edgeA_vtx[iedge * 2 + 0];
    const unsigned int i2 = edgeA_vtx[iedge * 2 + 1];
    ::glVertex2dv(vtxA_xy[i1].data());
    ::glVertex2dv(vtxA_xy[i2].data());
  }
  ::glEnd();
}


// ---------------------------------------

int main() {
  delfem2::glfw::CViewer3 viewer(1.5);
  delfem2::glfw::InitGLOld();
  viewer.OpenWindow();

  std::mt19937 rngeng(std::random_device{}());
  std::uniform_real_distribution<double> dist_m1p1(-1, +1);

  bool is_gjk = false;
  while (!glfwWindowShouldClose(viewer.window)) {
    std::vector<dfm2::CVec2d> vtxA_xy(10);
    for (auto &p: vtxA_xy) {
      p = {dist_m1p1(rngeng), dist_m1p1(rngeng)};
    }
    std::vector<unsigned int> edgeA_vtx;
    delfem2::ConvexHull2<delfem2::CVec2d>(edgeA_vtx, vtxA_xy);
    //
    std::vector<dfm2::CVec2d> vtxB_xy0(10);
    for (auto &p: vtxB_xy0) {
      p = {dist_m1p1(rngeng), dist_m1p1(rngeng)};
    }
    std::vector<unsigned int> edgeB_vtx;
    delfem2::ConvexHull2<delfem2::CVec2d>(edgeB_vtx, vtxB_xy0);
    //
    is_gjk = !is_gjk;
    if( is_gjk ){
      ::glfwSetWindowTitle(viewer.window, "GJK algorithm");
    }
    else{
      ::glfwSetWindowTitle(viewer.window, "SAT algorithm");
    }
    const double t_newsetting = glfwGetTime();
    while (true) {
      double t0 = glfwGetTime();
      if (t0 - t_newsetting > 5) { break; }
      double t1 = t0 * 0.1;
      std::vector<dfm2::CVec2d> vtxB_xy(vtxB_xy0.size());
      for (unsigned int ivtxB = 0; ivtxB < vtxB_xy0.size(); ++ivtxB) { // translate and rotate points
        double x0 = vtxB_xy0[ivtxB][0];
        double y0 = vtxB_xy0[ivtxB][1];
        double x1 = x0 + 2 * std::sin(3 * t1);
        double y1 = y0;
        double x2 = x1 * std::cos(t1) - y1 * std::sin(t1);
        double y2 = x1 * std::sin(t1) + y1 * std::cos(t1);
        vtxB_xy[ivtxB] = {x2, y2};
      }
      bool is_intersect = false;
      if( is_gjk ) {
        is_intersect = delfem2::IsIntersect_Points2_Points2_Gjk(vtxA_xy, vtxB_xy);
      }
      else {
        is_intersect = IsIntersect_Points2_Points2_Sat(vtxA_xy, vtxB_xy);
      }
      viewer.DrawBegin_oldGL();
      if (is_intersect) { ::glColor3d(1, 0, 0); }
      else { ::glColor3d(0, 0, 0); }
      DrawPointsAndTheirConvexhull(vtxA_xy, edgeA_vtx);
      DrawPointsAndTheirConvexhull(vtxB_xy, edgeB_vtx);
      glfwSwapBuffers(viewer.window);
      glfwPollEvents();
    }
  }
  glfwDestroyWindow(viewer.window);
  glfwTerminate();
  exit(EXIT_SUCCESS);
}
