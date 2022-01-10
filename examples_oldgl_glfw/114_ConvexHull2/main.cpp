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

#include "delfem2/glfw/viewer3.h"
#include "delfem2/glfw/util.h"
#include "delfem2/vec2.h"
#include "delfem2/geo_convhull2.h"

namespace dfm2 = delfem2;

// ---------------------------------------

int main() {
  delfem2::glfw::CViewer3 viewer(1.5);
  delfem2::glfw::InitGLOld();
  viewer.OpenWindow();

  std::mt19937 rngeng(std::random_device{}());
  std::uniform_real_distribution<double> dist_m1p1(-1, +1);

  std::vector<dfm2::CVec2d> vtx_xy(100);
  std::vector<unsigned int> edge_vtx;

  double time_last_update = -2;
  while (!glfwWindowShouldClose(viewer.window)) {
    const double time_now = glfwGetTime();
    if (time_now - time_last_update > 1) {
      for (auto &p: vtx_xy) {
        p.x = dist_m1p1(rngeng);
        p.y = dist_m1p1(rngeng);
      }
      delfem2::ConvexHull2<delfem2::CVec2d>(edge_vtx, vtx_xy);
      time_last_update = time_now;
    }

    viewer.DrawBegin_oldGL();

    ::glEnable(GL_BLEND);
    ::glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    ::glEnable(GL_CULL_FACE);
    ::glCullFace(GL_BACK);
    {
      ::glColor3d(0, 0, 0);
      ::glPointSize(3);
      ::glBegin(GL_POINTS);
      for (const auto &xy : vtx_xy) {
        ::glVertex2dv(xy.data());
      }
      ::glEnd();

      ::glColor3d(0, 0, 0);
      ::glBegin(GL_LINE_LOOP);
      for (unsigned int iedge = 0; iedge < edge_vtx.size() / 2; iedge++) {
        const unsigned int i1 = edge_vtx[iedge * 2 + 0];
        const unsigned int i2 = edge_vtx[iedge * 2 + 1];
        ::glVertex2dv(vtx_xy[i1].data());
        ::glVertex2dv(vtx_xy[i2].data());
      }
      ::glEnd();
    }

    glfwSwapBuffers(viewer.window);
    glfwPollEvents();
  }
  glfwDestroyWindow(viewer.window);
  glfwTerminate();
  exit(EXIT_SUCCESS);
}
