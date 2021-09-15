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
#  include <windows.h>
#endif
#define GL_SILENCE_DEPRECATION
#include <GLFW/glfw3.h>

#include "delfem2/glfw/viewer3.h"
#include "delfem2/glfw/util.h"
#include "delfem2/geoconvhull3.h"
#include "delfem2/vec3.h"


#if defined(_MSC_VER)
#  pragma warning( disable : 4100 )
#endif

namespace dfm2 = delfem2;

// ---------------------------------------

static void myGlVertex3d(
    const dfm2::CVec3d &v) {
  ::glVertex3d(v.x, v.y, v.z);
}

static void myGlVertex3d(
    unsigned int i,
    const std::vector<dfm2::CVec3d> &aV) {
  const dfm2::CVec3d &v = aV[i];
  ::glVertex3d(v.x, v.y, v.z);
}

int main() {
  delfem2::glfw::CViewer3 viewer;
  delfem2::glfw::InitGLOld();
  viewer.InitGL();

  viewer.camera.view_height = 1.5;

  std::vector<dfm2::CVec3d> vtx_xyz(100);
  std::vector<unsigned int> tri_vtx;

  double time_last_update = -2;
  while (!glfwWindowShouldClose(viewer.window)) {
    const double time_now = glfwGetTime();
    if (time_now - time_last_update > 1) {
      std::mt19937 rngeng(std::random_device{}());
      std::uniform_real_distribution<double> dist_m1p1(-1, +1);
      for (auto &p: vtx_xyz) {
        p.x = dist_m1p1(rngeng);
        p.y = dist_m1p1(rngeng);
        p.z = dist_m1p1(rngeng);
      }
      delfem2::ConvexHull<double>(tri_vtx, vtx_xyz);
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
      for (const auto &ixyz : vtx_xyz) {
        myGlVertex3d(ixyz);
      }
      ::glEnd();

      const size_t num_tri = tri_vtx.size() / 3;
      ::glColor4d(1, 1, 1, 0.8);
      ::glBegin(GL_TRIANGLES);
      for (unsigned int itri = 0; itri < num_tri; itri++) {
        const unsigned int i1 = tri_vtx[itri * 3 + 0];
        const unsigned int i2 = tri_vtx[itri * 3 + 1];
        const unsigned int i3 = tri_vtx[itri * 3 + 2];
        myGlVertex3d(i1, vtx_xyz);
        myGlVertex3d(i2, vtx_xyz);
        myGlVertex3d(i3, vtx_xyz);
      }
      ::glEnd();
      ::glColor3d(0, 0, 0);
      ::glBegin(GL_LINES);
      for (unsigned int itri = 0; itri < num_tri; itri++) {
        const unsigned int i1 = tri_vtx[itri * 3 + 0];
        const unsigned int i2 = tri_vtx[itri * 3 + 1];
        const unsigned int i3 = tri_vtx[itri * 3 + 2];
        myGlVertex3d(i1, vtx_xyz);
        myGlVertex3d(i2, vtx_xyz);
        myGlVertex3d(i2, vtx_xyz);
        myGlVertex3d(i3, vtx_xyz);
        myGlVertex3d(i3, vtx_xyz);
        myGlVertex3d(i1, vtx_xyz);
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
