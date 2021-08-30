/*
 * Copyright (c) 2019 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */

#include <vector>
#include <algorithm>
#if defined(_WIN32) // windows
#  define NOMINMAX   // to remove min,max macro
#  include <windows.h>  // this should come before glfw3.h
#endif
#define GL_SILENCE_DEPRECATION
#include <GLFW/glfw3.h>

#include "delfem2/points.h" // random uniform
#include "delfem2/srchbvh.h"
#include "delfem2/srchbv3aabb.h"
#include "delfem2/srchbv3sphere.h"
#include "delfem2/glfw/viewer3.h"
#include "delfem2/glfw/util.h"
#include "delfem2/opengl/old/funcs.h"

namespace dfm2 = delfem2;

// ------------------------------------
// input parameter for simulation
std::vector<double> aXYZ; // 3d points
std::vector<dfm2::CNodeBVH2> aNodeBVH;
std::vector<dfm2::CBV3_Sphere<double>> aAABB;
double cur_time = 0.0;
double p0[3];
unsigned int ip_nearest;

// ----------------------------------------

void myGlutDisplay() {
  ::glDisable(GL_LIGHTING);
  ::glPointSize(2);
  ::glBegin(GL_POINTS);
  ::glColor3d(0, 0, 0);
  for (size_t ip = 0; ip < aXYZ.size() / 3; ++ip) {
    ::glVertex3d(aXYZ[ip * 3 + 0], aXYZ[ip * 3 + 1], aXYZ[ip * 3 + 2]);
  }
  ::glEnd();
  //
  ::glPointSize(4);
  ::glBegin(GL_POINTS);
  ::glColor3d(1, 0, 0);
  ::glVertex3dv(p0);
  ::glEnd();
  //
  ::glColor3d(1, 0, 0);
  ::glBegin(GL_LINES);
  ::glVertex3dv(p0);
  ::glVertex3dv(aXYZ.data() + ip_nearest * 3);
  ::glEnd();
}

int main(
    [[maybe_unused]] int argc,
    [[maybe_unused]] char *argv[]) {
  {
    const double min_xyz[3] = {-1, -1, -1};
    const double max_xyz[3] = {+1, +1, +1};
    dfm2::CBV3d_AABB bb(min_xyz, max_xyz);
    {
      const unsigned int N = 1000;
      aXYZ.resize(N * 3);
      dfm2::Points_RandomUniform(aXYZ.data(),
                                 N, 3, min_xyz, max_xyz);
      // create duplicated points for debugging purpose
      srand(3);
      for (int iip = 0; iip < 10; ++iip) { // hash collision
        const unsigned int ip = N * (rand() / (RAND_MAX + 1.0));
        assert(ip < N);
        const double x0 = aXYZ[ip * 3 + 0];
        const double y0 = aXYZ[ip * 3 + 1];
        const double z0 = aXYZ[ip * 3 + 2];
        for (int itr = 0; itr < 2; itr++) {
          aXYZ.insert(aXYZ.begin(), z0);
          aXYZ.insert(aXYZ.begin(), y0);
          aXYZ.insert(aXYZ.begin(), x0);
        }
      }
    }
    std::vector<unsigned int> aSortedId;
    std::vector<std::uint32_t> aSortedMc;
    dfm2::SortedMortenCode_Points3(aSortedId, aSortedMc,
                                   aXYZ, min_xyz, max_xyz);
    dfm2::BVHTopology_Morton(aNodeBVH,
                             aSortedId, aSortedMc);
    dfm2::CLeafVolumeMaker_Point<dfm2::CBV3_Sphere<double>, double> lvm(aXYZ.data(), aXYZ.size() / 3);
    dfm2::BVH_BuildBVHGeometry(aAABB,
                               0, aNodeBVH, lvm);
    {
      dfm2::Check_MortonCode_Sort(aSortedId, aSortedMc, aXYZ, bb.bbmin, bb.bbmax);
      dfm2::Check_MortonCode_RangeSplit(aSortedMc);
      dfm2::Check_BVH(aNodeBVH, aXYZ.size() / 3);
    }
  }

  dfm2::glfw::CViewer3 viewer;
  dfm2::glfw::InitGLOld();
  viewer.InitGL();
  viewer.camera.view_height = 1.5;
  delfem2::opengl::setSomeLighting();

  while (!glfwWindowShouldClose(viewer.window)) {
    cur_time += 0.001;
    p0[0] = 3.0 * sin(cur_time * 1) - 1;
    p0[1] = 3.0 * sin(cur_time * 2) - 1;
    p0[2] = 3.0 * sin(cur_time * 3) - 1;
    // -----------
    double dist = -1;
    ip_nearest = 0;
    dfm2::BVH_IndPoint_NearestPoint(
        ip_nearest, dist,
        p0,
        0, aNodeBVH,
        aAABB);
    // -----------
    viewer.DrawBegin_oldGL();
    myGlutDisplay();
    viewer.SwapBuffers();
    glfwPollEvents();
  }

  glfwDestroyWindow(viewer.window);
  glfwTerminate();
  exit(EXIT_SUCCESS);
}


