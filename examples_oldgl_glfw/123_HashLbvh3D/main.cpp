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

#include "delfem2/msh_points.h" // random uniform
#include "delfem2/srch_bvh.h"
#include "delfem2/srch_bv3_aabb.h"
#include "delfem2/srch_bv3_sphere.h"
#include "delfem2/glfw/viewer3.h"
#include "delfem2/glfw/util.h"
#include "delfem2/opengl/old/funcs.h"

namespace dfm2 = delfem2;

// ----------------------------------

void myGlutDisplay(
    const std::vector<double> &aXYZ,
    unsigned int ip_nearest,
    const double p0[3]) {
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

int main() {
  // input parameter for simulation
  std::vector<double> vtx_xyz; // 3d points
  std::vector<dfm2::CNodeBVH2> bvh_nodes;
  std::vector<dfm2::CBV3_Sphere<double>> bounding_volumes;
  double cur_time = 0.0;
  double p0[3];
  unsigned int idx_vtx_nearest;

  {
    const double min_xyz[3] = {-1, -1, -1};
    const double max_xyz[3] = {+1, +1, +1};
    dfm2::CBV3d_AABB bb(min_xyz, max_xyz);
    {
      const unsigned int N = 1000;
      vtx_xyz.resize(N * 3);
      dfm2::Points_RandomUniform(
          vtx_xyz.data(),
          N, 3, min_xyz, max_xyz);
      // create duplicated points for debugging purpose
      srand(3);
      for (int iip = 0; iip < 10; ++iip) { // hash collision
        const auto ip = static_cast<unsigned int>(N * (rand() / (RAND_MAX + 1.0)));
        assert(ip < N);
        const double x0 = vtx_xyz[ip * 3 + 0];
        const double y0 = vtx_xyz[ip * 3 + 1];
        const double z0 = vtx_xyz[ip * 3 + 2];
        for (int itr = 0; itr < 2; itr++) {
          vtx_xyz.insert(vtx_xyz.begin(), z0);
          vtx_xyz.insert(vtx_xyz.begin(), y0);
          vtx_xyz.insert(vtx_xyz.begin(), x0);
        }
      }
    }
    std::vector<unsigned int> sorted_idx;
    std::vector<std::uint32_t> sorted_morton_codes;
    dfm2::SortedMortenCode_Points3(
        sorted_idx, sorted_morton_codes,
        vtx_xyz, min_xyz, max_xyz);
    dfm2::BVHTopology_Morton(
        bvh_nodes,
        sorted_idx, sorted_morton_codes);
    dfm2::CLeafVolumeMaker_Point<dfm2::CBV3_Sphere<double>, double> lvm(
        vtx_xyz.data(), vtx_xyz.size() / 3);
    dfm2::BVH_BuildBVHGeometry(
        bounding_volumes,
        0, bvh_nodes, lvm);
    {
      dfm2::Check_MortonCode_Sort(sorted_idx, sorted_morton_codes, vtx_xyz, bb.bbmin, bb.bbmax);
      dfm2::Check_MortonCode_RangeSplit(sorted_morton_codes);
      dfm2::Check_BVH(bvh_nodes, vtx_xyz.size() / 3);
    }
  }

  dfm2::glfw::CViewer3 viewer(1.5);
  //
  dfm2::glfw::InitGLOld();
  viewer.OpenWindow();
  delfem2::opengl::setSomeLighting();

  while (!glfwWindowShouldClose(viewer.window)) {
    cur_time += 0.001;
    p0[0] = 3.0 * sin(cur_time * 1) - 1;
    p0[1] = 3.0 * sin(cur_time * 2) - 1;
    p0[2] = 3.0 * sin(cur_time * 3) - 1;
    // -----------
    double dist = -1;
    idx_vtx_nearest = 0;
    dfm2::BVH_IndPoint_NearestPoint(
        idx_vtx_nearest, dist,
        p0,
        0, bvh_nodes,
        bounding_volumes);
    // -----------
    viewer.DrawBegin_oldGL();
    myGlutDisplay(vtx_xyz, idx_vtx_nearest, p0);
    viewer.SwapBuffers();
    glfwPollEvents();
  }

  glfwDestroyWindow(viewer.window);
  glfwTerminate();
  exit(EXIT_SUCCESS);
}


