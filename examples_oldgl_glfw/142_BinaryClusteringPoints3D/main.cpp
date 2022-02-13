/*
 * Copyright (c) 2019 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */

#include <vector>
#include <string>
#include <cstdlib>
#include <climits>
#include <random>
#if defined(_WIN32) // windows
#  define NOMINMAX   // to remove min,max macro
#  include <windows.h>  // this should come before glfw3.h
#endif
#define GL_SILENCE_DEPRECATION
#include <GLFW/glfw3.h>

#include "delfem2/mshmisc.h"
#include "delfem2/msh_io_ply.h"
#include "delfem2/msh_topology_uniform.h"
#include "delfem2/color.h"
#include "delfem2/msh_points.h"
#include "delfem2/clusterpoints.h"
#include "delfem2/glfw/viewer3.h"
#include "delfem2/glfw/util.h"
#include "delfem2/opengl/old/funcs.h"
#include "delfem2/opengl/old/mshuni.h"

namespace dfm2 = delfem2;

// -----------------------------

int main() {
  class CClusterData {
   public:
    std::vector<double> aXYZ; // center position of the cluster
    std::vector<double> aArea; // area of the cluster
    std::vector<double> aNorm; // normal of the cluster
    std::vector<unsigned int> psup_ind, psup; // connectivity of the cluster
    // below: data for visualization
    std::vector<float> aColor; // color of the cluster
    std::vector<unsigned int> map0c; // index of cluster for mesh points.
  };

  std::vector<CClusterData> aPointData;

  std::vector<unsigned int> aTri0;
  aPointData.resize(1);
  { // make level0
    {
      CClusterData &pd0 = aPointData[0];
      delfem2::Read_Ply(
          pd0.aXYZ, aTri0,
          std::filesystem::path(PATH_INPUT_DIR) / "bunny_2k.ply");
      dfm2::Normalize_Points3(pd0.aXYZ, 2.0);
    }
    { // make normal
      CClusterData &pd0 = aPointData[0];
      pd0.aNorm.resize(pd0.aXYZ.size());
      dfm2::Normal_MeshTri3D(
          pd0.aNorm.data(),
          pd0.aXYZ.data(), pd0.aXYZ.size() / 3,
          aTri0.data(), aTri0.size() / 3);
    }
    { // make area
      CClusterData &pd0 = aPointData[0];
      pd0.aArea.resize(pd0.aXYZ.size() / 3);
      dfm2::MassPoint_Tri3D(
		  pd0.aArea.data(),
		  1.0,                           
		  pd0.aXYZ.data(), pd0.aXYZ.size() / 3,
		  aTri0.data(), aTri0.size() / 3);
    }
    { // make psup
      CClusterData &pd0 = aPointData[0];
      dfm2::JArray_PSuP_MeshElem(
          pd0.psup_ind, pd0.psup,
          aTri0.data(), aTri0.size() / 3, 3,
          pd0.aXYZ.size() / 3);
    }
    {
      CClusterData &pd0 = aPointData[0];
      auto np0 = static_cast<unsigned int>(pd0.aXYZ.size() / 3);
      aPointData[0].map0c.resize(np0);
      for (unsigned int ip = 0; ip < np0; ++ip) {
        aPointData[0].map0c[ip] = ip;
      }
    }
  }

  for (unsigned int itr = 0; itr < 8; ++itr) {
    aPointData.resize(aPointData.size() + 1);
    const CClusterData &pd0 = aPointData[itr];
    CClusterData &pd1 = aPointData[itr + 1];
    std::vector<unsigned int> map01;
    dfm2::BinaryClustering_Points3d(
        pd1.aXYZ, pd1.aArea, pd1.aNorm, map01,
        pd0.aXYZ, pd0.aArea, pd0.aNorm, pd0.psup_ind, pd0.psup);
    dfm2::Clustering_Psup(
		pd1.psup_ind, pd1.psup,
		pd1.aXYZ.size() / 3,
		pd0.aXYZ.size() / 3,
		map01.data(), pd0.psup_ind.data(), pd0.psup.data());
    auto np0 = static_cast<unsigned int>(aPointData[0].aXYZ.size() / 3);
    pd1.map0c.resize(np0, UINT_MAX);
    for (unsigned int ip = 0; ip < np0; ++ip) {
      unsigned int ic0 = pd0.map0c[ip];
      assert(ic0 < map01.size());
      pd1.map0c[ip] = map01[ic0];
      assert(pd1.map0c[ip] < pd1.aXYZ.size() / 3);
    }

  }

  for (auto &pd : aPointData) {
    std::random_device rd;
    std::mt19937 eng(rd());
    std::uniform_real_distribution<double> dist(0, 1.0);
    const auto np = static_cast<unsigned int>(pd.aXYZ.size() / 3);
    pd.aColor.resize(np * 3);
    for (unsigned int ip = 0; ip < np; ++ip) {
      float *pc = pd.aColor.data() + ip * 3;
      dfm2::GetRGB_HSV(pc[0], pc[1], pc[2],
                       dist(eng), 1.0, 1.0);
    }
  }

  // -----------
  delfem2::glfw::CViewer3 viewer(1.5);
  //
  delfem2::glfw::InitGLOld();
  viewer.OpenWindow();
  //
  while (!glfwWindowShouldClose(viewer.window)) {
    for (const auto &pd: aPointData) {
      const CClusterData &dp0 = aPointData[0];
      for (unsigned int itr = 0; itr < 30; ++itr) {
        viewer.DrawBegin_oldGL();
        ::glBegin(GL_TRIANGLES);
        for (unsigned int it = 0; it < aTri0.size() / 3; ++it) {
          const unsigned int i0 = aTri0[it * 3 + 0];
          const unsigned int i1 = aTri0[it * 3 + 1];
          const unsigned int i2 = aTri0[it * 3 + 2];
          const unsigned int ic0 = pd.map0c[i0];
          assert(ic0 < pd.aColor.size() / 3);
          const unsigned int ic1 = pd.map0c[i1];
          assert(ic1 < pd.aColor.size() / 3);
          const unsigned int ic2 = pd.map0c[i2];
          assert(ic2 < pd.aColor.size() / 3);
          ::glColor3fv(pd.aColor.data() + ic0 * 3);
          ::glVertex3dv(dp0.aXYZ.data() + i0 * 3);
          ::glColor3fv(pd.aColor.data() + ic1 * 3);
          ::glVertex3dv(dp0.aXYZ.data() + i1 * 3);
          ::glColor3fv(pd.aColor.data() + ic2 * 3);
          ::glVertex3dv(dp0.aXYZ.data() + i2 * 3);
        }
        ::glEnd();
        ::glColor3d(0, 0, 0);
        dfm2::opengl::DrawMeshTri3D_Edge(
            aPointData[0].aXYZ.data(),
            aPointData[0].aXYZ.size() / 3,
            aTri0.data(),
            aTri0.size() / 3);
        viewer.SwapBuffers();
        glfwPollEvents();
      }
    }
    for (const auto &pd: aPointData) {
      for (unsigned int itr = 0; itr < 30; ++itr) {
        viewer.DrawBegin_oldGL();
        ::glColor3d(0, 0, 0);
        delfem2::opengl::DrawPoints3d_Psup(pd.aXYZ, pd.psup_ind, pd.psup);
        viewer.SwapBuffers();
        glfwPollEvents();
      }
    }
  }
  glfwDestroyWindow(viewer.window);
  glfwTerminate();
  exit(EXIT_SUCCESS);
}
