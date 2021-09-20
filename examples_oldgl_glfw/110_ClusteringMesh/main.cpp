/*
 * Copyright (c) 2019-2020 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */

#include <vector>
#include <cstdlib>
#include <random>
#include <filesystem>
#if defined(_WIN32) // windows
#  define NOMINMAX   // to remove min,max macro
#  include <windows.h>  // this should come before glfw3.h
#endif
#define GL_SILENCE_DEPRECATION
#include <GLFW/glfw3.h>

#include "delfem2/dijkstra.h"
#include "delfem2/points.h"
#include "delfem2/msh_io_ply.h"
#include "delfem2/mshuni.h"
#include "delfem2/color.h"
#include "delfem2/glfw/viewer3.h"
#include "delfem2/glfw/util.h"
#include "delfem2/opengl/old/funcs.h"
#include "delfem2/opengl/old/mshuni.h"
#include "delfem2/opengl/old/color.h"

namespace dfm2 = delfem2;

// ---------------------------

int main() {
  std::vector<double> vec_xyz;
  std::vector<unsigned int> vec_tri;

  delfem2::Read_Ply(
      vec_xyz, vec_tri,
      std::filesystem::path(PATH_INPUT_DIR) / "arm_16k.ply");
  delfem2::Normalize_Points3(vec_xyz);
  std::vector<unsigned int> tri_adjtri;
  ElSuEl_MeshElem(
      tri_adjtri,
      vec_tri.data(), vec_tri.size() / 3,
      delfem2::MESHELEM_TRI,
      vec_xyz.size() / 3);

  // ------
  std::mt19937 rdeng(std::random_device{}());
  std::uniform_int_distribution<unsigned int> ncluster_gen(1, 100);

  // above: data preparation
  // -----------------------
  // below: view

  delfem2::glfw::CViewer3 viewer;
  viewer.projection.view_height = 0.5;
  
  delfem2::glfw::InitGLOld();
  viewer.InitGL();
  delfem2::opengl::setSomeLighting();

  while (!glfwWindowShouldClose(viewer.window)) {
    const unsigned int ncluster = ncluster_gen(rdeng);
    std::vector<std::pair<int, dfm2::CColor> > aColor;
    for (unsigned int ic = 0; ic < ncluster; ++ic) {
      dfm2::CColor c;
      c.setRandomVividColor();
      aColor.emplace_back(2, c);
    }
    std::vector<unsigned int> vec_triangle_flag;
    dfm2::MeshClustering(
        vec_triangle_flag, ncluster, tri_adjtri,
        vec_tri.size() / 3);
    //
    for (unsigned int iframe = 0; iframe < 30; ++iframe) {
      viewer.DrawBegin_oldGL();
      ::glDisable(GL_LIGHTING);
      ::glColor3d(0, 0, 0);
      delfem2::opengl::DrawMeshTri3D_Edge(
          vec_xyz, vec_tri);
      delfem2::opengl::DrawMeshTri3DFlag_FaceNorm(
          vec_xyz, vec_tri, vec_triangle_flag, aColor);
      glfwSwapBuffers(viewer.window);
      glfwPollEvents();
    }
  }

  glfwDestroyWindow(viewer.window);
  glfwTerminate();
  exit(EXIT_SUCCESS);
}
