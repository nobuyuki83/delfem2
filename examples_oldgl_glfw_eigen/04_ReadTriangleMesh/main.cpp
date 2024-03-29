/*
 * Copyright (c) 2019 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */

#include <vector>
#include <cstdlib>
#include <fstream>
#include <tuple>
#if defined(_WIN32) // windows
#  define NOMINMAX   // to remove min,max macro
#  include <windows.h>  // should put before glfw3.h
#endif
#define GL_SILENCE_DEPRECATION
#include <GLFW/glfw3.h>

#include "delfem2/eigen/msh_io.h"
#include "delfem2/msh_io_obj.h"
#include "delfem2/glfw/viewer3.h"
#include "delfem2/glfw/util.h"
#include "delfem2/eigen_opengl/funcs.h"
#include "delfem2/opengl/old/mshuni.h"

namespace dfm2 = delfem2;

// ---------------------------------------

int main() {
  std::cout << "Available :SIMD Instructions: " << Eigen::SimdInstructionSetsInUse() << std::endl;

  const auto res = delfem2::eigen::ReadTriangleMeshObj(std::string(PATH_INPUT_DIR) + "/bunny_1k.obj");
  const auto& V = std::get<0>(res)*0.02;
  const auto& F = std::get<1>(res);
  std::cout << V.rows() << std::endl;
  std::cout << V.col(0).minCoeff() << " " << V.col(0).maxCoeff() << std::endl;

  delfem2::glfw::CViewer3 viewer;
  //
  delfem2::glfw::InitGLOld();
  viewer.OpenWindow();

  while (!glfwWindowShouldClose(viewer.window)) {
    viewer.DrawBegin_oldGL();
    ::glDisable(GL_LIGHTING);
    ::glColor3d(0.,0.,0.);
    delfem2::eigen_opengl::DrawMeshTri3_Edge_EigenMats(V,F);
    glfwSwapBuffers(viewer.window);
    glfwPollEvents();
  }
  glfwDestroyWindow(viewer.window);
  glfwTerminate();
  exit(EXIT_SUCCESS);
}
