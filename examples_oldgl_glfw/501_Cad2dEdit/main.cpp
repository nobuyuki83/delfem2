/*
 * Copyright (c) 2019-2021 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */

#include <cmath>
#include <stack>
#if defined(_WIN32) // windows
#  define NOMINMAX   // to remove min,max macro
#  include <windows.h>  // this should come before glfw3.h
#endif
#define GL_SILENCE_DEPRECATION
#include <GLFW/glfw3.h>

#include "delfem2/cad2_dtri2.h"
#include "delfem2/glfw/viewer3.h"
#include "delfem2/glfw/util.h"
#include "delfem2/opengl/old/funcs.h"
#include "delfem2/opengl/old/cad2dtriv2.h"

#ifndef M_PI
#  define M_PI 3.141592653589793
#endif
// -------------------------------------

int main() {
  class CCAD2DViewer : public delfem2::glfw::CViewer3 {
   public:
    CCAD2DViewer() : CViewer3(1.5) {
      const double poly[8] = {-1, -1, +1, -1, +1, +1, -1, +1};
      cad.AddPolygon(std::vector<double>(poly, poly + 8));
    }
    //
    void mouse_press(const float src[3], [[maybe_unused]] const float dir[3]) override {
      cad.Pick(src[0], src[1], 1.5);
    }
    void mouse_drag(const float src0[3], const float src1[3], [[maybe_unused]] const float dir[3]) override {
      cad.DragPicked(src1[0], src1[1], src0[0], src0[1]);
    }
    //
    void Draw() {
      DrawBegin_oldGL();
      delfem2::opengl::Draw_CCad2D(cad);
      SwapBuffers();
    }
    delfem2::CCad2D cad;
  };
  // --------------------
  CCAD2DViewer viewer;
  delfem2::glfw::InitGLOld();
  viewer.InitGL();
  delfem2::opengl::setSomeLighting();
  // --------------------
  while (!glfwWindowShouldClose(viewer.window)) {
    viewer.Draw();
    glfwPollEvents();
  }
  glfwDestroyWindow(viewer.window);
  glfwTerminate();
  exit(EXIT_SUCCESS);
}

