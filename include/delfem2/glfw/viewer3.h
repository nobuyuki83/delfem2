/*
 * Copyright (c) 2019 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */

#ifndef DFM2_OPENGL_GLFW_VIEWER3_H
#define DFM2_OPENGL_GLFW_VIEWER3_H

#include "delfem2/cam3_m4q.h" // for CNav3D_GLFW
#include "delfem2/dfm2_inline.h"
#include "delfem2/glfw/mouseinput.h"
#include <GLFW/glfw3.h>
#include <cstdio>
#include <iostream>

// ------------------------------------------------------

namespace delfem2{

namespace glfw {

class CViewer3 {
public:
  void DrawBegin_oldGL() const;

  void SwapBuffers() const;

  void InitGL();

  void ExitIfClosed() const;

  /**
   * @details for function override. Do nothing here
   */
  virtual void mouse_press(const float src[3], const float dir[3]) {}

  virtual void mouse_release() {}

  /**
   * @details for function override. Do nothing here
   */
  virtual void mouse_drag(const float src0[3], const float src1[3], const float dir[3]) {}

  virtual void key_press(int key, int mods) {
    if (key == GLFW_KEY_PAGE_UP) { camera.Scale(1.03); }
    if (key == GLFW_KEY_PAGE_DOWN) { camera.Scale(1.0 / 1.03); }
    if (key == GLFW_KEY_BACKSPACE) { camera.is_pars = !camera.is_pars; }
    if (key == GLFW_KEY_HOME) { camera.fovy *= 1.03; }
    if (key == GLFW_KEY_END) { camera.fovy *= 1.0 / 1.03; }
  }

  virtual void key_release(int key, int mods) {}

  virtual void mouse_wheel(double yoffset){}

public:
  GLFWwindow *window = nullptr;
  CMouseInput nav;
  delfem2::CCam3_OnAxisZplusLookOrigin<double> camera;
  double bgcolor[4] = {1, 1, 1, 1};
  unsigned int width = 640;
  unsigned int height = 480;

};

} // opengl
} // delfem2

#ifdef DFM2_HEADER_ONLY
  #include "delfem2/glfw/viewer3.cpp"
#endif

#endif /* viewer3_hpp */
