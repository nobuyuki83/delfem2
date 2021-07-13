/*
 * Copyright (c) 2019-2021 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */

#ifndef DFM2_GLFW_VIEWER2_H
#define DFM2_GLFW_VIEWER2_H

#include "delfem2/dfm2_inline.h"
#include "delfem2/glfw/mouseinput.h"
//
#include <cstdio>
#include <iostream>

// ------------------------------------------------------

struct GLFWwindow;

namespace delfem2{

namespace glfw {

class CViewer2 {
public:
  void InitGL();

  void DrawBegin_oldGL() const;

  void SwapBuffers() const;

  void ExitIfClosed() const;

  /**
   * @details for function override. Do nothing here
   */
  virtual void mouse_press(const float src[2]) {}

  /**
   * @details for function override. Do nothing here
   */
  virtual void mouse_drag(const float src0[2], const float src1[2]) {}

  virtual void key_press(int key, int mods) {}

  virtual void key_release(int key, int mods) {}

  void Mat4_MVP_OpenGL(float mMV[16], float mP[16], float asp) const;

public:
  GLFWwindow *window = nullptr;
  CMouseInput nav;
  std::string title;
  double bgcolor[4] = {1, 1, 1, 1};
  unsigned int width = 640;
  unsigned int height = 480;
  //
  float view_height = 1.f;
  float scale = 1.f;
  float trans[2] = {0.f, 0.f};
};

} // opengl
} // delfem2

#ifdef DFM2_HEADER_ONLY
  #include "delfem2/glfw/viewer2.cpp"
#endif

#endif /* viewer2_hpp */
