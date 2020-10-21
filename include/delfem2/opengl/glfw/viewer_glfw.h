/*
 * Copyright (c) 2019 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */

#ifndef DFM2_VIEWER_GLFW_H
#define DFM2_VIEWER_GLFW_H

#include <stdio.h>
#include <iostream>

#include "delfem2/opengl/glfw/cam_glfw.h" // for CNav3D_GLFW
#include "delfem2/dfm2_inline.h"

// ------------------------------------------------------

namespace delfem2{
namespace opengl{

class CViewer_GLFW{
public:
  void Init_oldGL();
  void DrawBegin_oldGL() const;
  void DrawEnd_oldGL() const;
  
  void Init_newGL();
  
  /**
   * @details for function override. Do nothing here
   */
  virtual void mouse_press(const float src[3], const float dir[3]) {}
  
  /**
   * @details for function override. Do nothing here
   */
  virtual void mouse_drag(const float src0[3], const float src1[3], const float dir[3]) {}
  
  
  virtual void key_press(int key, int mods) {}
  
  virtual void key_release(int key, int mods) {}
  
public:
  GLFWwindow* window;
  CNav3D_GLFW nav;
  double bgcolor[4] = {1,1,1,1};
};
  
}
}

#ifdef DFM2_HEADER_ONLY
# include "delfem2/opengl/glfw/viewer_glfw.cpp"
#endif

#endif /* glfw_viewer_hpp */
