/*
 * Copyright (c) 2019 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */

#ifndef DFM2_GLFW_VIEWER_H
#define DFM2_GLFW_VIEWER_H

#include <stdio.h>
#include <iostream>

#include "delfem2/opengl/glfw_cam.h" // for CNav3D_GLFW

// ------------------------------------------------------

namespace delfem2{
namespace opengl{

class CViewer_GLFW{
public:
  void Init_oldGL();
  void DrawBegin_oldGL();
  void DrawEnd_oldGL();
  
  void Init_newGL();
  
  virtual void mouse_press(const float src[3], const float dir[3]) {} // for function override. Do nothing here
  virtual void mouse_drag(const float src0[3], const float src1[3], const float dir[3]) {} // for function override. Do nothing here
  
public:
  GLFWwindow* window;
  CNav3D_GLFW nav;
};
  
}
}

#endif /* glfw_viewer_hpp */
