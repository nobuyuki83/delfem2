/*
 * Copyright (c) 2019 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */

#ifndef DFM2_GLFW_VIEWER3_H
#define DFM2_GLFW_VIEWER3_H

#include <cstdio>
#include <iostream>
#include <memory>  // unique_ptr
#include <functional>

#include "delfem2/cam_modelview.h"
#include "delfem2/cam_projection.h"
#include "delfem2/dfm2_inline.h"
#include "delfem2/glfw/mouseinput.h"

// ------------------------------------------------------

struct GLFWwindow;

namespace delfem2::glfw {

class CViewer3 {
 public:
  explicit CViewer3(double view_height=1) :
  projection(){
    projection = std::make_unique<Projection_LookOriginFromZplus<double>>(
        view_height, false);
    view_rotation = std::make_unique<ModelView_Trackball>();
  }
  
  void DrawBegin_oldGL() const;

  void SwapBuffers() const;

  virtual void InitGL();

  void ExitIfClosed() const;

  /**
   * @details for function override. Do nothing here
   */
  virtual void mouse_press(
      [[maybe_unused]] const float src[3],
      [[maybe_unused]] const float dir[3]) {}

  virtual void mouse_release() {}

  /**
   * @details for function override. Do nothing here
   */
  virtual void mouse_drag(
      [[maybe_unused]] const float src0[3],
      [[maybe_unused]] const float src1[3],
      [[maybe_unused]] const float dir[3]) {}

  virtual void key_press(
      [[maybe_unused]] int key,
      [[maybe_unused]] int mods) {}

  virtual void key_release(
      [[maybe_unused]] int key,
      [[maybe_unused]] int mods) {}

  virtual void mouse_wheel(
      [[maybe_unused]] double yoffset) {}

  [[nodiscard]] std::array<float,16> GetModelViewMatrix() const {
    CMat4f mv = view_rotation->GetMatrix();
    CMat4f mt = CMat4f::Translate(trans[0], trans[1], trans[2]);
    return (mt * mv).GetStlArray();
  }

  [[nodiscard]] std::array<float,16> GetProjectionMatrix() const;

  virtual void CursorPosition(double xpos, double ypos);

 public:
  GLFWwindow *window = nullptr;
  CMouseInput nav;
  double scale = 1.0;
  double trans[3] = {0,0,0};
  double bgcolor[4] = {1, 1, 1, 1};
  unsigned int width = 640;
  unsigned int height = 480;
  std::string window_title = "LearnOpenGL";
  std::unique_ptr<Projection<double>> projection;
  std::unique_ptr<delfem2::ModelView_Trackball> view_rotation;
  std::vector< std::function< void(int,int) >> keypress_callbacks;
};

} // namespace delfem2


#ifndef DFM2_STATIC_LIBRARY
#  include "delfem2/glfw/viewer3.cpp"
#endif

#endif // DFM2_GLFW_VIEWER3_H
