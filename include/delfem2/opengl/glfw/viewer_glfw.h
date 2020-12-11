/*
 * Copyright (c) 2019 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */

#ifndef DFM2_VIEWER_GLFW_H
#define DFM2_VIEWER_GLFW_H

#include "delfem2/cam3_m4q.h" // for CNav3D_GLFW
#include "delfem2/dfm2_inline.h"
#include <GLFW/glfw3.h>
#include <cstdio>
#include <iostream>

// ------------------------------------------------------

namespace delfem2{

/**
 * @brief class for 3D navigation (view rotation, scale & translation) for each GLFW window
 * @details this class should be compiled with OpenGL4.x
 */
class CNav3
{
public:
  CNav3(){
    ibutton = -1;
  }
  void MouseRay(
      float src[3],
      float dir[3],
      float asp,
      const float mMVP[16]) const
  {
    float mMVP_inv[16];
    ::delfem2::Inverse_Mat4(mMVP_inv, mMVP);
    const float ps[3] = {(float)mouse_x, (float)mouse_y, -1.0};
    const float pe[3] = {(float)mouse_x, (float)mouse_y, +1.0};
    float qs[3]; ::delfem2::Vec3_Vec3Mat4_AffineProjection(qs, ps, mMVP_inv);
    float qe[3]; ::delfem2::Vec3_Vec3Mat4_AffineProjection(qe, pe, mMVP_inv);
    src[0] = qs[0];
    src[1] = qs[1];
    src[2] = qs[2];
    dir[0] = qe[0] - qs[0];
    dir[1] = qe[1] - qs[1];
    dir[2] = qe[2] - qs[2];
  }

  void RayMouseMove(
      float src0[3],
      float src1[3],
      float dir0[3],
      float dir1[3],
      float asp,
      const float mMVP[16])
  {
    float mMVP_inv[16];
    ::delfem2::Inverse_Mat4(mMVP_inv, mMVP);
    const float p0s[3] = {(float)(mouse_x-dx), (float)(mouse_y-dy), -1.0};
    const float p0e[3] = {(float)(mouse_x-dx), (float)(mouse_y-dy), +1.0};
    const float p1s[3] = {(float)mouse_x, (float)mouse_y, -1.0};
    const float p1e[3] = {(float)mouse_x, (float)mouse_y, +1.0};
    float q0s[3]; ::delfem2::Vec3_Vec3Mat4_AffineProjection(q0s, p0s, mMVP_inv);
    float q0e[3]; ::delfem2::Vec3_Vec3Mat4_AffineProjection(q0e, p0e, mMVP_inv);
    float q1s[3]; ::delfem2::Vec3_Vec3Mat4_AffineProjection(q1s, p1s, mMVP_inv);
    float q1e[3]; ::delfem2::Vec3_Vec3Mat4_AffineProjection(q1e, p1e, mMVP_inv);
    //
    src0[0] = q0s[0];
    src0[1] = q0s[1];
    src0[2] = q0s[2];
    //
    src1[0] = q1s[0];
    src1[1] = q1s[1];
    src1[2] = q1s[2];
    //
    dir0[0] = q0e[0] - q0s[0];
    dir0[1] = q0e[1] - q0s[1];
    dir0[2] = q0e[2] - q0s[2];
    //
    dir1[0] = q1e[0] - q1s[0];
    dir1[1] = q1e[1] - q1s[1];
    dir1[2] = q1e[2] - q1s[2];
  }

public:
  int imodifier;
  int ibutton;
  double mouse_x, mouse_y;
  double dx, dy;
  double mouse_x_down, mouse_y_down;
};


namespace opengl{

class CViewer_GLFW{
public:
  void Init_oldGL();
  void DrawBegin_oldGL() const;
  void SwapBuffers() const;
  
  void Init_newGL();

  void ExitIfClosed() const;
  
  /**
   * @details for function override. Do nothing here
   */
  virtual void mouse_press(const float src[3], const float dir[3]) {}
  
  /**
   * @details for function override. Do nothing here
   */
  virtual void mouse_drag(const float src0[3], const float src1[3], const float dir[3]) {}

  virtual void key_press(int key, int mods) {
    if( key ==  GLFW_KEY_PAGE_UP  ){ camera.Scale(1.03); }
    if( key ==  GLFW_KEY_PAGE_DOWN ){ camera.Scale(1.0/1.03); }
    if( key == GLFW_KEY_BACKSPACE  ){ camera.is_pars = !camera.is_pars; }
    if( key ==  GLFW_KEY_HOME  ){ camera.fovy *= 1.03; }
    if( key ==  GLFW_KEY_END ){ camera.fovy *= 1.0/1.03; }
  }
  
  virtual void key_release(int key, int mods) {}
  
public:
  GLFWwindow* window = nullptr;
  CNav3 nav;
  delfem2::CCam3_OnAxisZplusLookOrigin<double> camera;
  double bgcolor[4] = {1,1,1,1};
  unsigned int width = 640;
  unsigned int height = 480;

};
  
}
}

#ifdef DFM2_HEADER_ONLY
# include "delfem2/opengl/glfw/viewer_glfw.cpp"
#endif

#endif /* glfw_viewer_hpp */
