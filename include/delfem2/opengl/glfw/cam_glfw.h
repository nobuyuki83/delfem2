/*
 * Copyright (c) 2019 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */


#ifndef DFM2_CAM_GLFW_H
#define DFM2_CAM_GLFW_H

#include "delfem2/cam3_m4q.h"
#include <stdio.h>
#include <stdlib.h>
#include <GLFW/glfw3.h>


/**
 * @brief class for 3D navigation (rotation translation) for each GLFW window
 * @details this class should be compiled with OpenGL4.x
 */
class CNav3D_GLFW
{
public:
  CNav3D_GLFW(){
    ibutton = -1;
  }
  void Mouse(GLFWwindow *window, int button, int action, int mods){
    int win_w, win_h;
    {
      // frame buffer size might be x2 larger than the window size.
      glfwGetWindowSize(window, &win_w, &win_h);
      /*
      glfwGetFramebufferSize(window, &win_w, &win_h);
      std::cout << win_w << " " << win_h << std::endl;
      GLint viewport[4];
      ::glGetIntegerv(GL_VIEWPORT,viewport);
      win_w = viewport[2];
      win_h = viewport[3];
       */
    }
    imodifier = mods;
    double x, y;  glfwGetCursorPos (window, &x,&y);
//    std::cout << " pos: " << x << " " << y << " " << win_w << " " << win_h << std::endl;
    mouse_x = (2.0*x-win_w)/win_w;
    mouse_y = (win_h-2.0*y)/win_h;
    if( action == 0 ){
      ibutton = -1;
    }
    else if( action == 1 ){ // mouse down
      ibutton = button;
      mouse_x_down = mouse_x;
      mouse_y_down = mouse_y;
    }
  }
  void Motion(GLFWwindow *window, double x, double y){

    int win_w, win_h;
    {
//      glfwGetFramebufferSize(window, &win_w, &win_h);
      glfwGetWindowSize(window, &win_w, &win_h);
    }
    const double mov_end_x = (2.0*x-win_w)/win_w;
    const double mov_end_y = (win_h-2.0*y)/win_h;
    dx = mov_end_x - mouse_x;
    dy = mov_end_y - mouse_y;
    if( ibutton == -1 ){
    }
    else{
      if(      imodifier == 4  ){
        camera.Rot_Camera(dx, dy);
      }
      else if( imodifier == 1 ){
        camera.Pan_Camera(dx, dy);
      }
    }
    mouse_x = mov_end_x;
    mouse_y = mov_end_y;
  }

  void Key(GLFWwindow* window, int key, int scancode, int action, int mods)
  {
    if( key ==  GLFW_KEY_PAGE_UP && action == GLFW_PRESS ){ camera.Scale(1.03); }
    if( key ==  GLFW_KEY_PAGE_DOWN && action == GLFW_PRESS ){ camera.Scale(1.0/1.03); }
  }

  /**
   * @brief make 4x4 affine matrix for view transformation.
   * @details the matrix strage format is *column-major*. Applying trasformation need to multply vector from *left hand side*.
   * A 3D point is transfromed with this affine matrix and then a cube [-1,+1, -1,+1, -1,+1] is looked from -Z directoin.
   * To look from +Z direction, The transformation needs a mirror transformation in XY plane.
   * We separate mMV and  mP because of the light (light position should not be affected by the modelview transform).
   * @param[out] mMV modelview matrix (column major order)
   * @param[out] mP  projection matrix (column major order)
   * @param[in] window glfw window handler
   */
  void Mat4_MVP_OpenGL(
      float mMV[16],
      float mP[16],
      GLFWwindow* window) const
  {
    int width, height;
    glfwGetWindowSize(window, &width, &height);
    float asp = width / (float) height;
    camera.Mat4_AffineTransProjection(mP, asp, 10); // project space into cube [-1,+1,-1,+1,-1,+1] and view from -Z
    camera.Mat4_AffineTransModelView(mMV);
  }


  void MouseRay(
      float src[3],
      float dir[3],
      GLFWwindow* window) const
  {
    float mMVP_inv[16];
    {
      float mMV[16], mP[16];
      this->Mat4_MVP_OpenGL(mMV, mP, window);
      float mMVP[16];
      ::delfem2::MatMat4(mMVP, mMV,mP);
      ::delfem2::Inverse_Mat4(mMVP_inv, mMVP);
    }
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
      GLFWwindow* window) const
  {
    float mMVP_inv[16];
    {
      float mMV[16], mP[16];
      this->Mat4_MVP_OpenGL(mMV, mP, window);
      float mMVP[16];
      ::delfem2::MatMat4(mMVP, mMV,mP);
      ::delfem2::Inverse_Mat4(mMVP_inv, mMVP);
    }
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
    src1[0] = q1s[0];
    src1[1] = q1s[1];
    src1[2] = q1s[2];
    dir0[0] = q0e[0] - q0s[0];
    dir0[1] = q0e[1] - q0s[1];
    dir0[2] = q0e[2] - q0s[2];
  }

public:
  int imodifier;
  int ibutton;
  delfem2::CCamera<double> camera;
  double mouse_x, mouse_y;
  double dx, dy;
  double mouse_x_down, mouse_y_down;
};

#endif /* utility_glfw_h */
