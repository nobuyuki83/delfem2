/*
 * Copyright (c) 2019 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */


#ifndef DFM2_CAM_GLFW_H
#define DFM2_CAM_GLFW_H

#include <stdio.h>
#include <stdlib.h>

#include <GLFW/glfw3.h>

#include "delfem2/camera.h"


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
   * @param[out] mMV modelview matrix (column major order)
   * @param[out] mP  projection matrix (column major order)
   * @param[in] window glfw window handler
   */
  void Matrix_MVP(float mMV[16],
                  float mP[16],
                  GLFWwindow* window) const
  {
    float asp;
    int width, height;
    glfwGetWindowSize(window, &width, &height);
    asp = width / (float) height;
    camera.Affine4f_Projection(mP, asp, 10);
    camera.Affine4f_ModelView(mMV);
  }
  void PosMouse2D(float& x, float& y,
                  GLFWwindow* window) const
  {
    float mMV[16], mP[16]; this->Matrix_MVP(mMV, mP, window);
    const float sp0[3] = {(float)mouse_x, (float)mouse_y,0.0};
    float src_pick[3];
    delfem2::screenUnProjection(src_pick,
                       sp0, mMV,mP);
    x = src_pick[0];
    y = src_pick[1];
  }
  void MouseRay(float src[3], float dir[3],
                GLFWwindow* window) const
  {
    float mMV[16], mP[16]; this->Matrix_MVP(mMV, mP, window);
    {
      const float sp0[3] = {(float)mouse_x, (float)mouse_y,0.0};
      delfem2::screenUnProjection(src,
                         sp0, mMV,mP);
    }
    {
      const float dir0[3] = {0.0, 0.0, +1.0};
      delfem2::screenUnProjectionDirection(dir,
                                           dir0, mMV,mP);
    }
  }
  void RayMouseMove(float src0[3], float src1[3], float dir[3],
                    GLFWwindow* window) const
  {
    float mMV[16], mP[16]; this->Matrix_MVP(mMV, mP, window);
    {
      const float sp0[3] = {(float)(mouse_x-dx), (float)(mouse_y-dy),0.0};
      delfem2::screenUnProjection(src0,
                         sp0, mMV,mP);
    }
    {
      const float sp1[3] = {(float)mouse_x, (float)mouse_y,0.0};
      delfem2::screenUnProjection(src1,
                         sp1, mMV,mP);
    }
    {
      const float dir0[3] = {0.0, 0.0, -1.0};
      delfem2::screenUnProjectionDirection(dir,
                                  dir0, mMV,mP);
    }
  }
  void PosMove2D(float& x0, float& y0,
                 float& x1, float& y1,
                 GLFWwindow* window) const
  {
    float mMV[16], mP[16]; this->Matrix_MVP(mMV, mP, window);
    {
      const float sp0[3] = {(float)(mouse_x-dx), (float)(mouse_y-dy),0.0};
      float src0[3];
      delfem2::screenUnProjection(src0,
                         sp0, mMV,mP);
      x0 = src0[0]; y0 = src0[1];
    }
    {
      const float sp1[3] = {(float)mouse_x, (float)mouse_y,0.0};
      float src1[3];
      delfem2::screenUnProjection(src1,
                                  sp1, mMV,mP);
      x1 = src1[0]; y1 = src1[1];
    }
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
