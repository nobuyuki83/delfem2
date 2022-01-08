/*
 * Copyright (c) 2019 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */


#include <cassert>
#include <vector>
#include <random>
#if defined(_WIN32) // windows
#  define NOMINMAX   // to remove min,max macro
#  include <windows.h>  // this should come before glfw3.h
#endif
#define GL_SILENCE_DEPRECATION
#include <GLFW/glfw3.h>

#include "delfem2/vec3.h"
#include "delfem2/vec2.h"
#include "delfem2/geo_curve_cubic.h"
#include "delfem2/opengl/old/mshuni.h"
#include "delfem2/glfw/viewer3.h"
#include "delfem2/glfw/util.h"

namespace dfm2 = delfem2;

void Draw(
    const std::vector< dfm2::CVec3d > &cps,
    const std::vector< dfm2::CVec3d > &sample){
  ::glPointSize(5);
  ::glColor3d(0,0,0);
  ::glBegin(GL_POINTS);
  for(const auto& cp : cps){
    ::glVertex3dv(cp.data());
  }
  ::glEnd();
  //
  ::glLineWidth(1);
  ::glColor3d(0,0,0);
  ::glBegin(GL_LINE_STRIP);
  for(const auto& cp : cps){
    ::glVertex3dv(cp.data());
  }
  ::glEnd();
  //
  ::glLineWidth(1);
  ::glColor3d(1,0,0);
  ::glBegin(GL_LINE_STRIP);
  for(const auto& p : sample){
    ::glVertex3dv(p.data());
  }
  ::glEnd();
}

int main()
{
  std::vector< dfm2::CVec3d > cps = {
      {-0.4, 0, 0},
      {-0.06, 0, 0},
      {-0.02, 0.45, 0.5},
      {-0.34, 0.45, 0.5},
      {-0.31, 1, 0.2},
//      {0, 1, 0.0},
      {0.31, 1, 0.2},
      {0.34, 0.45, 0.5},
      {0.02, 0.45, 0.5},
      {0.06, 0, 0},
      {0.4, 0, 0}};
  std::vector< dfm2::CVec3d > sample;
  delfem2::Polyline_CubicBSplineCurve(
      sample,
      1001, cps);
  //
  delfem2::glfw::CViewer3 viewer(1);
  //
  delfem2::glfw::InitGLOld();
  viewer.OpenWindow();
  //
  while (!glfwWindowShouldClose(viewer.window)) {
    //
    viewer.DrawBegin_oldGL();
    Draw(cps,sample);
    ::glPushMatrix();
    ::glTranslated(0,0.4,0);
    Draw(cps,sample);
    ::glPopMatrix();
    glfwSwapBuffers(viewer.window);
    glfwPollEvents();
  }
  glfwDestroyWindow(viewer.window);
  glfwTerminate();
  exit(EXIT_SUCCESS);
}
