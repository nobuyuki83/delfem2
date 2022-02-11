/*
 * Copyright (c) 2019 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */

#include <iostream>
#include <vector>
#include <cmath>
#include <random>
#if defined(_WIN32) // windows
#  define NOMINMAX   // to remove min,max macro
#  include <windows.h>  // this should come before glfw3.h
#endif
#define GL_SILENCE_DEPRECATION
#include <GLFW/glfw3.h>

#include "delfem2/fdm_array2.h"
#include "delfem2/colormap.h"
#include "delfem2/glfw/viewer3.h"
#include "delfem2/glfw/util.h"

void PressureUpdate_StaggerdGrid2(
    FdmArray2<double> &press,
    int ni_grid,
    int nj_grid,
    const FdmArray2<double> &velou,
    const FdmArray2<double> &velov,
    double h,
    double rho,
    double c,
    double dt) {
  assert(press.ni == ni_grid && press.nj == nj_grid);
  assert(velou.ni == ni_grid + 1 && velou.nj == nj_grid);
  assert(velov.ni == ni_grid && velov.nj == nj_grid + 1);
  double t0 = - dt * rho * c * c;
  for (int jg = 0; jg < nj_grid; jg++) {
    for (int ig = 0; ig < ni_grid; ig++) {
      double du = velou(ig + 1, jg) - velou(ig, jg);
      double dv = velov(ig, jg + 1) - velov(ig, jg);
      press(ig, jg) += t0 * (du + dv) / h;
    }
  }
}

void VelocityUpdate_StaggeredGrid2(
    FdmArray2<double> &velou,
    FdmArray2<double> &velov,
    int ni_grid,
    int nj_grid,
    const FdmArray2<double> &press,
    double dt,
    double rho,
    double h) {
  assert(velou.ni == ni_grid + 1 && velou.nj == nj_grid);
  assert(velov.ni == ni_grid && velov.nj == nj_grid + 1);
  assert(press.ni == ni_grid && press.nj == nj_grid);
  double scale = dt / rho * h * h;
  for (int ig = 0; ig < ni_grid - 1; ig++) {
    for (int jg = 0; jg < nj_grid; jg++) {
      double p1 = press(ig + 1, jg);
      double p0 = press(ig + 0, jg);
      velou(ig + 1, jg) -= scale * (p1 - p0);
    }
  }
  for (int ig = 0; ig < ni_grid; ig++) {
    for (int jg = 0; jg < nj_grid - 1; jg++) {
      double p1 = press(ig, jg + 1);
      double p0 = press(ig, jg + 0);
      velov(ig, jg + 1) -= scale * (p1 - p0);
    }
  }
}

void glutMyDisplay(
    int ni_grid,
    int nj_grid,
    double h,
    const FdmArray2<double> &velou,
    const FdmArray2<double> &velov,
    const FdmArray2<double> &press,
    double scale,
    const std::tuple<float,float,float> colormap[],
    unsigned int ncolormap) {
  assert(velou.ni == ni_grid + 1 && velou.nj == nj_grid);
  assert(velov.ni == ni_grid && velov.nj == nj_grid + 1);
  glClear(GL_COLOR_BUFFER_BIT);
  double pmin = -0.5;
  double pmax = +0.5;
  { // quad for pressure
    ::glBegin(GL_QUADS);
    for (int jg = 0; jg < nj_grid; jg++) {
      for (int ig = 0; ig < ni_grid; ig++) {
        double p = press(ig, jg);
        int ic = static_cast<int>(ncolormap*(p-pmin)/(pmax-pmin));
        if( ic < 0 ){ ic = 0; }
        if( ic >= (int)ncolormap ){ ic = (int)ncolormap-1; }
        ::glColor3f(
            std::get<0>(colormap[ic]),
            std::get<1>(colormap[ic]),
            std::get<2>(colormap[ic]));
        ::glVertex2d((ig + 0) * h, (jg + 0) * h);
        ::glVertex2d((ig + 1) * h, (jg + 0) * h);
        ::glVertex2d((ig + 1) * h, (jg + 1) * h);
        ::glVertex2d((ig + 0) * h, (jg + 1) * h);
      }
    }
    ::glEnd();
  }

  { // draw velocity
    ::glColor3d(0, 0, 0);
    ::glBegin(GL_LINES);
    for (int ig = 0; ig < ni_grid; ig++) {
      for (int jg = 0; jg < nj_grid; jg++) {
        const double p[2] = {(ig + 0.5) * h, (jg + 0.5) * h};
        double u0 = velou(ig + 0, jg);
        double u1 = velou(ig + 1, jg);
        const double u = (u0 + u1) * 0.5;
        double v0 = velov(ig, jg + 0);
        double v1 = velov(ig, jg + 1);
        const double v = (v0 + v1) * 0.5;
        ::glVertex2d(p[0], p[1]);
        ::glVertex2d(
            p[0] + u * scale,
            p[1] + v * scale);
      }
    }
    ::glEnd();
  }
}

int main() {
  const unsigned int ni_grid = 32;
  const unsigned int nj_grid = 40;
  FdmArray2<double> velou(ni_grid + 1, nj_grid, 0.0);
  FdmArray2<double> velov(ni_grid, nj_grid + 1, 0.0);
  FdmArray2<double> press(ni_grid, nj_grid, 0.00);
  // -----------
  const double h = 1.0 / 32;
  const double dt = 0.1;
  const double rho = 1.0;
  const double c = 10.0;

  constexpr auto& colormap = delfem2::colormap_plasma<float>;
  unsigned int ncolormap = sizeof(colormap) / sizeof(colormap[0]);

  delfem2::glfw::CViewer3 viewer;
  delfem2::glfw::InitGLOld();
  viewer.OpenWindow();

  for(unsigned int iframe=0;;++iframe){
    // pressure source
    press(ni_grid/2, nj_grid/2) = sin(iframe*dt);
    //
    PressureUpdate_StaggerdGrid2(
        press,
        ni_grid, nj_grid, velou, velov,
        h, rho, c, dt);
    VelocityUpdate_StaggeredGrid2(
        velou, velov,
        ni_grid, nj_grid, press,
        dt, rho, h);
    // ----
    viewer.DrawBegin_oldGL();
    glutMyDisplay(
        ni_grid, nj_grid, h,
        velou, velov, press,
        100,
        colormap, ncolormap);
    viewer.SwapBuffers();
    glfwPollEvents();
    if (glfwWindowShouldClose(viewer.window)) goto EXIT;
  }
  EXIT:
  glfwDestroyWindow(viewer.window);
  glfwTerminate();
  exit(EXIT_SUCCESS);
}
