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
#include "delfem2/fdm_stablefluids.h"
#include "delfem2/glfw/viewer3.h"
#include "delfem2/glfw/util.h"

void EnforceNonSlipBoundary_StaggeredGrid2(
    FdmArray2<double> &velou,
    FdmArray2<double> &velov,
    int nx_grid,
    int ny_grid) {
  assert(velou.ni == nx_grid + 1 && velou.nj == ny_grid);
  assert(velov.ni == nx_grid && velov.nj == ny_grid + 1);
  for (int jg = 0; jg < ny_grid; jg++) {
    velou(0, jg) = 0;
    velou(nx_grid, jg) = 0;
  }
  for (int ig = 0; ig < nx_grid; ig++) {
    velov(ig, 0) = 0;
    velov(ig, ny_grid) = 0;
  }
}

void AssignGravity_StaggeredGrid2(
    FdmArray2<double> &velou,
    FdmArray2<double> &velov,
    int ni_grid,
    int nj_grid,
    const double gravity[2],
    double dt) {
  assert(velou.ni == ni_grid + 1 && velou.nj == nj_grid);
  assert(velov.ni == ni_grid && velov.nj == nj_grid + 1);
  for (int jg = 0; jg < nj_grid + 0; jg++) {
    for (int ig = 0; ig < ni_grid + 1; ig++) {
      velou(ig, jg) += gravity[0] * dt;
    }
  }
  for (int jg = 0; jg < nj_grid + 1; jg++) {
    for (int ig = 0; ig < ni_grid + 0; ig++) {
      velov(ig, jg) += gravity[1] * dt;
    }
  }
}

void glutMyDisplay(
    int ni_grid,
    int nj_grid,
    double h,
    const FdmArray2<double> &velou,
    const FdmArray2<double> &velov,
    const FdmArray2<double> &press) {
  assert(velou.ni == ni_grid + 1 && velou.nj == nj_grid);
  assert(velov.ni == ni_grid && velov.nj == nj_grid + 1);
  glClear(GL_COLOR_BUFFER_BIT);
  { // quad for pressure
    ::glBegin(GL_QUADS);
    for (int jg = 0; jg < nj_grid; jg++) {
      for (int ig = 0; ig < ni_grid; ig++) {
        double p = press(ig, jg);
        glColor4d(p > 0, 0.0, p < 0, 0.8);
        ::glVertex2d((ig + 0) * h, (jg + 0) * h);
        ::glVertex2d((ig + 1) * h, (jg + 0) * h);
        ::glVertex2d((ig + 1) * h, (jg + 1) * h);
        ::glVertex2d((ig + 0) * h, (jg + 1) * h);
      }
    }
    ::glEnd();
  }

  { // draw velocity
    ::glColor3d(0, 1, 1);
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
        ::glVertex2d(p[0] + u, p[1] + v);
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
  FdmArray2<double> press(ni_grid, nj_grid);
  FdmArray2<double> divag(ni_grid, nj_grid);
  FdmArray2<std::array<double, 2>> velou_tmp(ni_grid + 1, nj_grid);
  FdmArray2<std::array<double, 2>> velov_tmp(ni_grid, nj_grid + 1);
  // -----------
  const double h = 1.0 / 32;
  const double dt = 0.02;
  const double rho = 1.0;
  const double gravity[2] = {0, -0.01};
  std::mt19937 rndeng(std::random_device{}());
  std::uniform_real_distribution<> dist(-1.0, 1.0);

  delfem2::glfw::CViewer3 viewer;
  delfem2::glfw::InitGLOld();
  viewer.OpenWindow();

  int iframe = -1;
  while (true) {
    iframe = (iframe + 1) % 100;
    if (iframe == 0) {
      for (auto &u: velou.v) { u = dist(rndeng); }
      for (auto &v: velov.v) { v = dist(rndeng); }
    }
    // ----
    AssignGravity_StaggeredGrid2(
        velou, velov,
        ni_grid, nj_grid, gravity, dt);
    EnforceNonSlipBoundary_StaggeredGrid2(
        velou, velov,
        ni_grid, nj_grid);
    Divergence_StaggerdGrid2(
        divag,
        ni_grid, nj_grid, velou, velov, h);
    SolvePoissionEquationOnGrid2ByGaussSeidelMethod(
        press,
        divag, rho / dt * h * h, ni_grid, nj_grid, 200);
    SubstructPressureGradient_StaggeredGrid2(
        velou, velov,
        ni_grid, nj_grid, dt / (rho * h), press);
    AdvectionSemiLagrangian_StaggeredGrid2(
        velou, velov, velou_tmp, velov_tmp,
        ni_grid, nj_grid, dt, h);
    // ----
    viewer.DrawBegin_oldGL();
    glutMyDisplay(
        ni_grid, nj_grid, h,
        velou, velov, press);
    viewer.SwapBuffers();
    glfwPollEvents();
    if (glfwWindowShouldClose(viewer.window)) goto EXIT;
  }
  EXIT:
  glfwDestroyWindow(viewer.window);
  glfwTerminate();
  exit(EXIT_SUCCESS);
}
