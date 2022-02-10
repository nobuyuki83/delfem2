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

void EnforceNonSlipBoundary_CellCenteredGrid2(
    FdmArray2<std::array<double,2>> &velo,
    int nx_grid,
    int ny_grid) {
  assert(velo.ni == nx_grid && velo.nj == ny_grid);
  for (int jg = 0; jg < ny_grid; jg++) {
    velo(0, jg)[0] = 0;
    velo(nx_grid-1, jg)[0] = 0;
  }
  for (int ig = 0; ig < nx_grid; ig++) {
    velo(ig, 0)[1] = 0;
    velo(ig, ny_grid-1)[1] = 0;
  }
}

void AssignGravity_CellCnteredGrid2(
    FdmArray2<std::array<double,2>> &velo,
    int ni_grid,
    int nj_grid,
    const double gravity[2],
    double dt) {
  assert(velo.ni == ni_grid && velo.nj == nj_grid);
  for (int jg = 0; jg < nj_grid; jg++) {
    for (int ig = 0; ig < ni_grid; ig++) {
      velo(ig, jg)[0] += gravity[0] * dt;
      velo(ig, jg)[1] += gravity[1] * dt;
    }
  }
}

void glutMyDisplay(
    int ni_grid,
    int nj_grid,
    double h,
    const FdmArray2<std::array<double,2>> &velo,
    const FdmArray2<double> &press) {
  assert(velo.ni == ni_grid && velo.nj == nj_grid);
  assert(press.ni == ni_grid && press.nj == nj_grid);

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
        const std::array<double,2>& v = velo(ig, jg);
        ::glVertex2d(p[0], p[1]);
        ::glVertex2d(p[0] + v[0], p[1] + v[1]);
      }
    }
    ::glEnd();
  }
}

int main() {
  const unsigned int ni_grid = 32;
  const unsigned int nj_grid = 40;
  FdmArray2<std::array<double,2>> velo(ni_grid, nj_grid, {0.0,0.0});
  FdmArray2<double> press(ni_grid, nj_grid);
  FdmArray2<double> divag(ni_grid, nj_grid);
  FdmArray2<std::array<double, 2>> velo_tmp(ni_grid, nj_grid);
  // -----------
  const double h = 1.0 / 32;
  const double dt = 0.05;
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
      for (auto &u: velo.v) { u = {dist(rndeng), dist(rndeng)}; }
    }
    // ----
    AssignGravity_CellCnteredGrid2(
        velo,
        ni_grid, nj_grid, gravity, dt);
    EnforceNonSlipBoundary_CellCenteredGrid2(
        velo,
        ni_grid, nj_grid);
    Divergence_CellCenteredGrid2(
        divag,
        ni_grid, nj_grid, velo, h);
    SolvePoissionEquationOnGrid2ByGaussSeidelMethod(
        press,
        divag, rho / dt * h * h, ni_grid, nj_grid, 200);
    SubstructPressureGradient_CellCenteredGrid2(
        velo,
        ni_grid, nj_grid, dt / (rho * h), press);
    AdvectionSemiLagrangian_CellCenteredGrid2(
        velo, velo_tmp,
        ni_grid, nj_grid, dt, h);
    // ----
    viewer.DrawBegin_oldGL();
    glutMyDisplay(
        ni_grid, nj_grid, h,
        velo, press);
    viewer.SwapBuffers();
    glfwPollEvents();
    if (glfwWindowShouldClose(viewer.window)) goto EXIT;
  }
  EXIT:
  glfwDestroyWindow(viewer.window);
  glfwTerminate();
  exit(EXIT_SUCCESS);
}
