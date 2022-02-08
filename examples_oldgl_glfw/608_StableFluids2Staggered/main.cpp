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

#include "delfem2/glfw/viewer3.h"
#include "delfem2/glfw/util.h"

template<typename T>
class Array2 {
 public:
  Array2(int ni_, int nj_) : ni(ni_), nj(nj_) { v.resize(ni * nj); }
  Array2(int ni_, int nj_, T v_) : ni(ni_), nj(nj_) { v.resize(ni * nj, v_); }

  const T &operator()(int i, int j) const {
    assert(i >= 0 && i < ni && j >= 0 && j < nj);
    return v[i + ni * j];
  }

  T &operator()(int i, int j) {
    assert(i >= 0 && i < ni && j >= 0 && j < nj);
    return v[i + ni * j];
  }

  // Clamped Fetch
  T ClampedFetch(
      int i,
      int j) {
    auto mymax = [](int i, int j) { return i > j ? i : j; };
    auto mymin = [](int i, int j) { return i < j ? i : j; };
    i = mymax(0, i);
    j = mymax(0, j);
    i = mymin(i, ni - 1);
    j = mymin(j, nj - 1);
    assert(i >= 0 && i < ni);
    assert(j >= 0 && j < nj);
    return v[ni * j + i];
  }
 public:
  int ni, nj;
  std::vector<T> v;
};

template<int ndim>
double LinearInterpolationOnGrid2(
    const Array2<std::array<double, ndim>> &data,
    unsigned int idim,
    double x,
    double y) {
  unsigned int ni = data.ni;
  unsigned int nj = data.nj;
  auto mymin = [](double a, double b) { return a < b ? a : b; };
  auto mymax = [](double a, double b) { return a > b ? a : b; };
  double x1 = mymax(0.0, mymin((double) ni, x));
  double y1 = mymax(0.0, mymin((double) nj, y));
  auto i = static_cast<unsigned int>(mymin(x1, (double) ni - 2 - 1.0e-10));
  auto j = static_cast<unsigned int>(mymin(y1, (double) nj - 2 - 1.0e-10));
  assert(i >= 0 && i < ni - 2);
  assert(j >= 0 && j < nj - 2);
  double v00 = data(i + 0, j + 0)[idim];
  double v10 = data(i + 1, j + 0)[idim];
  double v01 = data(i + 0, j + 1)[idim];
  double v11 = data(i + 1, j + 1)[idim];
  double rx = x - i;
  double ry = y - j;
  return (1 - rx) * (1 - ry) * v00 + rx * (1 - ry) * v10 + (1 - rx) * ry * v01 + rx * ry * v11;
}

// Gauss-Seidel Iteration
void SolvePoissionEquationOnGrid2ByGaussSeidelMethod(
    Array2<double> &p,
    const Array2<double> &d,
    double scale,
    unsigned int nx,
    unsigned int ny,
    unsigned int num_iteration) {
  for (unsigned int k = 0; k < num_iteration; k++) {
    double t = 0;
    for (unsigned int j = 0; j < ny; j++) {
      for (unsigned int i = 0; i < nx; i++) {
        const double p0 = p(i, j);
        double p21 = p.ClampedFetch(i + 1, j + 0);
        double p01 = p.ClampedFetch(i - 1, j + 0);
        double p12 = p.ClampedFetch(i + 0, j + 1);
        double p10 = p.ClampedFetch(i + 0, j - 1);
        p(i, j) = (p21 + p01 + p12 + p10 - d(i, j) * scale) / 4.0;
        const double p1 = p(i, j);
        t += fabs(p0 - p1);
      }
    }
    if (k % 100 == 0) {
      std::cout << "itr : " << k << " " << t << std::endl;
    }
  }
}

void Divergence_StaggerdGrid2(
    Array2<double> &divag,
    int nx_grid,
    int ny_grid,
    const Array2<double> &velou,
    const Array2<double> &velov,
    double h) {
  assert(divag.ni == nx_grid && divag.nj == ny_grid);
  assert(velou.ni == nx_grid + 1 && velou.nj == ny_grid);
  assert(velov.ni == nx_grid && velov.nj == ny_grid + 1);
  for (int jg = 0; jg < ny_grid; jg++) {
    for (int ig = 0; ig < nx_grid; ig++) {
      double du = velou(ig + 1, jg) - velou(ig, jg);
      double dv = velov(ig, jg + 1) - velov(ig, jg);
      divag(ig, jg) = (du + dv) / h;
    }
  }
}

void EnforceBoundary(
    Array2<double> &velou,
    Array2<double> &velov,
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

void SubstructPressureGradientStaggeredGrid2(
    Array2<double> &velou,
    Array2<double> &velov,
    int nx_grid,
    int ny_grid,
    double h,
    double dt,
    double rho,
    const Array2<double> &press) {
  double dtmp1 = dt / (rho * h);
  for (int ig = 0; ig < nx_grid - 1; ig++) {
    for (int jg = 0; jg < ny_grid; jg++) {
      double p1 = press(ig + 1, jg);
      double p0 = press(ig + 0, jg);
      velou(ig + 1, jg) -= dtmp1 * (p1 - p0);
    }
  }
  for (int ig = 0; ig < nx_grid; ig++) {
    for (int jg = 0; jg < ny_grid - 1; jg++) {
      double p1 = press(ig, jg + 1);
      double p0 = press(ig, jg + 0);
      velov(ig, jg + 1) -= dtmp1 * (p1 - p0);
    }
  }
}

void CompAdvectionSemiLagrangianStaggeredGrid2(
    Array2<double> &velou,
    Array2<double> &velov,
    Array2<std::array<double, 2>> &velou_tmp,
    Array2<std::array<double, 2>> &velov_tmp,
    int nx_grid,
    int ny_grid,
    double dt) {
  assert(nx_grid > 0 && ny_grid > 0);
  for (int jg = 0; jg < ny_grid + 0; jg++) {
    for (int ig = 0; ig < nx_grid + 1; ig++) {
      double v01 = velov.ClampedFetch(ig - 1, jg + 0);
      double v11 = velov.ClampedFetch(ig + 0, jg + 0);
      double v02 = velov.ClampedFetch(ig - 1, jg + 1);
      double v12 = velov.ClampedFetch(ig + 0, jg + 1);
      velou_tmp(ig, jg) = {
          velou(ig, jg),
          (v01 + v11 + v02 + v12) * 0.25};
    }
  }

  for (int ig = 0; ig < nx_grid + 0; ig++) {
    for (int jg = 0; jg < ny_grid + 1; jg++) {
      double u10 = velou.ClampedFetch(ig + 0, jg - 1);
      double u20 = velou.ClampedFetch(ig + 1, jg - 1);
      double u11 = velou.ClampedFetch(ig + 0, jg + 0);
      double u21 = velou.ClampedFetch(ig + 1, jg + 0);
      velov_tmp(ig, jg) = {
          (u10 + u20 + u11 + u21) * 0.25,
          velov(ig, jg)};
    }
  }

  for (int jg = 0; jg < ny_grid + 0; jg++) {
    for (int ig = 0; ig < nx_grid + 1; ig++) {
      const std::array<double, 2> &velo = velou_tmp(ig, jg);
      const double p[2] = {ig - velo[0] * dt * nx_grid, jg - velo[1] * dt * ny_grid};
      const double u = LinearInterpolationOnGrid2<2>(
          velou_tmp, 0, p[0], p[1]);
      velou(ig, jg) = u;
    }
  }

  for (int jg = 0; jg < ny_grid + 1; jg++) {
    for (int ig = 0; ig < nx_grid + 0; ig++) {
      const std::array<double, 2> &velo = velov_tmp(ig, jg);
      const double p[2] = {ig - velo[0] * dt * nx_grid, jg - velo[1] * dt * ny_grid};
      const double v = LinearInterpolationOnGrid2<2>(
          velov_tmp, 1, p[0], p[1]);
      velov(ig, jg) = v;
    }
  }
}

void AssignGravity(
    Array2<double> &velou,
    Array2<double> &velov,
    int nx_grid,
    int ny_grid,
    const double gravity[2],
    double dt) {
  assert(velou.ni == nx_grid + 1 && velou.nj == ny_grid);
  assert(velov.ni == nx_grid && velov.nj == ny_grid + 1);
  for (int jg = 0; jg < ny_grid + 0; jg++) {
    for (int ig = 0; ig < nx_grid + 1; ig++) {
      velou(ig, jg) += gravity[0] * dt;
    }
  }
  for (int jg = 0; jg < ny_grid + 1; jg++) {
    for (int ig = 0; ig < nx_grid + 0; ig++) {
      velov(ig, jg) += gravity[1] * dt;
    }
  }

}

void glutMyDisplay(
    int nx_grid,
    int ny_grid,
    double h,
    const Array2<double> &velou,
    const Array2<double> &velov,
    const Array2<double> &press) {
  glClear(GL_COLOR_BUFFER_BIT);
  { // quad for pressure
    ::glBegin(GL_QUADS);
    for (int jg = 0; jg < ny_grid; jg++) {
      for (int ig = 0; ig < nx_grid; ig++) {
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
    for (int ig = 0; ig < nx_grid; ig++) {
      for (int jg = 0; jg < ny_grid; jg++) {
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
  const unsigned int nx_grid = 32;
  const unsigned int ny_grid = 40;
  Array2<double> velou(nx_grid + 1, ny_grid, 0.0);
  Array2<double> velov(nx_grid, ny_grid + 1, 0.0);
  Array2<double> press(nx_grid, ny_grid);
  Array2<double> divag(nx_grid, ny_grid);
  Array2<std::array<double, 2>> velou_tmp(nx_grid + 1, ny_grid);
  Array2<std::array<double, 2>> velov_tmp(nx_grid, ny_grid + 1);
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
      for (auto &u: velou.v) { u = dist(rndeng) * 1; }
      for (auto &v: velov.v) { v = dist(rndeng) * 1; }
    }
    // ----
    AssignGravity(
        velou, velov,
        nx_grid, ny_grid, gravity, dt);
    EnforceBoundary(
        velou, velov,
        nx_grid, ny_grid);
    Divergence_StaggerdGrid2(
        divag,
        nx_grid, ny_grid, velou, velov, h);
    SolvePoissionEquationOnGrid2ByGaussSeidelMethod(
        press,
        divag, rho / dt * h * h, nx_grid, ny_grid, 200);
    SubstructPressureGradientStaggeredGrid2(
        velou, velov,
        nx_grid, ny_grid, h, dt, rho, press);
    CompAdvectionSemiLagrangianStaggeredGrid2(
        velou, velov, velou_tmp, velov_tmp,
        nx_grid, ny_grid, dt);
    // ----
    viewer.DrawBegin_oldGL();
    glutMyDisplay(
        nx_grid, ny_grid, h,
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
