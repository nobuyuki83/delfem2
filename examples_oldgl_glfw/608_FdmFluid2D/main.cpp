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

template<typename VAL>
VAL mymax(VAL i, VAL j) { return (i > j ? i : j); }

template<typename VAL>
VAL mymin(VAL i, VAL j) { return (i > j ? j : i); }

// Clamped Fetch
double ClampedFetch(
    double *x,
    unsigned int i,
    unsigned int j,
    unsigned int nw,
    unsigned int nh) {
  i = mymin(mymax((unsigned int) 0, i), nw - 1);
  j = mymin(mymax((unsigned int) 0, j), nh - 1);
  return x[i + j * nw];
}

// Gauss-Seidel Iteration
void SolvePoissionEquationOnGrid2ByGaussSeidelMethod(
    std::vector<double> &p,
    const std::vector<double> &d,
    unsigned int nx,
    unsigned int ny,
    unsigned int K0) {
  for (unsigned int k = 0; k < K0; k++) {
    double t = 0;
    for (unsigned int i = 0; i < nx; i++) {
      for (unsigned int j = 0; j < ny; j++) {
        const double p0 = p[nx * j + i];
        double p21 = ClampedFetch(p.data(), i + 1, j + 0, nx, ny);
        double p01 = ClampedFetch(p.data(), i - 1, j + 0, nx, ny);
        double p12 = ClampedFetch(p.data(), i + 0, j + 1, nx, ny);
        double p10 = ClampedFetch(p.data(), i + 0, j - 1, nx, ny);
        p[nx * j + i] = (p21 + p01 + p12 + p10 - d[nx * j + i]) / 4.0;
        const double p1 = p[nx * j + i];
        t += fabs(p0 - p1);
      }
    }
    if (k % 100 == 0) {
      std::cout << "itr : " << k << " " << t << std::endl;
    }
  }
}

void Divergence_StaggerdGrid2(
    std::vector<double> &divag,
    unsigned int nx_grid,
    unsigned int ny_grid,
    const std::vector<double> &velou,
    const std::vector<double> &velov,
    double h) {
  assert(velou.size() == (nx_grid + 1) * ny_grid);
  assert(velov.size() == nx_grid * (ny_grid + 1));
  for (unsigned int ig = 0; ig < nx_grid; ig++) {
    for (unsigned int jg = 0; jg < ny_grid; jg++) {
      double du = velou[(nx_grid + 1) * jg + ig + 1] - velou[(nx_grid + 1) * jg + ig];
      double dv = velov[nx_grid * (jg + 1) + ig] - velov[nx_grid * jg + ig];
      divag[ig + jg * nx_grid] = (du + dv) / h;
    }
  }
}

void EnforceBoundary(
    std::vector<double> &velou,
    std::vector<double> &velov,
    unsigned int nx_grid,
    unsigned int ny_grid) {
  assert(velou.size() == (nx_grid + 1) * ny_grid);
  for (unsigned int jg = 0; jg < ny_grid; jg++) {
    velou[(nx_grid + 1) * jg] = 0;
    velou[(nx_grid + 1) * jg + nx_grid] = 0;
  }
  assert(velov.size() == nx_grid * (ny_grid + 1));
  for (unsigned int ig = 0; ig < nx_grid; ig++) {
    velov[ig] = 0;
    velov[nx_grid * nx_grid + ig] = 0;
  }
}

void CompPressureGaussSidel(
    std::vector<double> &divergence,
    std::vector<double> &press,
    unsigned int nx_grid,
    unsigned int ny_grid,
    double h,
    double rho,
    double dt) {
  const double dtmp1 = rho / dt * h * h;
  for (unsigned int jg = 0; jg < ny_grid; jg++) {
    for (unsigned int ig = 0; ig < nx_grid; ig++) {
      divergence[nx_grid * jg + ig] *= dtmp1;
    }
  }
  SolvePoissionEquationOnGrid2ByGaussSeidelMethod(
      press,
      divergence, nx_grid, ny_grid, 1000);
}

void SubtractPressure(
    std::vector<double> &velou,
    std::vector<double> &velov,
    unsigned int ngrid,
    double h,
    double dt,
    double rho,
    const std::vector<double> &press) {
  double dtmp1 = dt / (rho * h);
  for (unsigned int ig = 0; ig < ngrid - 1; ig++) {
    for (unsigned int jg = 0; jg < ngrid; jg++) {
      double p1 = press[ngrid * jg + (ig + 1)];
      double p0 = press[ngrid * jg + (ig + 0)];
      velou[(ngrid + 1) * jg + (ig + 1)] -= dtmp1 * (p1 - p0);
    }
  }
  for (unsigned int ig = 0; ig < ngrid; ig++) {
    for (unsigned int jg = 0; jg < ngrid - 1; jg++) {
      double p1 = press[ngrid * (jg + 1) + ig];
      double p0 = press[ngrid * (jg + 0) + ig];
      velov[ngrid * (jg + 1) + ig] -= dtmp1 * (p1 - p0);
    }
  }
}

double LinearInterpolation(
    const double *d,
    unsigned int ndim,
    unsigned int idim,
    unsigned int nw,
    unsigned int nh,
    double x,
    double y) {
  x = mymax(0.0, mymin((double) nw, x));
  y = mymax(0.0, mymin((double) nh, y));
  int i = static_cast<int>(mymin(x, (double) nw - 2));
  int j = static_cast<int>(mymin(y, (double) nh - 2));
  double v0 = d[(i + j * nw) * ndim + idim];
  double v1 = d[(i + 1 + j * nw) * ndim + idim];
  double v2 = d[(i + (j + 1) * nw) * ndim + idim];
  double v3 = d[(i + 1 + (j + 1) * nw) * ndim + idim];
  double r0 = (i + 1 - x) * (j + 1 - y);
  double r1 = (x - i) * (j + 1 - y);
  double r2 = (i + 1 - x) * (y - j);
  double r3 = (x - i) * (y - j);
  return r0 * v0 + r1 * v1 + r2 * v2 + r3 * v3;
}

void CompAdvectionSemiLagrangian(
    std::vector<double> &velou,
    std::vector<double> &velov,
    std::vector<double> &velou_tmp,
    std::vector<double> &velov_tmp,
    unsigned int nx_grid,
    unsigned int ny_grid,
    double dt) {
  for (unsigned int jg = 0; jg < ny_grid + 0; jg++) {
    for (unsigned int ig = 0; ig < nx_grid + 1; ig++) {
      double v01 = ClampedFetch(velov.data(), ig - 1, jg + 0, nx_grid, ny_grid + 1);
      double v11 = ClampedFetch(velov.data(), ig + 0, jg + 0, nx_grid, ny_grid + 1);
      double v02 = ClampedFetch(velov.data(), ig - 1, jg + 1, nx_grid, ny_grid + 1);
      double v12 = ClampedFetch(velov.data(), ig + 0, jg + 1, nx_grid, ny_grid + 1);
      velou_tmp[((nx_grid + 1) * jg + ig) * 2 + 0] = velou[(nx_grid + 1) * jg + ig];
      velou_tmp[((nx_grid + 1) * jg + ig) * 2 + 1] = (v01 + v11 + v02 + v12) * 0.25;
    }
  }

  for (unsigned int jg = 0; jg < ny_grid + 1; jg++) {
    for (unsigned int ig = 0; ig < nx_grid + 0; ig++) {
      double v10 = ClampedFetch(velou.data(), ig + 0, jg - 1, nx_grid + 1, ny_grid);
      double v20 = ClampedFetch(velou.data(), ig + 1, jg - 1, nx_grid + 1, ny_grid);
      double v11 = ClampedFetch(velou.data(), ig + 0, jg + 0, nx_grid + 1, ny_grid);
      double v21 = ClampedFetch(velou.data(), ig + 1, jg + 0, nx_grid + 1, ny_grid);
      velov_tmp[(nx_grid * jg + ig) * 2 + 0] = (v10 + v20 + v11 + v21) * 0.25;
      velov_tmp[(nx_grid * jg + ig) * 2 + 1] = velov[nx_grid * jg + ig];
    }
  }

  for (unsigned int jg = 0; jg < ny_grid + 0; jg++) {
    for (unsigned int ig = 0; ig < nx_grid + 1; ig++) {
      const double *velo = velou_tmp.data() + (ig + jg * (nx_grid + 1)) * 2;
      const double p[2] = {ig - velo[0] * dt * nx_grid, jg - velo[1] * dt * ny_grid};
      const double u = LinearInterpolation(velou_tmp.data(), 2, 0, nx_grid + 1, ny_grid, p[0], p[1]);
      velou[(nx_grid + 1) * jg + ig] = u;
    }
  }

  for (unsigned int jg = 0; jg < ny_grid + 1; jg++) {
    for (unsigned int ig = 0; ig < nx_grid + 0; ig++) {
      const double *velo = velov_tmp.data() + (ig + jg * nx_grid) * 2;
      const double p[2] = {ig - velo[0] * dt * nx_grid, jg - velo[1] * dt * ny_grid};
      const double v = LinearInterpolation(velov_tmp.data(), 2, 1, nx_grid, ny_grid + 1, p[0], p[1]);
      velov[nx_grid * jg + ig] = v;
    }
  }
}

void AssignGravity(
    std::vector<double> &velou,
    std::vector<double> &velov,
    unsigned int nx_grid,
    unsigned int ny_grid,
    const double gravity[2],
    double dt) {
  assert(velou.size() == ny_grid * (nx_grid + 1));
  assert(velov.size() == (ny_grid + 1) * nx_grid);
  for (unsigned int jg = 0; jg < ny_grid + 0; jg++) {
    for (unsigned int ig = 0; ig < nx_grid + 1; ig++) {
      velou[(nx_grid + 1) * jg + ig] += gravity[0] * dt;
    }
  }
  for (unsigned int jg = 0; jg < ny_grid + 1; jg++) {
    for (unsigned int ig = 0; ig < nx_grid + 0; ig++) {
      velov[nx_grid * jg + ig] += gravity[1] * dt;
    }
  }

}

void glutMyDisplay(
    unsigned int nx_grid,
    unsigned int ny_grid,
    double h,
    const std::vector<double> &velou,
    const std::vector<double> &velov,
    const std::vector<double> &press) {
  glClear(GL_COLOR_BUFFER_BIT);
  { // quad for pressure
    ::glBegin(GL_QUADS);
    for (unsigned int jg = 0; jg < ny_grid; jg++) {
      for (unsigned int ig = 0; ig < nx_grid; ig++) {
        double p = press[nx_grid * jg + ig];
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
    for (unsigned int ig = 0; ig < nx_grid; ig++) {
      for (unsigned int jg = 0; jg < ny_grid; jg++) {
        const double p[2] = {(ig + 0.5) * h, (jg + 0.5) * h};
        double u0 = velou[(nx_grid + 1) * jg + ig + 0];
        double u1 = velou[(nx_grid + 1) * jg + ig + 1];
        const double u = (u0 + u1) * 0.5 * 3;
        double v0 = velov[nx_grid * jg + ig];
        double v1 = velov[nx_grid * (jg + 1) + ig];
        const double v = (v0 + v1) * 0.5 * 3;
        ::glVertex2d(p[0], p[1]);
        ::glVertex2d(p[0] + u, p[1] + v);
      }
    }
    ::glEnd();
  }
}

int main() {
  const unsigned int nx_grid = 32;
  const unsigned int ny_grid = 32;
  std::vector<double> velou((nx_grid + 1) * ny_grid, 0.0);  // (ngrid+1)*ngrid　
  std::vector<double> velov(nx_grid * (ny_grid + 1), 0.0);  // ngrid*(ngrid+1)
  std::vector<double> press(nx_grid * ny_grid, 0.0);  // ngrid*ngrid
  std::vector<double> divag(nx_grid * ny_grid, 0.0);  // ngrid*ngrid;
  std::vector<double> velou_tmp((nx_grid + 1) * ny_grid * 2, 0.0);  // (ngrid+1)*ngrid*2
  std::vector<double> velov_tmp(nx_grid * (ny_grid + 1) * 2, 0.0);  // ngrid*(ngrid+1)*2
  // -----------
  const double h = 1.0 / 32;
  const double dt = 0.1;
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
      for (unsigned int i = 0; i < (nx_grid + 1) * ny_grid * 2; i++) { velou[i] = dist(rndeng); }
      for (unsigned int i = 0; i < (ny_grid + 1) * nx_grid * 2; i++) { velov[i] = dist(rndeng); }
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
    CompPressureGaussSidel(
        divag, press,
        nx_grid, ny_grid, h, rho, dt);
    SubtractPressure(
        velou, velov,
        nx_grid, h, dt, rho, press);
    CompAdvectionSemiLagrangian(
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
