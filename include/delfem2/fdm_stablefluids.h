//
// Created by Nobuyuki Umetani on 2022/02/10.
//

#ifndef DFM2_FDM_STABLEFLUIDS_H_
#define DFM2_FDM_STABLEFLUIDS_H_

#include <cassert>
#include <array>

#include "delfem2/fdm_array2.h"

// Gauss-Seidel Iteration
void SolvePoissionEquationOnGrid2ByGaussSeidelMethod(
    FdmArray2<double> &p,
    const FdmArray2<double> &d,
    double scale,
    int nx,
    int ny,
    unsigned int num_iteration) {
  assert(p.ni == nx && p.nj == ny);
  assert(d.ni == nx && d.nj == ny);
  for (unsigned int k = 0; k < num_iteration; k++) {
    double t = 0;
    for (int j = 0; j < ny; j++) {
      for (int i = 0; i < nx; i++) {
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
    FdmArray2<double> &divag,
    int ni_grid,
    int nj_grid,
    const FdmArray2<double> &velou,
    const FdmArray2<double> &velov,
    double h) {
  assert(divag.ni == ni_grid && divag.nj == nj_grid);
  assert(velou.ni == ni_grid + 1 && velou.nj == nj_grid);
  assert(velov.ni == ni_grid && velov.nj == nj_grid + 1);
  for (int jg = 0; jg < nj_grid; jg++) {
    for (int ig = 0; ig < ni_grid; ig++) {
      double du = velou(ig + 1, jg) - velou(ig, jg);
      double dv = velov(ig, jg + 1) - velov(ig, jg);
      divag(ig, jg) = (du + dv) / h;
    }
  }
}

void Divergence_CellCenteredGrid2(
    FdmArray2<double> &divag,
    int ni_grid,
    int nj_grid,
    const FdmArray2<std::array<double, 2>> &velo,
    double h) {
  assert(divag.ni == ni_grid && divag.nj == nj_grid);
  assert(velo.ni == ni_grid && velo.nj == nj_grid);
  for (int jg = 1; jg < nj_grid-1; jg++) {
    for (int ig = 1; ig < ni_grid-1; ig++) {
      double right_u = velo(ig + 1, jg)[0];
      double left_u = velo(ig - 1, jg)[0];
      double up_v = velo(ig, jg + 1)[1];
      double down_v = velo(ig, jg - 1)[1];
      divag(ig, jg) = (left_u - right_u + down_v - up_v) / h;
    }
  }
}

void SubstructPressureGradient_StaggeredGrid2(
    FdmArray2<double> &velou,
    FdmArray2<double> &velov,
    int ni_grid,
    int nj_grid,
    double scale,
    const FdmArray2<double> &press) {
  assert(velou.ni == ni_grid + 1 && velou.nj == nj_grid);
  assert(velov.ni == ni_grid && velov.nj == nj_grid + 1);
  assert(press.ni == ni_grid && press.nj == nj_grid);
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

void SubstructPressureGradient_CellCenteredGrid2(
    FdmArray2<std::array<double, 2>> &velo,
    int ni_grid,
    int nj_grid,
    double scale,
    const FdmArray2<double> &press) {
  assert(velo.ni == ni_grid && velo.nj == nj_grid);
  assert(press.ni == ni_grid && press.nj == nj_grid);
  for (int ig = 1; ig < ni_grid-1; ig++) {
    for (int jg = 1; jg < nj_grid-1; jg++) {
      double right = press(ig + 1, jg);
      double left = press(ig - 1, jg);
      double up = press(ig, jg + 1);
      double down = press(ig, jg - 1);
      velo(ig, jg)[0] -= scale * (left - right) * 0.5;
      velo(ig, jg)[1] -= scale * (down - up) * 0.5;
    }
  }
}

void AdvectionSemiLagrangian_StaggeredGrid2(
    FdmArray2<double> &velou,
    FdmArray2<double> &velov,
    FdmArray2<std::array<double, 2>> &velou_tmp,
    FdmArray2<std::array<double, 2>> &velov_tmp,
    int ni_grid,
    int nj_grid,
    double dt,
    double h) {
  assert(velou.ni == ni_grid + 1 && velou.nj == nj_grid);
  assert(velov.ni == ni_grid && velov.nj == nj_grid + 1);
  assert(ni_grid > 0 && nj_grid > 0);
  for (int jg = 0; jg < nj_grid + 0; jg++) {
    for (int ig = 0; ig < ni_grid + 1; ig++) {
      double v01 = velov.ClampedFetch(ig - 1, jg + 0);
      double v11 = velov.ClampedFetch(ig + 0, jg + 0);
      double v02 = velov.ClampedFetch(ig - 1, jg + 1);
      double v12 = velov.ClampedFetch(ig + 0, jg + 1);
      velou_tmp(ig, jg) = {
          velou(ig, jg),
          (v01 + v11 + v02 + v12) * 0.25};
    }
  }

  for (int ig = 0; ig < ni_grid + 0; ig++) {
    for (int jg = 0; jg < nj_grid + 1; jg++) {
      double u10 = velou.ClampedFetch(ig + 0, jg - 1);
      double u20 = velou.ClampedFetch(ig + 1, jg - 1);
      double u11 = velou.ClampedFetch(ig + 0, jg + 0);
      double u21 = velou.ClampedFetch(ig + 1, jg + 0);
      velov_tmp(ig, jg) = {
          (u10 + u20 + u11 + u21) * 0.25,
          velov(ig, jg)};
    }
  }

  for (int jg = 0; jg < nj_grid + 0; jg++) {
    for (int ig = 0; ig < ni_grid + 1; ig++) {
      const std::array<double, 2> &velo = velou_tmp(ig, jg);
      const double p[2] = {ig - velo[0] * dt / h, jg - velo[1] * dt / h};
      const double u = LinearInterpolationOnGrid2<2>(velou_tmp, 0, p[0], p[1]);
      velou(ig, jg) = u;
    }
  }

  for (int jg = 0; jg < nj_grid + 1; jg++) {
    for (int ig = 0; ig < ni_grid + 0; ig++) {
      const std::array<double, 2> &velo = velov_tmp(ig, jg);
      const double p[2] = {ig - velo[0] * dt / h, jg - velo[1] * dt / h};
      const double v = LinearInterpolationOnGrid2<2>(velov_tmp, 1, p[0], p[1]);
      velov(ig, jg) = v;
    }
  }
}

void AdvectionSemiLagrangian_CellCenteredGrid2(
    FdmArray2<std::array<double, 2>> &velo,
    FdmArray2<std::array<double, 2>> &velo_tmp,
    int ni_grid,
    int nj_grid,
    double dt,
    double h) {
  assert(velo.ni == ni_grid && velo.nj == nj_grid);
  assert(velo_tmp.ni == ni_grid && velo_tmp.nj == nj_grid);
  assert(ni_grid > 0 && nj_grid > 0);
  velo_tmp.v = velo.v;
  for (int jg = 0; jg < nj_grid; jg++) {
    for (int ig = 0; ig < ni_grid; ig++) {
      const std::array<double, 2> &vold = velo_tmp(ig, jg);
      const double p[2] = {ig - vold[0] * dt / h, jg - vold[1] * dt / h};
      const double u = LinearInterpolationOnGrid2<2>(velo_tmp, 0, p[0], p[1]);
      const double v = LinearInterpolationOnGrid2<2>(velo_tmp, 1, p[0], p[1]);
      velo(ig, jg) = {u, v};
    }
  }
}

#endif //DFM2_FDM_STABLEFLUIDS_H_
