//
// Created by Nobuyuki Umetani on 2022/02/01.
//


#ifndef DFM2_FEM_INVERTIBLEFEM_H
#define DFM2_FEM_INVERTIBLEFEM_H

#include <functional>

#include "delfem2/geo3_v23m34q.h"
#include "delfem2/mat3.h"
#include "delfem2/geo_tet.h"
#include "delfem2/svd3.h"

namespace delfem2 {

// dPdF[k][l][i][j] = dP_kl / dF_ij
double DiffPiolaKirchhoff1st(
    delfem2::CMat3d &P0,
    double dPdF[3][3][3][3],
    const delfem2::CMat3d &U0,
    const delfem2::CMat3d &S0,
    const delfem2::CMat3d &V0,
    const double diff[9][3][3],
    const std::function<double(double[3], double[3][3], double, double, double)> &neohook) {
  namespace dfm2 = delfem2;
  double dWdl0[3], ddWddl0[3][3];
  const double W0 = neohook(dWdl0, ddWddl0, S0(0, 0), S0(1, 1), S0(2, 2));
  const dfm2::CMat3d T0(dWdl0[0], dWdl0[1], dWdl0[2]);
  P0 = U0 * (T0 * V0.transpose());
  for (unsigned i = 0; i < 3; ++i) {
    for (unsigned j = 0; j < 3; ++j) {
      const dfm2::CMat3d dU = -U0 * dfm2::CMat3d::Spin(
          std::array<double, 3>{diff[0][i][j], diff[1][i][j], diff[2][i][j]}.data());
      const dfm2::CMat3d dV = -V0 * dfm2::CMat3d::Spin(
          std::array<double, 3>{diff[6][i][j], diff[7][i][j], diff[8][i][j]}.data());
      const dfm2::CMat3d dS(diff[3][i][j], diff[4][i][j], diff[5][i][j]);
      const double t0 = ddWddl0[0][0] * dS(0, 0) + ddWddl0[0][1] * dS(1, 1) + ddWddl0[0][2] * dS(2, 2);
      const double t1 = ddWddl0[1][0] * dS(0, 0) + ddWddl0[1][1] * dS(1, 1) + ddWddl0[1][2] * dS(2, 2);
      const double t2 = ddWddl0[2][0] * dS(0, 0) + ddWddl0[2][1] * dS(1, 1) + ddWddl0[2][2] * dS(2, 2);
      const dfm2::CMat3d dT(t0, t1, t2);
      const dfm2::CMat3d dP = dU * T0 * V0.transpose()
          + U0 * dT * V0.transpose()
          + U0 * (T0 * dV.transpose());
      for (unsigned k = 0; k < 3; ++k) {
        for (unsigned l = 0; l < 3; ++l) {
          dPdF[k][l][i][j] = dP(k, l);
        }
      }
    }
  }
  return W0;
}

double WdWddW_InvertibleFEM(
    double dW[4][3],
    double ddW[4][4][3][3],
    const double Pos[4][3],
    const double pos[4][3],
    const std::function<double(double[3], double[3][3], double, double, double)> &neohook) {
  namespace dfm2 = delfem2;
  const dfm2::CMat3d F0 = DeformationGradientOfTet(
      Pos[0], Pos[1], Pos[2], Pos[3],
      pos[0], pos[1], pos[2], pos[3]);

  double dFdu[4][3];
  delfem2::DiffDeformationGradientOfTet(
      dFdu,
      Pos[0], Pos[1], Pos[2], Pos[3]);

  const auto[U0, S0, V0] = dfm2::Svd3<dfm2::CMat3d>(F0, 30);

  double diff[9][3][3];
  Svd3Differential(
      diff,
      U0, S0, V0);

  dfm2::CMat3d P0;
  double dPdF[3][3][3][3];
  double W0 = DiffPiolaKirchhoff1st(
      P0, dPdF,
      U0, S0, V0, diff, neohook);

  for (unsigned int ino = 0; ino < 4; ++ino) {
    for (unsigned int i = 0; i < 3; ++i) {
      dW[ino][i] = P0(i, 0) * dFdu[ino][0] + P0(i, 1) * dFdu[ino][1] + P0(i, 2) * dFdu[ino][2];
    }
    for (unsigned int jno = 0; jno < 4; ++jno) {
      for (unsigned int i = 0; i < 3; ++i) {
        for (unsigned int j = 0; j < 3; ++j) {
          ddW[ino][jno][i][j] = dPdF[i][0][j][0] * dFdu[ino][0] * dFdu[jno][0] +
              dPdF[i][0][j][1] * dFdu[ino][0] * dFdu[jno][1] +
              dPdF[i][0][j][2] * dFdu[ino][0] * dFdu[jno][2] +
              dPdF[i][1][j][0] * dFdu[ino][1] * dFdu[jno][0] +
              dPdF[i][1][j][1] * dFdu[ino][1] * dFdu[jno][1] +
              dPdF[i][1][j][2] * dFdu[ino][1] * dFdu[jno][2] +
              dPdF[i][2][j][0] * dFdu[ino][2] * dFdu[jno][0] +
              dPdF[i][2][j][1] * dFdu[ino][2] * dFdu[jno][1] +
              dPdF[i][2][j][2] * dFdu[ino][2] * dFdu[jno][2];
        }
      }
    }
  }

  return W0;
}

}


#endif /* DFM2_FEM_INVERTIBLEFEM_H */