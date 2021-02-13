#ifndef DFM2_CVBUNDLEADJUST_H
#define DFM2_CVBUNDLEADJUST_H

#include "delfem2/geo3_v23m34q.h"
#include <vector>

namespace delfem2 {

void ReprojectionEnergyPoint(
    double F[2],
    double dF[2][7],
    //
    const double p0[2],
    const double K0inv[9],
    float Z0_point,
    const double t01[3],
    const double K1Rt01[9],
    const double p1[2]) {
  namespace dfm2 = delfem2;
  F[0] = 0.0;
  F[1] = 0.0;
  for (int i = 0; i < 2 * 7; ++i) { (&dF[0][0])[i] = 0.0; }
  //
  double x0 = p0[0];
  double y0 = p0[1];
  const double P0[3] = {
      K0inv[0] * x0 + K0inv[1] * y0 + K0inv[2],
      K0inv[3] * x0 + K0inv[4] * y0 + K0inv[5],
      K0inv[6] * x0 + K0inv[7] * y0 + K0inv[8]};
  const double P0a[3] = {
      P0[0] * (Z0_point / P0[2]),
      P0[1] * (Z0_point / P0[2]),
      Z0_point};
  const double P1[3] = {
      P0a[0] - t01[0],
      P0a[1] - t01[1],
      P0a[2] - t01[2]};
  double p1a[3];
  dfm2::MatVec3(p1a, K1Rt01, P1);
  const double p1b[2] = {p1a[0] / p1a[2], p1a[1] / p1a[2]};
  F[0] = p1b[0] - p1[0];
  F[1] = p1b[1] - p1[1];
  { // dEdZ
    const double dP0a[3] = {P0[0] / P0[2], P0[1] / P0[2], 1.0};
    double dp1a[3];
    dfm2::MatVec3(dp1a, K1Rt01, dP0a);
    const double dp1b[2] = {
        dp1a[0] / p1a[2] - p1a[0] / (p1a[2] * p1a[2]) * dp1a[2],
        dp1a[1] / p1a[2] - p1a[1] / (p1a[2] * p1a[2]) * dp1a[2]};
    dF[0][0] = dp1b[0];
    dF[1][0] = dp1b[1];
  }
  const double tmpXY[2][3] = {
      {-1 / p1a[2], 0,           +p1a[0] / (p1a[2] * p1a[2])},
      {0,           -1 / p1a[2], +p1a[1] / (p1a[2] * p1a[2])}};
  dfm2::MatTVec3(&dF[0][1], K1Rt01, tmpXY[0]);
  dfm2::MatTVec3(&dF[1][1], K1Rt01, tmpXY[1]);
  { // dEdr
    double tmp3[3];
    dfm2::MatTVec3(tmp3, K1Rt01, tmpXY[0]);
    dfm2::Cross3(&dF[0][4], P1, tmp3);
    dfm2::MatTVec3(tmp3, K1Rt01, tmpXY[1]);
    dfm2::Cross3(&dF[1][4], P1, tmp3);
  }
}

double ReprojectionEnergy(
    std::vector<double> &F,
    std::vector<double> &dF,
    //
    const std::vector<double> &aPnt0,
    const double K0inv[9],
    const std::vector<float> &aZ0,
    const double R01[9],
    const double t01[3],
    const double K1[9],
    const std::vector<double> &aPnt1) {
  namespace dfm2 = delfem2;
  const unsigned int np = aPnt0.size() / 2;
  assert(aPnt1.size() == 2 * np);
  double K1Rt01[9];
  dfm2::MatMatT3(K1Rt01, K1, R01);
  F.resize(np * 2);
  dF.resize(np * 2 * 7);
  double eng = 0.0;
  for (unsigned int ipc = 0; ipc < np; ++ipc) {
    double f0i[2], df0i[2][7];
    ReprojectionEnergyPoint(
        f0i, df0i,
        aPnt0.data() + ipc * 2,
        K0inv,
        aZ0[ipc],
        t01,
        K1Rt01,
        aPnt1.data() + ipc * 2);
    eng += f0i[0] * f0i[0] + f0i[1] * f0i[1];
    F[ipc * 2 + 0] = f0i[0];
    F[ipc * 2 + 1] = f0i[1];
    for (int i = 0; i < 7; ++i) {
      dF[ipc * 14 + i + 0] = df0i[0][i];
      dF[ipc * 14 + i + 7] = df0i[1][i];
    }
  }
  return eng;
}

void BundleAdjust_ImgPair_Stochastic(
    std::vector<float> &aZ0,
    double R01[9],
    double t01[3],
    //
    const std::vector<double> &aPnt0,
    const double K0inv[9],
    const double K1[9],
    const std::vector<double> &aPnt1,
    std::mt19937 &rdeng) {
  namespace dfm2 = delfem2;
  std::uniform_real_distribution<float> dist_m1p1(-1, 1);
  std::vector<double> F, dF;
  {
    const float E_A = ReprojectionEnergy(
        F, dF,
        aPnt0, K0inv, aZ0, R01, t01, K1, aPnt1);
    dfm2::CVec3d t01b(
        t01[0] + dist_m1p1(rdeng) * 0.01,
        t01[1] + dist_m1p1(rdeng) * 0.01,
        t01[2] + dist_m1p1(rdeng) * 0.01);
    dfm2::CMat3d dR;
    dR.SetRotMatrix_BryantAngle(
        dist_m1p1(rdeng) * 0.01,
        dist_m1p1(rdeng) * 0.01,
        dist_m1p1(rdeng) * 0.01);
    const dfm2::CMat3d R01b = dfm2::CMat3d(R01) * dR;
    const float E_B = ReprojectionEnergy(
        F, dF,
        aPnt0, K0inv, aZ0, R01b.mat, t01b.p, K1, aPnt1);
    if (E_B < E_A) { // fails. put it back..
      std::cout << E_B << std::endl;
      t01b.CopyTo(t01);
      R01b.CopyTo(R01);
    }
  }
  {
    double K1Rt01[9];
    dfm2::MatMatT3(K1Rt01, K1, R01);
    const int npc = aPnt0.size() / 2;
    for (int ipc = 0; ipc < npc; ++ipc) {
      const float z1a = aZ0[ipc];
      double FA[2], dFA[2][7];
      ReprojectionEnergyPoint(
          FA, dFA,
          aPnt0.data() + ipc * 2,
          K0inv, z1a, t01, K1Rt01,
          aPnt1.data() + ipc * 2);
      float z1b = z1a + dist_m1p1(rdeng) * 0.003f;
      if (z1b > 0) { z1b = 0.; }
      double FB[2], dFB[2][7];
      ReprojectionEnergyPoint(
          FB, dFB,
          aPnt0.data() + ipc * 2,
          K0inv, z1b, t01, K1Rt01,
          aPnt1.data() + ipc * 2);
      if (FB[0] * FB[0] + FB[1] * FB[1] < FA[0] * FA[0] + FA[1] * FA[1]) {
        aZ0[ipc] = z1b;
      }
    }
  }
}

void BundleAdjust_CheckGrad(
    const std::vector<double> &aPnt0,
    const double *K0inv,
    const std::vector<float> &aZ0,
    double *R01,
    double *t01,
    const double *K1,
    const std::vector<double> &aPnt1)
 {
  const unsigned int np = aPnt0.size() / 2;
  assert(aPnt1.size() == np * 2);
  namespace dfm2 = delfem2;
  double K1Rt01[9];
  dfm2::MatMatT3(K1Rt01, K1, R01);
  for (unsigned int ip = 0; ip < np; ++ip) {
    double eps = 1.0e-5;
    //
    double FA[2], dFA[2][7];
    ReprojectionEnergyPoint(
        FA, dFA,
        aPnt0.data() + ip * 2,
        K0inv, aZ0[ip], t01, K1Rt01,
        aPnt1.data() + ip * 2);
    {
      double FB[2], dFB[2][7];
      ReprojectionEnergyPoint(
          FB, dFB,
          aPnt0.data() + ip * 2,
          K0inv, aZ0[ip] + eps, t01, K1Rt01,
          aPnt1.data() + ip * 2);
      {
        double d0 = (FB[0] - FA[0]) / eps;
        double r0 = fabs(dFA[0][0] - d0) / (1.e-3 + fabs(dFA[0][0]));
        double d1 = (FB[1] - FA[1]) / eps;
        double r1 = fabs(dFA[1][0] - d1) / (1.e-3 + fabs(dFA[1][0]));
        std::cout << ip << " --> " << r0 << " " << r1 << std::endl;
      }
    }
    for (int idimt = 0; idimt < 3; ++idimt) {
      double t01_B[3] = {t01[0], t01[1], t01[2]};
      t01_B[idimt] += eps;
      double FB[2], dFB[2][7];
      ReprojectionEnergyPoint(
          FB, dFB,
          aPnt0.data() + ip * 2,
          K0inv, aZ0[ip], t01_B, K1Rt01,
          aPnt1.data() + ip * 2);
      {
        double d0 = (FB[0] - FA[0]) / eps;
        double r0 = fabs(dFA[0][idimt + 1] - d0) / (1.e-3 + fabs(dFA[0][idimt + 1]));
        double d1 = (FB[1] - FA[1]) / eps;
        double r1 = fabs(dFA[1][idimt + 1] - d1) / (1.e-3 + fabs(dFA[1][idimt + 1]));
//        std::cout << ip << " --> " << d0 << " " << d1 << "  " << dFA[0][idimt+1] << " " << dFA[1][idimt+1] << std::endl;
        std::cout << ip << " --> " << r0 << " " << r1 << std::endl;
      }
    }
    for (int idimt = 0; idimt < 3; ++idimt) {
      double dw[3] = {0, 0, 0};
      dw[idimt] += eps;
      double dR[9];
      dfm2::Mat3_Rotation_Cartesian(dR, dw);
      double K1Rt01_B[9];
      dfm2::MatMatT3(K1Rt01_B, K1Rt01, dR);
      double FB[2], dFB[2][7];
      ReprojectionEnergyPoint(
          FB, dFB,
          aPnt0.data() + ip * 2,
          K0inv, aZ0[ip], t01, K1Rt01_B,
          aPnt1.data() + ip * 2);
      {
        double d0 = (FB[0] - FA[0]) / eps;
        double r0 = fabs(dFA[0][idimt + 4] - d0) / (1.e-3 + fabs(dFA[0][idimt + 4]));
        double d1 = (FB[1] - FA[1]) / eps;
        double r1 = fabs(dFA[1][idimt + 4] - d1) / (1.e-3 + fabs(dFA[1][idimt + 4]));
//        std::cout << ip << " --> " << d0 << " " << d1 << "  " << dFA[0][idimt+4] << " " << dFA[1][idimt+4] << std::endl;
        std::cout << ip << " --> " << r0 << " " << r1 << std::endl;
      }
    }
  }
}

}

#endif