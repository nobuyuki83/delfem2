/*
 * Copyright (c) 2019 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */

#include "delfem2/femmitc3.h"

#include <cmath>
#include <cassert>

DFM2_INLINE void delfem2::WdWddW_PlateBendingMITC3(
    double &W,
    double dW[3][3],
    double ddW[3][3][3][3],
    const double C[3][2], // initial XY position
    const double u[3][3], // z displacement + xy axis rotation
    double thk,
    double lambda,
    double myu) {
  namespace lcl = ::delfem2::femutil;
  const double Area = lcl::TriArea2D(C[0], C[1], C[2]);
  const double Gd[3][3] = { // undeformed edge vector
      {C[1][0] - C[0][0], C[1][1] - C[0][1], 0.0},
      {C[2][0] - C[0][0], C[2][1] - C[0][1], 0.0},
      {0.0, 0.0, 0.5 * thk}};

  double Gu[3][3]; // inverse of Gd
  {
    lcl::Cross3D(Gu[0], Gd[1], Gd[2]);
    const double invtmp1 = 1.0 / lcl::Dot3D(Gu[0], Gd[0]);
    Gu[0][0] *= invtmp1;
    Gu[0][1] *= invtmp1;
    Gu[0][2] *= invtmp1;
    //
    lcl::Cross3D(Gu[1], Gd[2], Gd[0]);
    const double invtmp2 = 1.0 / lcl::Dot3D(Gu[1], Gd[1]);
    Gu[1][0] *= invtmp2;
    Gu[1][1] *= invtmp2;
    Gu[1][2] *= invtmp2;
    //
    lcl::Cross3D(Gu[2], Gd[0], Gd[1]);
    const double invtmp3 = 1.0 / lcl::Dot3D(Gu[2], Gd[2]);
    Gu[2][0] *= invtmp3;
    Gu[2][1] *= invtmp3;
    Gu[2][2] *= invtmp3;
  }

  const double GuGu2[4] = {
      lcl::Dot3D(Gu[0], Gu[0]), // rr 0
      lcl::Dot3D(Gu[1], Gu[1]), // ss 1
      lcl::Dot3D(Gu[0], Gu[1]), // sr 2
      lcl::Dot3D(Gu[2], Gu[2]), // tt 3
  };

  {
    const double CnstA[3][3] = { // {rr,ss,sr} x {rr,ss,sr}
        {lambda * GuGu2[0] * GuGu2[0] + 2 * myu * (GuGu2[0] * GuGu2[0]), // 00(0):00(0) 00(0):00(0)
         lambda * GuGu2[0] * GuGu2[1] + 2 * myu * (GuGu2[2] * GuGu2[2]), // 00(0):11(1) 01(2):01(2)
         lambda * GuGu2[0] * GuGu2[2] + 2 * myu * (GuGu2[0] * GuGu2[2])},// 00(0):01(2) 00(0):01(2)
        {lambda * GuGu2[1] * GuGu2[0] + 2 * myu * (GuGu2[2] * GuGu2[2]), // 11(1):00(0) 01(2):01(2)
         lambda * GuGu2[1] * GuGu2[1] + 2 * myu * (GuGu2[1] * GuGu2[1]), // 11(1):11(1) 11(1):11(1)
         lambda * GuGu2[1] * GuGu2[2] + 2 * myu * (GuGu2[1] * GuGu2[2])},// 11(1):01(2) 11(1):01(2)
        {lambda * GuGu2[2] * GuGu2[0] + 2 * myu * (GuGu2[0] * GuGu2[2]), // 01(2):00(0) 00(0):01(2)
         lambda * GuGu2[2] * GuGu2[1] + 2 * myu * (GuGu2[2] * GuGu2[1]), // 01(2):11(1) 11(1):01(2)
         lambda * GuGu2[2] * GuGu2[2]
             + 1 * myu * (GuGu2[0] * GuGu2[1] + GuGu2[2] * GuGu2[2])} // 01(2):01(2) 00(0):11(1) 01(2):01(2)
    };
    const double CnstB[3][3] = { // {rr,ss,sr} x {rr,ss,sr}
        {1.0 * CnstA[0][0], 1.0 * CnstA[0][1], 2.0 * CnstA[0][2]},
        {1.0 * CnstA[1][0], 1.0 * CnstA[1][1], 2.0 * CnstA[1][2]},
        {2.0 * CnstA[2][0], 2.0 * CnstA[2][1], 4.0 * CnstA[2][2]},
    };
    const double EA0t = (Gd[0][0] * (u[1][2] - u[0][2]) - Gd[0][1] * (u[1][1] - u[0][1])) * 0.5 * thk;
    const double EA1t = (Gd[1][0] * (u[2][2] - u[0][2]) - Gd[1][1] * (u[2][1] - u[0][1])) * 0.5 * thk;
    const double EA2t = (Gd[0][0] * (u[2][2] - u[0][2]) - Gd[0][1] * (u[2][1] - u[0][1])
        + Gd[1][0] * (u[1][2] - u[0][2]) - Gd[1][1] * (u[1][1] - u[0][1])) * 0.25 * thk;
    ////
    for (int iGauss = 0; iGauss < 2; ++iGauss) {
      const double t0 = (iGauss == 0) ? -1.0 / sqrt(3) : +1.0 / sqrt(3);
      const double wGauss = Area * thk / 2.0;
      const double E[3] = {t0 * EA0t, t0 * EA1t, t0 * EA2t};
      const double dE[3][3][3] = {
          {{0, +Gd[0][1] * 0.5 * thk * t0, -Gd[0][0] * 0.5 * thk * t0},
           {0, -Gd[0][1] * 0.5 * thk * t0, +Gd[0][0] * 0.5 * thk * t0},
           {0, 0, 0}},
          {{0, +Gd[1][1] * 0.5 * thk * t0, -Gd[1][0] * 0.5 * thk * t0},
           {0, 0, 0},
           {0, -Gd[1][1] * 0.5 * thk * t0, +Gd[1][0] * 0.5 * thk * t0}},
          {{0, +(Gd[0][1] + Gd[1][1]) * 0.25 * thk * t0, -(Gd[0][0] + Gd[1][0]) * 0.25 * thk * t0},
           {0, -Gd[1][1] * 0.25 * thk * t0, +Gd[1][0] * 0.25 * thk * t0},
           {0, -Gd[0][1] * 0.25 * thk * t0, +Gd[0][0] * 0.25 * thk * t0}}};
      ////
      const double SB[3] = {
          CnstB[0][0] * E[0] + CnstB[0][1] * E[1] + CnstB[0][2] * E[2],
          CnstB[1][0] * E[0] + CnstB[1][1] * E[1] + CnstB[1][2] * E[2],
          CnstB[2][0] * E[0] + CnstB[2][1] * E[1] + CnstB[2][2] * E[2]};
      W += wGauss * 0.5 * (E[0] * SB[0] + E[1] * SB[1] + E[2] * SB[2]);
      for (int ino = 0; ino < 3; ++ino) {
        for (int idof = 0; idof < 3; ++idof) {
          dW[ino][idof] += wGauss * (SB[0] * dE[0][ino][idof] + SB[1] * dE[1][ino][idof] + SB[2] * dE[2][ino][idof]);
        }
      }
      for (int ino = 0; ino < 3; ++ino) {
        for (int jno = 0; jno < 3; ++jno) {
          for (int idof = 0; idof < 3; ++idof) {
            for (int jdof = 0; jdof < 3; ++jdof) {
              const double dtmp =
                  +dE[0][ino][idof] * CnstB[0][0] * dE[0][jno][jdof]
                      + dE[0][ino][idof] * CnstB[0][1] * dE[1][jno][jdof]
                      + dE[0][ino][idof] * CnstB[0][2] * dE[2][jno][jdof]
                      + dE[1][ino][idof] * CnstB[1][0] * dE[0][jno][jdof]
                      + dE[1][ino][idof] * CnstB[1][1] * dE[1][jno][jdof]
                      + dE[1][ino][idof] * CnstB[1][2] * dE[2][jno][jdof]
                      + dE[2][ino][idof] * CnstB[2][0] * dE[0][jno][jdof]
                      + dE[2][ino][idof] * CnstB[2][1] * dE[1][jno][jdof]
                      + dE[2][ino][idof] * CnstB[2][2] * dE[2][jno][jdof];
              ddW[ino][jno][idof][jdof] += wGauss * dtmp;
            }
          }
        }
      }
    }
  }
  {
    const double CnstA[2][2] = { // {rt,st} x {rt,st}
        {myu * GuGu2[0] * GuGu2[3],  // rt*rt -> rr(0):tt(3)
         myu * GuGu2[2] * GuGu2[3]}, // st*rt -> sr(2):tt(3)
        {myu * GuGu2[2] * GuGu2[3],  // rt*st -> rs(2):tt(3)
         myu * GuGu2[1] * GuGu2[3]}  // st*st -> ss(1):tt(3)
    };
    const double CnstB[2][2] = {
        {4.0 * CnstA[0][0], 2.0 * CnstA[0][1]},
        {2.0 * CnstA[1][0], 4.0 * CnstA[1][1]}};
    const double Ert_01 =
        0.5 * thk * (u[1][0] - u[0][0] + 0.5 * Gd[0][0] * (u[0][2] + u[1][2]) - 0.5 * Gd[0][1] * (u[0][1] + u[1][1]));
    const double Ert_12 =
        0.5 * thk * (u[1][0] - u[0][0] + 0.5 * Gd[0][0] * (u[1][2] + u[2][2]) - 0.5 * Gd[0][1] * (u[1][1] + u[2][1]));
    const double Est_12 =
        0.5 * thk * (u[2][0] - u[0][0] + 0.5 * Gd[1][0] * (u[1][2] + u[2][2]) - 0.5 * Gd[1][1] * (u[1][1] + u[2][1]));
    const double Est_20 =
        0.5 * thk * (u[2][0] - u[0][0] + 0.5 * Gd[1][0] * (u[2][2] + u[0][2]) - 0.5 * Gd[1][1] * (u[2][1] + u[0][1]));
    const double dErt_01[3][3] = {
        {-0.5 * thk, -0.25 * thk * Gd[0][1], +0.25 * thk * Gd[0][0]},
        {+0.5 * thk, -0.25 * thk * Gd[0][1], +0.25 * thk * Gd[0][0]},
        {0, 0, 0}};
    const double dEst_20[3][3] = {
        {-0.5 * thk, -0.25 * thk * Gd[1][1], +0.25 * thk * Gd[1][0]},
        {0, 0, 0},
        {+0.5 * thk, -0.25 * thk * Gd[1][1], +0.25 * thk * Gd[1][0]}};
    const double dEC[3][3] = {
        {0, +0.25 * thk * Gd[0][1] - 0.25 * thk * Gd[1][1], -0.25 * thk * Gd[0][0] + 0.25 * thk * Gd[1][0]},
        {0, +0.25 * thk * Gd[1][1], -0.25 * thk * Gd[1][0]},
        {0, -0.25 * thk * Gd[0][1], +0.25 * thk * Gd[0][0]}};
    const double CE = (Ert_12 - Ert_01) - (Est_12 - Est_20);
    const double pGauss[3][2] = {{0.5, 0.0}, {0.5, 0.5}, {0.0, 0.5}};
    for (int iGauss = 0; iGauss < 3; ++iGauss) {
      const double r = pGauss[iGauss][0];
      const double s = pGauss[iGauss][1];
      const double wGauss = Area * thk / 3.0;
      const double E[2] = {Ert_01 + CE * s, Est_20 - CE * r};
      double dE[2][3][3];
      for (int ino = 0; ino < 3; ++ino) {
        for (int idof = 0; idof < 3; ++idof) {
          dE[0][ino][idof] = dErt_01[ino][idof] + dEC[ino][idof] * s;
          dE[1][ino][idof] = dEst_20[ino][idof] - dEC[ino][idof] * r;
        }
      }
      const double SB[2] = {
          CnstB[0][0] * E[0] + CnstB[0][1] * E[1],
          CnstB[1][0] * E[0] + CnstB[1][1] * E[1]};
      W += wGauss * 0.5 * (SB[0] * E[0] + SB[1] * E[1]);
      for (int ino = 0; ino < 3; ++ino) {
        for (int idof = 0; idof < 3; ++idof) {
          dW[ino][idof] += wGauss * (SB[0] * dE[0][ino][idof] + SB[1] * dE[1][ino][idof]);
        }
      }
      for (int ino = 0; ino < 3; ++ino) {
        for (int jno = 0; jno < 3; ++jno) {
          for (int idof = 0; idof < 3; ++idof) {
            for (int jdof = 0; jdof < 3; ++jdof) {
              const double dtmp =
                  +dE[0][ino][idof] * CnstB[0][0] * dE[0][jno][jdof]
                      + dE[0][ino][idof] * CnstB[0][1] * dE[1][jno][jdof]
                      + dE[1][ino][idof] * CnstB[1][0] * dE[0][jno][jdof]
                      + dE[1][ino][idof] * CnstB[1][1] * dE[1][jno][jdof];
              ddW[ino][jno][idof][jdof] += wGauss * dtmp;
            }
          }
        }
      }
    }
  }
}

void delfem2::MassLumped_ShellPlateBendingMITC3(
    double *aM,
    double rho, double thick,
    const double *aXY, unsigned int nXY,
    const unsigned int *aTri, unsigned int nTri) {
  const unsigned int nDoF = nXY * 3;
  for (unsigned int i = 0; i < nDoF; ++i) { aM[i] = 0.0; }
  for (unsigned int it = 0; it < nTri; ++it) {
    const unsigned int i0 = aTri[it * 3 + 0];
    assert(i0 < nXY);
    const unsigned int i1 = aTri[it * 3 + 1];
    assert(i1 < nXY);
    const unsigned int i2 = aTri[it * 3 + 2];
    assert(i2 < nXY);
    const double *p0 = aXY + i0 * 2;
    const double *p1 = aXY + i1 * 2;
    const double *p2 = aXY + i2 * 2;
    const double a012 = delfem2::femutil::TriArea2D(p0, p1, p2);
    double m0 = a012 / 3.0 * rho * thick;
    double m1 = a012 / 3.0 * rho * thick * thick * thick / 12.0;
    double m2 = m1;
    aM[i0 * 3 + 0] += m0;
    aM[i1 * 3 + 0] += m0;
    aM[i2 * 3 + 0] += m0;
    aM[i0 * 3 + 1] += m1;
    aM[i1 * 3 + 1] += m1;
    aM[i2 * 3 + 1] += m1;
    aM[i0 * 3 + 2] += m2;
    aM[i1 * 3 + 2] += m2;
    aM[i2 * 3 + 2] += m2;
  }
}



