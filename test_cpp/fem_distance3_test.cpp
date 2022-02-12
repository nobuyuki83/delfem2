/*
 * Copyright (c) 2019 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */

#include <random>

#include "gtest/gtest.h" // need to be defiend in the beginning
//
#include "delfem2/fem_distance3.h"
#include "delfem2/mshmisc.h"
#include "delfem2/sampling.h"
#include "delfem2/vec2.h"

namespace dfm2 = delfem2;

TEST(fem_distance3, distancetri2d3d) {
  std::random_device randomDevice;
  std::mt19937 randomEng(randomDevice());
  std::uniform_real_distribution<double> dist_01(0, 1);
  //
  for (unsigned int itr = 0; itr < 10000; ++itr) {
    double P[3][2];  // undeformed triangle vertex positions
    delfem2::Fill2dArrayWithRandomValue<3,2>(P, dist_01, randomEng);
    if (dfm2::Distance2(P[0], P[1]) < 0.1) { continue; }
    if (dfm2::Distance2(P[0], P[2]) < 0.1) { continue; }
    if (dfm2::Distance2(P[1], P[2]) < 0.1) { continue; }
    double p[3][3]; //  deformed triangle vertex positions)
    delfem2::Fill2dArrayWithRandomValue<3,3>(p, dist_01, randomEng);
    // --------------
    double C[3], dCdp[3][9];
    dfm2::CdC_LengthTriEdges23(C, dCdp, P, p);
    for (int ine = 0; ine < 3; ++ine) {
      for (int idim = 0; idim < 3; ++idim) {
        double eps = 1.0e-6;
        double p1[3][3];
        for (int i = 0; i < 9; ++i) { (&p1[0][0])[i] = (&p[0][0])[i]; }
        p1[ine][idim] += eps;
        double C1[3], dCdp1[3][9];
        dfm2::CdC_LengthTriEdges23(C1, dCdp1, P, p1);
        for (int jdim = 0; jdim < 3; ++jdim) {
          const double val0 = (C1[jdim] - C[jdim]) / eps;
          const double val1 = dCdp[jdim][ine * 3 + idim];
          EXPECT_NEAR(val0, val1, 1.0e-2 * (fabs(val1) + 1.0));
        }
      }
    }
  }
}

TEST(fem_distance3, WdWddW_SquareLengthLineseg3D) {
  std::mt19937 rndeng(std::random_device{}());
  std::uniform_real_distribution<double> dist_m1p1(-1, +1);
  std::uniform_real_distribution<double> dist_12(+1, +2);
  std::uniform_real_distribution<double> dist_01(+0, +1);
  double eps = 1.0e-5;
  //
  for (int itr = 0; itr < 10000; ++itr) {
    const double stiff_stretch = dist_12(rndeng);
    const dfm2::CVec3d P[2] = {
        dfm2::RandomVec<3>(dist_01, rndeng),
        dfm2::RandomVec<3>(dist_01, rndeng) };
    if ((P[0] - P[1]).norm() < 0.1) { continue; }
    dfm2::CVec3d dW_dP[2];
    dfm2::CMat3d ddW_ddP[2][2];
    const double L0 = 1.0;
    double W = WdWddW_SquareLengthLineseg3D(
        dW_dP, ddW_ddP,
        stiff_stretch, P, L0);
    // -----
    const dfm2::CVec3d dP[2] = {
        dfm2::RandomVec<3>(dist_01, rndeng, eps),
        dfm2::RandomVec<3>(dist_01, rndeng, eps) };
    const dfm2::CVec3d p[2] = {P[0] + dP[0], P[1] + dP[1]};
    double w;
    dfm2::CVec3d dw_dP[2];
    {
      dfm2::CMat3d ddw_ddP[2][2];
      w = WdWddW_SquareLengthLineseg3D(
          dw_dP, ddw_ddP,
          stiff_stretch, p, L0);
    }
    {
      const double val0 = (w - W) / eps;
      const double val1 = (+dW_dP[0].dot(dP[0]) + dW_dP[1].dot(dP[1])) / eps;
      EXPECT_NEAR(val0, val1, 1.0e-3 * (1 + fabs(val1)));
    }
    {
      const dfm2::CVec3d val0 = (dw_dP[0] - dW_dP[0]) / eps;
      const dfm2::CVec3d val1 = (+ddW_ddP[0][0] * dP[0] + ddW_ddP[0][1] * dP[1]) / eps;
      EXPECT_LT((val0 - val1).norm(), 1.0e-3 * (1 + val1.norm()));
    }
    {
      const dfm2::CVec3d val0 = (dw_dP[1] - dW_dP[1]) / eps;
      const dfm2::CVec3d val1 = (+ddW_ddP[1][0] * dP[0] + ddW_ddP[1][1] * dP[1]) / eps;
      EXPECT_LT((val0 - val1).norm(), 1.0e-3 * (1 + val1.norm()));
    }
  }
}


