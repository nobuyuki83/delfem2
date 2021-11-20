//
// Created by Nobuyuki Umetani on 2021/11/20.
//


#include <random>

#include "gtest/gtest.h" // need to be defiend in the beginning

#include "delfem2/vec2.h"
#include "delfem2/fem_stvk.h"

bool RandomTri2(double P[3][2]) {
  namespace dfm2 = delfem2;
  std::mt19937 randomEng(std::random_device{}());
  std::uniform_real_distribution<double> dist_m5p5(-1, 1);
  for (int ip = 0; ip < 3; ++ip) {
    for (int idim = 0; idim < 2; ++idim) {
      P[ip][idim] = dist_m5p5(randomEng);
    }
  }
  if (dfm2::Distance2(P[0], P[1]) < 0.1) { return false; }
  if (dfm2::Distance2(P[1], P[2]) < 0.1) { return false; }
  if (dfm2::Distance2(P[0], P[2]) < 0.1) { return false; }
  double a0 = dfm2::Area_Tri2(P[0], P[1], P[2]);
  if (fabs(a0) < 0.2) return false;
  return true;
}

TEST(fem_stvk, Check_CdC_TriStrain) {
  namespace dfm2 = delfem2;
  std::mt19937 randomEng(std::random_device{}());
  std::uniform_real_distribution<double> dist_m5p5(-5, 5);
// -----
  for (int itr = 0; itr < 200; ++itr) {
    double P[3][2];
    if (!RandomTri2(P)) { continue; }
    const double p[3][3] = {
      {dist_m5p5(randomEng), dist_m5p5(randomEng), dist_m5p5(randomEng)},
      {dist_m5p5(randomEng), dist_m5p5(randomEng), dist_m5p5(randomEng)},
      {dist_m5p5(randomEng), dist_m5p5(randomEng), dist_m5p5(randomEng)}};
    const double eps = 1.0e-5;
    double C[3], dCdp[3][9];
    dfm2::CdC_StVK(C, dCdp, P, p);
    for (int ine = 0; ine < 3; ++ine) {
      for (int idim = 0; idim < 3; ++idim) {
        double p1[3][3];
        for (int i = 0; i < 9; ++i) { (&p1[0][0])[i] = (&p[0][0])[i]; }
        p1[ine][idim] += eps;
        double C1[3], dCdp1[3][9];
        dfm2::CdC_StVK(C1, dCdp1, P, p1);
        EXPECT_NEAR((C1[0] - C[0]) / eps, dCdp[0][ine * 3 + idim], 1.0e-2);
        EXPECT_NEAR((C1[1] - C[1]) / eps, dCdp[1][ine * 3 + idim], 1.0e-2);
        EXPECT_NEAR((C1[2] - C[2]) / eps, dCdp[2][ine * 3 + idim], 1.0e-2);
      }
    }
  }
}

TEST(fem_stvk, Check_CdC_EnergyStVK) {
  namespace dfm2 = delfem2;
  std::mt19937 randomEng(std::random_device{}());
  std::uniform_real_distribution<double> dist_01(0, 1);
  std::uniform_real_distribution<double> dist_12(1, 2);
  //
  for (int itr = 0; itr < 1000; ++itr) {
    double P[3][2];
    if (!RandomTri2(P)) { continue; }
    const double p[3][3] = { //  deformed triangle vertex positions)
      {dist_01(randomEng), dist_01(randomEng), dist_01(randomEng)},
      {dist_01(randomEng), dist_01(randomEng), dist_01(randomEng)},
      {dist_01(randomEng), dist_01(randomEng), dist_01(randomEng)}};
    const double lambda = dist_12(randomEng);
    const double myu = dist_12(randomEng);
    // ---------------------
    double C, dCdp[9];
    dfm2::CdC_EnergyStVK(C, dCdp, P, p, lambda, myu);
    for (int ine = 0; ine < 3; ++ine) {
      for (int idim = 0; idim < 3; ++idim) {
        double eps = 1.0e-8;
        double p1[3][3];
        for (int i = 0; i < 9; ++i) { (&p1[0][0])[i] = (&p[0][0])[i]; }
        p1[ine][idim] += eps;
        double C1, dCdp1[9];
        dfm2::CdC_EnergyStVK(C1, dCdp1, P, p1, lambda, myu);
        EXPECT_NEAR(
          (C1 - C) / eps,
          dCdp[ine * 3 + idim],
          1.0e-3 * (1 + fabs(dCdp[ine * 3 + idim])));
      }
    }
  }
}
