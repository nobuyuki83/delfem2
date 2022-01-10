/*
 * Copyright (c) 2019 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */

#include <random>

#include "gtest/gtest.h" // need to be defiend in the beginning
//
#include "delfem2/lsmats.h" // this is necessary for merge operation in fem-related headers
#include "delfem2/geo_tri.h"
#include "delfem2/fem_discreteshell.h"
#include "delfem2/fempoisson.h"
#include "delfem2/mshmisc.h"

namespace dfm2 = delfem2;

// --------------------------------------

TEST(femem2, poisson_quad) {
  std::random_device randomDevice;
  std::mt19937 randomEng(randomDevice());
  std::uniform_real_distribution<double> dist(0.01, 3);
  double lx = dist(randomEng);
  double ly = dist(randomEng);
  double em0[4][4];
  dfm2::EMat_Poission2_QuadOrth(em0, lx, ly);
  for (unsigned int ngauss = 1; ngauss < 3; ++ngauss) {
    double em1[4][4];
    dfm2::EMat_Poisson2_QuadOrth_GaussInt(em1, lx, ly, ngauss);
    double diff = 0.0;
    for (unsigned int i = 0; i < 16; ++i) {
      double v0 = (&em0[0][0])[i];
      double v1 = (&em1[0][0])[i];
      diff += (v0 - v1) * (v0 - v1);
    }
    EXPECT_NEAR(diff, 0.0, 1.0e-10);
  }
}


TEST(pbd, Check_CdC_DiscreteShell) {
  std::mt19937 randomEng(std::random_device{}());
  std::uniform_real_distribution<double> dist_01(0, 1);
  const double eps = 1.0e-5;
  for (unsigned int itr = 0; itr < 10000; ++itr) {
    const dfm2::CVec3d p[4] = {
        dfm2::CVec3d::Random(dist_01, randomEng),
        dfm2::CVec3d::Random(dist_01, randomEng),
        dfm2::CVec3d::Random(dist_01, randomEng),
        dfm2::CVec3d::Random(dist_01, randomEng)};
    if (dfm2::Distance(p[0], p[1]) < 0.1) { continue; }
    if (dfm2::Distance(p[0], p[2]) < 0.1) { continue; }
    if (dfm2::Distance(p[0], p[3]) < 0.1) { continue; }
    if (dfm2::Distance(p[1], p[2]) < 0.1) { continue; }
    if (dfm2::Distance(p[1], p[3]) < 0.1) { continue; }
    if (dfm2::Distance(p[2], p[3]) < 0.1) { continue; }
    if (dfm2::Area_Tri3(p[0], p[2], p[3]) < 0.01) { continue; }
    if (dfm2::Area_Tri3(p[1], p[2], p[3]) < 0.01) { continue; }
    double C;
    dfm2::CVec3d dC[4];
    CdC_DiscreteShell(C, dC, p[0], p[1], p[2], p[3]);
    for (int ino = 0; ino < 4; ++ino) {
      for (int idim = 0; idim < 3; ++idim) {
        dfm2::CVec3d p1[4] = {p[0], p[1], p[2], p[3]};
        p1[ino][idim] += eps;
        double C1;
        dfm2::CVec3d dC1[4];
        CdC_DiscreteShell(C1, dC1, p1[0], p1[1], p1[2], p1[3]);
        const double val0 = (C1 - C) / eps;
        const double val1 = dC[ino][idim];
        EXPECT_NEAR(val0, val1, (1.0 + fabs(val1)) * 1.5e-2);
      }
    }
  }
}