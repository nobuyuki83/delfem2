/*
 * Copyright (c) 2019 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */

#include <random>

#include "gtest/gtest.h" // need to be defiend in the beginning
//
#include "delfem2/sampling.h"
#include "delfem2/geo3_v23m34q.h"
#include "delfem2/vec3_funcs.h"
#include "delfem2/geo_tri.h"
#include "delfem2/femmips_geo3.h"
#include "delfem2/mshmisc.h"

namespace dfm2 = delfem2;

TEST(fem_mips, MIPS) {
  std::mt19937 randomEng(std::random_device{}());
  std::uniform_real_distribution<double> dist_01(0, 1);
  const double eps = 1.0e-5;
  //
  for (unsigned int itr = 0; itr < 10000; ++itr) {
    double P[3][3];
    delfem2::Fill2dArrayWithRandomValue<3,3>(P, dist_01, randomEng);
    if (dfm2::Distance3(P[0], P[1]) < 0.1) { continue; }
    if (dfm2::Distance3(P[0], P[2]) < 0.1) { continue; }
    if (dfm2::Distance3(P[1], P[2]) < 0.1) { continue; }
    if (dfm2::Area_Tri3(P[0], P[1], P[2]) < 0.01) { continue; }
    double p[3][3];
    {
      dfm2::CMat3d m;
      m.SetRotMatrix_Cartesian(0.3, 1.0, 0.5);
      for (int ino = 0; ino < 3; ++ino) {
        auto vec = m.MatVec(P[ino]);
        p[ino][0] = vec[0];
        p[ino][1] = vec[1];
        p[ino][2] = vec[2];
      }
    }
    if (dfm2::Area_Tri3(p[0], p[1], p[2]) < 0.01) { continue; }
    double E, dE[3][3], ddE[3][3][3][3];
    dfm2::WdWddW_MIPS(
        E, dE, ddE,
        p, P);
    for (int ino = 0; ino < 3; ++ino) {
      for (int idim = 0; idim < 3; ++idim) {
        double c1[3][3] = {
            {p[0][0], p[0][1], p[0][2]},
            {p[1][0], p[1][1], p[1][2]},
            {p[2][0], p[2][1], p[2][2]}};
        c1[ino][idim] += eps;
        double E1, dE1[3][3], ddE1[3][3][3][3];
        dfm2::WdWddW_MIPS(E1, dE1, ddE1,
                          c1, P);
        {
          const double val0 = (E1 - E) / eps;
          const double val1 = dE[ino][idim];
          EXPECT_NEAR(val0, val1, 1.0e-2 * (1 + fabs(val1)));
        }
        for (int jno = 0; jno < 3; ++jno) {
          for (int jdim = 0; jdim < 3; ++jdim) {
            const double val0 = (dE1[jno][jdim] - dE[jno][jdim]) / eps;
            const double val1 = ddE[jno][ino][jdim][idim];
            EXPECT_NEAR(val0, val1, 3.e-2 * (1 + fabs(val1)));
          }
        }
      }
    }
  }
}
