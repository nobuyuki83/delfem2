//
// Created by Nobuyuki Umetani on 2021-11-19.
//

#include <random>

#include "gtest/gtest.h" // need to be defiend in the beginning

#include "delfem2/fem_rod2.h"
#include "delfem2/vec2.h"
#include "delfem2/sampling.h"

TEST(fem_rod2_straight, check_WdWddW) {
  namespace dfm2 = delfem2;
  std::mt19937 rndeng(std::random_device{}());
  std::uniform_real_distribution<double> dist_m1p1(-1., +1.);
  std::uniform_real_distribution<double> dist_12(+1., +2.);
  for (int nitr = 0; nitr < 100; ++nitr) {
    double aP0[3][2];
    delfem2::Fill2dArrayWithRandomValue<3,2>(aP0, dist_m1p1, rndeng);
    if (delfem2::Distance2(aP0[0], aP0[1]) < 0.1) { continue; }
    if (delfem2::Distance2(aP0[1], aP0[2]) < 0.1) { continue; }
    if (delfem2::Distance2(aP0[0], aP0[2]) < 0.1) { continue; }
    if ((dfm2::CVec2d(aP0[2]) - dfm2::CVec2d(aP0[1])).dot(
        dfm2::CVec2d(aP0[1]) - dfm2::CVec2d(aP0[0])) < 0.1) { continue; }
//
    double aL[2] = {1., 1.};
    double stiffA = dist_12(rndeng);
    double stiffB = dist_12(rndeng);
    double W0, dW0[3][2], ddW0[3][3][2][2];
    dfm2::WdWddW_Rod2(
        W0, dW0, ddW0,
        aP0, aL, stiffA, stiffA, stiffB);
    double eps = 1.0e-5;
    for (unsigned int ip = 0; ip < 3; ++ip) {
      for (unsigned int idim = 0; idim < 2; ++idim) {
        double aP1[3][2] = {
            {aP0[0][0], aP0[0][1]},
            {aP0[1][0], aP0[1][1]},
            {aP0[2][0], aP0[2][1]}};
        aP1[ip][idim] += eps;
        double W1, dW1[3][2], ddW1[3][3][2][2];
        dfm2::WdWddW_Rod2(W1, dW1, ddW1, aP1, aL, stiffA, stiffA, stiffB);
        double v0 = (W1 - W0) / eps;
        double v1 = dW0[ip][idim];
        EXPECT_NEAR(v0, v1, 2.0e-2 * (1 + fabs(v1)));
      }
    }
  }
}