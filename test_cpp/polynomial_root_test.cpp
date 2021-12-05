//
// Created by Nobuyuki Umetani on 2021/12/05.
//

#include "gtest/gtest.h" // need to be defined in the beginning

#include <random>
#include <vector>

#include "delfem2/polynomial_root.h"

TEST(polynomial_root, test0) {
  std::mt19937 rndeng(std::random_device{}());
  std::uniform_real_distribution<double> dist_m1p1(-1, +1);
  for (unsigned int itr = 0; itr < 1000; ++itr) {
    const double a = dist_m1p1(rndeng);
    const double b = dist_m1p1(rndeng);
    const double c = dist_m1p1(rndeng);
    auto roots0 = delfem2::RootsInRange_QuadraticPolynomial(0., 1., a, b, c);
    int nroot1 = delfem2::NumberOfRootsInRange_QuadraticPolynomial(0., 1., a, b, c);
    EXPECT_EQ(roots0.size(), nroot1);
  }
}