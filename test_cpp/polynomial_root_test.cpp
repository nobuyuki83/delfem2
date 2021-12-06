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
    double strum[3][3];
    delfem2::StrumSequence_QuadraticPolynomial(strum, a,b,c);
    int nroot1 = delfem2::StrumNumber<3>(0,strum) - delfem2::StrumNumber<3>(1,strum);
    EXPECT_EQ(roots0.size(), nroot1);
  }
}



TEST(polynomial_root, cubic_3roots) {
  std::mt19937 rndeng(std::random_device{}());
  std::uniform_real_distribution<double> dist_m1p1(-1, +1);
  for (unsigned int itr = 0; itr < 1000; ++itr) {
    const double a = dist_m1p1(rndeng);
    const double x0 = dist_m1p1(rndeng) * 2;
    const double x1 = dist_m1p1(rndeng) * 2;
    const double x2 = dist_m1p1(rndeng) * 2;
    const double b = -a * (x0+x1+x2);
    const double c = a * (x0*x1+x1*x2+x2*x0);
    const double d = -a * x0 * x1 * x2;
    double strum[4][4];
    delfem2::StrumSequence_CubicPolynomial(strum, a, b, c, d);
    int nroot1 = delfem2::StrumNumber<4>(0,strum) - delfem2::StrumNumber<4>(1,strum);
    {
      int nroot0 = 0;
      if (x0 > 0 && x0 < 1) { nroot0++; }
      if (x1 > 0 && x1 < 1) { nroot0++; }
      if (x2 > 0 && x2 < 1) { nroot0++; }
      EXPECT_EQ(nroot0, nroot1);
    }
    {
      std::vector<std::pair<double, double>> intvls = delfem2::RootInterval_StrumSequence<4>(0, 1, strum);
      EXPECT_EQ(intvls.size(), nroot1);
      for(auto intvl : intvls){
        bool b0 = (intvl.first-x0)*(intvl.second-x0) < 0;
        bool b1 = (intvl.first-x1)*(intvl.second-x1) < 0;
        bool b2 = (intvl.first-x2)*(intvl.second-x2) < 0;
        EXPECT_TRUE( b0 || b1 || b2 );
        int nr = delfem2::StrumNumber<4>(intvl.first,strum) - delfem2::StrumNumber<4>(intvl.second,strum);
        EXPECT_EQ(nr,1);
      }
    }
  }
}

TEST(polynomial_root, quad_2roots) {
  std::mt19937 rndeng(std::random_device{}());
  std::uniform_real_distribution<double> dist_m1p1(-1, +1);
  for (unsigned int itr = 0; itr < 1000; ++itr) {
    const double a = dist_m1p1(rndeng);
    const double x0 = dist_m1p1(rndeng) * 2;
    const double x1 = dist_m1p1(rndeng) * 2;
    const double b = -a * (x0+x1);
    const double c = a * x0 * x1;
    double strum[3][3];
    delfem2::StrumSequence_QuadraticPolynomial(strum, a,b,c);
    int nroot1 = delfem2::StrumNumber<3>(0,strum) - delfem2::StrumNumber<3>(1,strum);
    {
      int nroot0 = 0;
      if (x0 > 0 && x0 < 1) { nroot0++; }
      if (x1 > 0 && x1 < 1) { nroot0++; }
      EXPECT_EQ(nroot0, nroot1);
    }
    {
      std::vector<std::pair<double, double>> intvls = delfem2::RootInterval_StrumSequence<3>(0, 1, strum);
      EXPECT_EQ(intvls.size(), nroot1);
      for(auto intvl : intvls){
        bool b0 = (intvl.first-x0)*(intvl.second-x0) < 0;
        bool b1 = (intvl.first-x1)*(intvl.second-x1) < 0;
        EXPECT_TRUE( b0 || b1 );
        int nr = delfem2::StrumNumber<3>(intvl.first,strum) - delfem2::StrumNumber<3>(intvl.second,strum);
        EXPECT_EQ(nr,1);
      }
    }
  }
}