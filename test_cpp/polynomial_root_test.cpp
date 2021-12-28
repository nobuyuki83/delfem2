//
// Created by Nobuyuki Umetani on 2021/12/05.
//

#include "gtest/gtest.h" // need to be defined in the beginning

#include <random>
#include <vector>
#include <array>

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
    delfem2::StrumSequenceOfPolynomial<3>(strum, std::array<double,3>{c,b,a}.data());
    int nroot1 = delfem2::StrumNumber<3>(0, strum) - delfem2::StrumNumber<3>(1, strum);
    EXPECT_EQ(roots0.size(), nroot1);
  }
}

TEST(polynomial_root, quad_2roots) {
  std::mt19937 rndeng(std::random_device{}());
  std::uniform_real_distribution<double> dist_m1p1(-1, +1);
  for (unsigned int itr = 0; itr < 1000; ++itr) {
    const double a = dist_m1p1(rndeng);
    const double x0 = dist_m1p1(rndeng) * 2;
    const double x1 = dist_m1p1(rndeng) * 2;
    const double b = -a * (x0 + x1);
    const double c = a * x0 * x1;
    double strum[3][3];
    delfem2::StrumSequenceOfPolynomial(strum, std::array<double,3>{c,b,a}.data());
    int nroot1 = delfem2::StrumNumber<3>(0, strum) - delfem2::StrumNumber<3>(1, strum);
    {
      int nroot0 = 0;
      if (x0 > 0 && x0 < 1) { nroot0++; }
      if (x1 > 0 && x1 < 1) { nroot0++; }
      EXPECT_EQ(nroot0, nroot1);
    }
    {
      std::vector<std::pair<double, double>> intvls = delfem2::RootInterval_StrumSequence<3>(0, 1, strum);
      EXPECT_EQ(intvls.size(), nroot1);
      for (auto intvl: intvls) {
        const bool b0 = (intvl.first - x0) * (intvl.second - x0) < 0;
        const bool b1 = (intvl.first - x1) * (intvl.second - x1) < 0;
        EXPECT_TRUE(b0 || b1);
        int nr = delfem2::StrumNumber<3>(intvl.first, strum) - delfem2::StrumNumber<3>(intvl.second, strum);
        EXPECT_EQ(nr, 1);
        auto coeff = std::array<double, 3>{c, b, a};
        const double v0 = delfem2::Eval_Polynomial(intvl.first, coeff.data(), 3);
        const double v1 = delfem2::Eval_Polynomial(intvl.second, coeff.data(), 3);
        EXPECT_LT(v0 * v1, 0.);
        const double y0 = delfem2::RootInInterval_Bisection(
          intvl.first, intvl.second,
          coeff.data(), coeff.size(), 15);
        if (b0) {
          EXPECT_NEAR(y0, x0, 1.0e-4);
        } else {
          EXPECT_NEAR(y0, x1, 1.0e-4);
        }
      }
    }
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
    const double b = -a * (x0 + x1 + x2);
    const double c = a * (x0 * x1 + x1 * x2 + x2 * x0);
    const double d = -a * x0 * x1 * x2;
    const auto coeff = std::array<double, 4>{d, c, b, a};
    double strum[4][4];
    delfem2::StrumSequenceOfPolynomial<4>(strum, std::array<double,4>{d,c,b,a}.data());
    int nroot1 = delfem2::StrumNumber<4>(0, strum) - delfem2::StrumNumber<4>(1, strum);
    { // number of roots
      int nroot0 = 0;
      if (x0 > 0 && x0 < 1) { nroot0++; }
      if (x1 > 0 && x1 < 1) { nroot0++; }
      if (x2 > 0 && x2 < 1) { nroot0++; }
      EXPECT_EQ(nroot0, nroot1);
    }
    {
      std::vector<std::pair<double, double>> intvls = delfem2::RootInterval_StrumSequence<4>(0, 1, strum);
      EXPECT_EQ(intvls.size(), nroot1);
      for (auto intvl: intvls) {
        bool b0 = (intvl.first - x0) * (intvl.second - x0) < 0;
        bool b1 = (intvl.first - x1) * (intvl.second - x1) < 0;
        bool b2 = (intvl.first - x2) * (intvl.second - x2) < 0;
        EXPECT_TRUE(b0 || b1 || b2);
        // the number of roots should be 1 in the interval
        int nr = delfem2::StrumNumber<4>(intvl.first, strum) - delfem2::StrumNumber<4>(intvl.second, strum);
        EXPECT_EQ(nr, 1);
        double v0 = delfem2::Eval_Polynomial(intvl.first, coeff.data(), 4);
        double v1 = delfem2::Eval_Polynomial(intvl.second, coeff.data(), 4);
        EXPECT_LT(v0 * v1, 0.);
        double y0 = delfem2::RootInInterval_Bisection(
          intvl.first, intvl.second,
          coeff.data(), coeff.size(), 15);
        if (b0) {
          EXPECT_NEAR(y0, x0, 1.0e-4);
        } else if (b1) {
          EXPECT_NEAR(y0, x1, 1.0e-4);
        } else {
          EXPECT_NEAR(y0, x2, 1.0e-4);
        }
      }
    }
  }
}

TEST(polynomial_root, quartic_4roots) {
  std::mt19937 rndeng(std::random_device{}());
  std::uniform_real_distribution<double> dist_m1p1(-1, +1);
  for (unsigned int itr = 0; itr < 1000; ++itr) {
    const double a = dist_m1p1(rndeng);
    const double x0 = dist_m1p1(rndeng) * 2;
    const double x1 = dist_m1p1(rndeng) * 2;
    const double x2 = dist_m1p1(rndeng) * 2;
    const double x3 = dist_m1p1(rndeng) * 2;
    const double b = -a * (x0 + x1 + x2 + x3);
    const double c = +a * (x0 * x1 + x0 * x2 + x0 * x3 + x1 * x2 + x1 * x3 + x2 * x3);
    const double d = -a * (x0 * x1 * x2 + x1 * x2 * x3 + x2 * x3 * x0 + x3 * x0 * x1);
    const double e = +a * x0 * x1 * x2 * x3;
    const auto coeff = std::array<double, 5>{e, d, c, b, a};
    double strum[5][5];
    delfem2::StrumSequenceOfPolynomial<5>(strum, coeff.data());
    int nroot1 = delfem2::StrumNumber<5>(0, strum) - delfem2::StrumNumber<5>(1, strum);
    {
      int nroot0 = 0;
      if (x0 > 0 && x0 < 1) { nroot0++; }
      if (x1 > 0 && x1 < 1) { nroot0++; }
      if (x2 > 0 && x2 < 1) { nroot0++; }
      if (x3 > 0 && x3 < 1) { nroot0++; }
      EXPECT_EQ(nroot0, nroot1);
    }
    {
      std::vector<std::pair<double, double>> intvls = delfem2::RootInterval_StrumSequence<5>(0, 1, strum);
      EXPECT_EQ(intvls.size(), nroot1);
      for (auto intvl: intvls) {
        bool b0 = (intvl.first - x0) * (intvl.second - x0) < 0;
        bool b1 = (intvl.first - x1) * (intvl.second - x1) < 0;
        bool b2 = (intvl.first - x2) * (intvl.second - x2) < 0;
        bool b3 = (intvl.first - x3) * (intvl.second - x3) < 0;
        EXPECT_TRUE(b0 || b1 || b2 || b3);
        const int nr = delfem2::StrumNumber<5>(intvl.first, strum) - delfem2::StrumNumber<5>(intvl.second, strum);
        EXPECT_EQ(nr, 1);
        const double v0 = delfem2::Eval_Polynomial(intvl.first, coeff.data(), 5);
        const double v1 = delfem2::Eval_Polynomial(intvl.second, coeff.data(), 5);
        EXPECT_LT(v0 * v1, 0.);
      }
    }
  }
}

TEST(polynomial_root, quintic_5roots) {
  std::mt19937 rndeng(std::random_device{}());
  std::uniform_real_distribution<double> dist_m1p1(-1, +1);
  for (unsigned int itr = 0; itr < 1000; ++itr) {
    const double a = dist_m1p1(rndeng);
    const double x0 = dist_m1p1(rndeng) * 2;
    const double x1 = dist_m1p1(rndeng) * 2;
    const double x2 = dist_m1p1(rndeng) * 2;
    const double x3 = dist_m1p1(rndeng) * 2;
    const double x4 = dist_m1p1(rndeng) * 2;
    const double b = -a * (x0 + x1 + x2 + x3 + x4);
    const double c = +a * (
      + x0 * (x1 + x2 + x3 + x4)
      + x1 * (x2 + x3 + x4)
      + x2 * (x3 + x4)
      + x3 * x4 );
    const double d = -a * (
      + x0 * x1 * (x2 + x3 + x4)
      + x0 * x2 * (x3 + x4)
      + x0 * x3 * x4
      + x1 * (x2 * x3 + x2 * x4 + x3 * x4)
      + x2 * x3 * x4);
    const double e = +a * (
      + x1 * x2 * x3 * x4
      + x0 * x2 * x3 * x4
      + x0 * x1 * x3 * x4
      + x0 * x1 * x2 * x4
      + x0 * x1 * x2 * x3);
    const double f = -a * x0 * x1 * x2 * x3 * x4;
    const auto coeff = std::array<double, 6>{f, e, d, c, b, a};
    // F(x)
    // = a*(x-x0)*(x-x1)*(x-x2)*(x-x3)*(x-x4)
    // = a * x^5 + b * x^4 + c * x^3 + d * x^2 + e * x + f
    double strum[6][6];
    delfem2::StrumSequenceOfPolynomial<6>(strum, coeff.data() );
    int nroot1 = delfem2::StrumNumber<6>(0, strum) - delfem2::StrumNumber<6>(1, strum);
    {  // check the number of roots
      int nroot0 = 0;
      if (x0 > 0 && x0 < 1) { nroot0++; }
      if (x1 > 0 && x1 < 1) { nroot0++; }
      if (x2 > 0 && x2 < 1) { nroot0++; }
      if (x3 > 0 && x3 < 1) { nroot0++; }
      if (x4 > 0 && x4 < 1) { nroot0++; }
      EXPECT_EQ(nroot0, nroot1);
    }
    {
      std::vector<std::pair<double, double>> intvls = delfem2::RootInterval_StrumSequence<6>(0, 1, strum);
      EXPECT_EQ(intvls.size(), nroot1);
      for (auto intvl: intvls) {
        const bool b0 = (intvl.first - x0) * (intvl.second - x0) < 0;
        const bool b1 = (intvl.first - x1) * (intvl.second - x1) < 0;
        const bool b2 = (intvl.first - x2) * (intvl.second - x2) < 0;
        const bool b3 = (intvl.first - x3) * (intvl.second - x3) < 0;
        const bool b4 = (intvl.first - x4) * (intvl.second - x4) < 0;
        EXPECT_TRUE(b0 || b1 || b2 || b3 || b4);
        const int nr = delfem2::StrumNumber<6>(intvl.first, strum) - delfem2::StrumNumber<6>(intvl.second, strum);
        EXPECT_EQ(nr, 1);
        const double v0 = delfem2::Eval_Polynomial(intvl.first, coeff.data(), 6);
        const double v1 = delfem2::Eval_Polynomial(intvl.second, coeff.data(), 6);
        EXPECT_LT(v0 * v1, 0.);
      }
    }
  }
}
