//
// Created by Nobuyuki Umetani on 2021/12/09.
//

#include <random>

#include "gtest/gtest.h"

#include "delfem2/vec2.h"
#include "delfem2/geo_bspline.h"
#include "delfem2/geo_bezier_quadratic.h"
#include "delfem2/geo_bezier_cubic.h"

TEST(bspline, quadratic) {
  namespace dfm2 = delfem2;
  std::mt19937 rndeng(std::random_device{}());
  std::uniform_real_distribution<double> dist_01(0, 1);
  for(unsigned int itr=0;itr<1000;++itr) {
    const dfm2::CVec2d p0 = dfm2::CVec2d(dist_01(rndeng), dist_01(rndeng));
    const dfm2::CVec2d p1 = dfm2::CVec2d(dist_01(rndeng), dist_01(rndeng));
    const dfm2::CVec2d p2 = dfm2::CVec2d(dist_01(rndeng), dist_01(rndeng));
    {
      const auto p0a = delfem2::Sample_QuadraticBsplineCurve<dfm2::CVec2d>(0, {p0, p1, p2});
      const auto p0b = delfem2::Sample_BsplineCurve<dfm2::CVec2d, 3>(0, {p0, p1, p2});
      assert((p0a - p0).norm() < 1.0e-10);
      assert((p0b - p0).norm() < 1.0e-10);
    }
    {
      const auto p2a = delfem2::Sample_QuadraticBsplineCurve<dfm2::CVec2d>(1, {p0, p1, p2});
      const auto p2b = delfem2::Sample_BsplineCurve<dfm2::CVec2d, 3>(1, {p0, p1, p2});
      assert((p2a - p2).norm() < 1.0e-10);
      assert((p2b - p2).norm() < 1.0e-10);
    }
    int ndiv = 10;
    for(unsigned int i=0;i<ndiv+1;++i) {
      double t = double(i)/double(ndiv);
      const dfm2::CVec2d pt = delfem2::PointOnQuadraticBezierCurve(t, p0, p1, p2);
      const auto pta = delfem2::Sample_QuadraticBsplineCurve<dfm2::CVec2d>(t, {p0, p1, p2});
      const auto ptb = delfem2::Sample_BsplineCurve<dfm2::CVec2d,3>(t, {p0, p1, p2});
      EXPECT_LT( (pt-pta).norm(), 1.0e-10 );
      EXPECT_LT( (pt-ptb).norm(), 1.0e-10 );
    }
  }
  for(unsigned int itr=0;itr<1000;++itr) {
    const dfm2::CVec2d p0 = dfm2::CVec2d(dist_01(rndeng), dist_01(rndeng));
    const dfm2::CVec2d p1 = dfm2::CVec2d(dist_01(rndeng), dist_01(rndeng));
    const dfm2::CVec2d p2 = dfm2::CVec2d(dist_01(rndeng), dist_01(rndeng));
    const dfm2::CVec2d p3 = dfm2::CVec2d(dist_01(rndeng), dist_01(rndeng));
    {
      const auto p0a = delfem2::Sample_QuadraticBsplineCurve<dfm2::CVec2d>(0, {p0, p1, p2, p3});
      EXPECT_LT((p0a - p0).norm(), 1.0e-10);
      const auto p3a = delfem2::Sample_QuadraticBsplineCurve<dfm2::CVec2d>(2, {p0, p1, p2, p3});
      EXPECT_LT((p3a - p3).norm(), 1.0e-10);
    }
    double t0 = dist_01(rndeng)*1.6 + 0.2;
    {
      const auto p0a = delfem2::Sample_BsplineCurve<dfm2::CVec2d,3>(t0, {p0, p1, p2, p3});
      const auto p0b = delfem2::Sample_QuadraticBsplineCurve<dfm2::CVec2d>(t0, {p0, p1, p2, p3});
      EXPECT_LT( (p0a-p0b).norm(), 1.0e-10 );
    }
    double eps = 1.0e-5;
    {
      double eps = 1.0e-5;
      const auto p0a = delfem2::Sample_BsplineCurve<dfm2::CVec2d,3>(0, {p0, p1, p2, p3});
      const auto p0b = delfem2::Sample_BsplineCurve<dfm2::CVec2d,3>(0+eps, {p0, p1, p2, p3});
      const auto t0 = delfem2::Sample_BsplineDerivative<dfm2::CVec2d,3>(0,{p0, p1, p2, p3});
      EXPECT_LT( ((p0b-p0a)/eps-t0).norm(), 1.0e-4 );
    }
    {
      double eps = 1.0e-5;
      const auto p0a = delfem2::Sample_BsplineCurve<dfm2::CVec2d,3>(2.-eps, {p0, p1, p2, p3});
      const auto p0b = delfem2::Sample_BsplineCurve<dfm2::CVec2d,3>(2, {p0, p1, p2, p3});
      const auto t0 = delfem2::Sample_BsplineDerivative<dfm2::CVec2d,3>(2,{p0, p1, p2, p3});
      EXPECT_LT( ((p0b-p0a)/eps-t0).norm(), 1.0e-4 );
    }
    {
      const auto p0a = delfem2::Sample_BsplineCurve<dfm2::CVec2d,3>(t0-eps, {p0, p1, p2, p3});
      const auto p0b = delfem2::Sample_BsplineCurve<dfm2::CVec2d,3>(t0, {p0, p1, p2, p3});
      const auto s0 = delfem2::Sample_BsplineDerivative<dfm2::CVec2d,3>(t0,{p0, p1, p2, p3});
      EXPECT_LT( ((p0b-p0a)/eps-s0).norm(), 1.0e-4 );
    }
  }
}

TEST(bspline, cubic) {
  namespace dfm2 = delfem2;
  std::mt19937 rndeng(std::random_device{}());
  std::uniform_real_distribution<double> dist_01(0, 1);
  for (unsigned int itr = 0; itr < 1000; ++itr) {
    const dfm2::CVec2d p0 = dfm2::CVec2d(dist_01(rndeng), dist_01(rndeng));
    const dfm2::CVec2d p1 = dfm2::CVec2d(dist_01(rndeng), dist_01(rndeng));
    const dfm2::CVec2d p2 = dfm2::CVec2d(dist_01(rndeng), dist_01(rndeng));
    const dfm2::CVec2d p3 = dfm2::CVec2d(dist_01(rndeng), dist_01(rndeng));
    {
      const auto p0b = delfem2::Sample_BsplineCurve<dfm2::CVec2d, 4>(0, {p0, p1, p2, p3});
      assert((p0b - p0).norm() < 1.0e-10);
    }
    {
      const auto p3b = delfem2::Sample_BsplineCurve<dfm2::CVec2d, 4>(1, {p0, p1, p2, p3});
      assert((p3b - p3).norm() < 1.0e-10);
    }
    int ndiv = 10;
    for(unsigned int i=0;i<ndiv+1;++i) {
      double t = double(i)/double(ndiv);
      const dfm2::CVec2d pt = delfem2::PointOnCubicBezierCurve(t, p0, p1, p2, p3);
      const auto ptb = delfem2::Sample_BsplineCurve<dfm2::CVec2d,4>(t, {p0, p1, p2, p3});
      EXPECT_LT( (pt-ptb).norm(), 1.0e-10 );
    }
    {
      double eps = 1.0e-5;
      const auto p0a = delfem2::Sample_BsplineCurve<dfm2::CVec2d,4>(0, {p0, p1, p2, p3});
      const auto p0b = delfem2::Sample_BsplineCurve<dfm2::CVec2d,4>(0+eps, {p0, p1, p2, p3});
      const auto t = delfem2::Sample_BsplineDerivative<dfm2::CVec2d,4>(0,{p0, p1, p2, p3});
      EXPECT_LT( ((p0b-p0a)/eps-t).norm(), 1.0e-4 );
    }
    {
      double eps = 1.0e-5;
      const auto p0a = delfem2::Sample_BsplineCurve<dfm2::CVec2d,4>(1.-eps, {p0, p1, p2, p3});
      const auto p0b = delfem2::Sample_BsplineCurve<dfm2::CVec2d,4>(1, {p0, p1, p2, p3});
      const auto t = delfem2::Sample_BsplineDerivative<dfm2::CVec2d,4>(1,{p0, p1, p2, p3});
      EXPECT_LT( ((p0b-p0a)/eps-t).norm(), 1.0e-4 );
    }
    {
      double t0 = dist_01(rndeng)*0.8 + 0.1;
      double eps = 1.0e-5;
      const auto p0a = delfem2::Sample_BsplineCurve<dfm2::CVec2d,4>(t0-eps, {p0, p1, p2, p3});
      const auto p0b = delfem2::Sample_BsplineCurve<dfm2::CVec2d,4>(t0, {p0, p1, p2, p3});
      const auto t = delfem2::Sample_BsplineDerivative<dfm2::CVec2d,4>(t0,{p0, p1, p2, p3});
      EXPECT_LT( ((p0b-p0a)/eps-t).norm(), 1.0e-4 );
    }
  }
}

