//
// Created by Nobuyuki Umetani on 2021/12/09.
//

#include <random>

#include "gtest/gtest.h"

#include "delfem2/vec2.h"
#include "delfem2/geo_bspline.h"
#include "delfem2/geo_bezier_quadratic.h"

TEST(bspline, test0) {
  namespace dfm2 = delfem2;
  std::mt19937 rndeng(std::random_device{}());
  std::uniform_real_distribution<double> dist_01(0, 1);
  for(unsigned int itr=0;itr<1000;++itr) {
    const dfm2::CVec2d p0 = dfm2::CVec2d(dist_01(rndeng), dist_01(rndeng));
    const dfm2::CVec2d p1 = dfm2::CVec2d(dist_01(rndeng), dist_01(rndeng));
    const dfm2::CVec2d p2 = dfm2::CVec2d(dist_01(rndeng), dist_01(rndeng));
    const dfm2::CVec2d p0a = delfem2::Sample_QuadraticBsplineCurve<dfm2::CVec2d>(0, {p0, p1, p2});
    assert((p0a-p0).norm()<1.0e-10);
    const dfm2::CVec2d p2a = delfem2::Sample_QuadraticBsplineCurve<dfm2::CVec2d>(1, {p0, p1, p2});
    assert((p2a-p2).norm()<1.0e-10);
    int ndiv = 10;
    for(unsigned int i=0;i<ndiv+1;++i) {
      double t = double(i)/double(ndiv);
      const dfm2::CVec2d pt = delfem2::PointOnQuadraticBezierCurve(t, p0, p1, p2);
      const dfm2::CVec2d pta = delfem2::Sample_QuadraticBsplineCurve<dfm2::CVec2d>(t, {p0, p1, p2});
      EXPECT_LT( (pt-pta).norm(), 1.0e-10 );
    }
  }
  for(unsigned int itr=0;itr<1000;++itr) {
    const dfm2::CVec2d p0 = dfm2::CVec2d(dist_01(rndeng), dist_01(rndeng));
    const dfm2::CVec2d p1 = dfm2::CVec2d(dist_01(rndeng), dist_01(rndeng));
    const dfm2::CVec2d p2 = dfm2::CVec2d(dist_01(rndeng), dist_01(rndeng));
    const dfm2::CVec2d p3 = dfm2::CVec2d(dist_01(rndeng), dist_01(rndeng));
    const dfm2::CVec2d p0a = delfem2::Sample_QuadraticBsplineCurve<dfm2::CVec2d>(0, {p0, p1, p2,p3});
    assert((p0a - p0).norm() < 1.0e-10);
    const dfm2::CVec2d p3a = delfem2::Sample_QuadraticBsplineCurve<dfm2::CVec2d>(1, {p0, p1, p2,p3});
    assert((p3a - p3).norm() < 1.0e-10);
  }
}



