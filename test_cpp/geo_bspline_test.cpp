//
// Created by Nobuyuki Umetani on 2021/12/09.
//

#include <random>

#include "gtest/gtest.h"

#include "delfem2/vec2.h"
#include "delfem2/geo_bspline.h"
#include "delfem2/geo_bezier_quadratic.h"
#include "delfem2/geo_bezier_cubic.h"
#include "delfem2/geo_polyline2.h"

TEST(bspline, quadratic_bezier) {
  namespace dfm2 = delfem2;
  std::mt19937 rndeng(std::random_device{}());
  std::uniform_real_distribution<double> dist_01(0, 1);
  for (unsigned int itr = 0; itr < 1000; ++itr) {
    const dfm2::CVec2d p0 = dfm2::CVec2d(dist_01(rndeng), dist_01(rndeng));
    const dfm2::CVec2d p1 = dfm2::CVec2d(dist_01(rndeng), dist_01(rndeng));
    const dfm2::CVec2d p2 = dfm2::CVec2d(dist_01(rndeng), dist_01(rndeng));
    int ndiv = 10;
    for (unsigned int i = 0; i < ndiv + 1; ++i) {
      double t = double(i) / double(ndiv);
      const dfm2::CVec2d pt = delfem2::PointOnQuadraticBezierCurve(t, p0, p1, p2);
      const auto pta = delfem2::Sample_QuadraticBsplineCurve<dfm2::CVec2d>(t, {p0, p1, p2});
      const auto ptb = delfem2::Sample_BsplineCurve<dfm2::CVec2d, 3>(t, {p0, p1, p2});
      EXPECT_LT((pt - pta).norm(), 1.0e-10);
      EXPECT_LT((pt - ptb).norm(), 1.0e-10);
    }
  }
}

// compare quadratic-specific function and general function for BSpline
TEST(bspline, quadratic_general) {
  namespace dfm2 = delfem2;
  std::mt19937 rndeng(std::random_device{}());
  std::uniform_real_distribution<double> dist_01(0, 1);
  for(unsigned int itr=0;itr<1000;++itr) {
    std::vector<dfm2::CVec2d> poly = {
      dfm2::CVec2d(dist_01(rndeng), dist_01(rndeng)),
      dfm2::CVec2d(dist_01(rndeng), dist_01(rndeng))};
    for (int ip = 0; ip < 10; ++ip) {  //  increasing number of points
      poly.emplace_back(dist_01(rndeng), dist_01(rndeng));
      const auto t_end = static_cast<double>(poly.size()-2);
      double t0 = t_end * 0.8 * dist_01(rndeng) + 0.1;
      const double eps = 1.0e-5;
      {  // match start
        const auto p0a = delfem2::Sample_QuadraticBsplineCurve<dfm2::CVec2d>(0, poly);
        EXPECT_LT((p0a - poly[0]).norm(), 1.0e-10);
        const auto p0b = delfem2::Sample_BsplineCurve<dfm2::CVec2d, 3>(0 + eps, poly);
        const auto dp0 = delfem2::Sample_BsplineDerivative<dfm2::CVec2d, 3>(0, poly);
        EXPECT_LT(((p0b - p0a) / eps - dp0).norm(), 1.0e-4);
      }
      {  // match end
        const auto pea = delfem2::Sample_QuadraticBsplineCurve<dfm2::CVec2d>(t_end, poly);
        EXPECT_LT((pea - poly[poly.size()-1]).norm(), 1.0e-10);
        const auto peb = delfem2::Sample_BsplineCurve<dfm2::CVec2d, 3>(t_end - eps, poly);
        const auto dpe = delfem2::Sample_BsplineDerivative<dfm2::CVec2d, 3>(t_end, poly);
        EXPECT_LT(((pea - peb) / eps - dpe).norm(), 1.0e-4);
      }
      {  // match general
        const auto p0a = delfem2::Sample_BsplineCurve<dfm2::CVec2d, 3>(t0, poly);
        const auto p0b = delfem2::Sample_QuadraticBsplineCurve<dfm2::CVec2d>(t0, poly);
        EXPECT_LT((p0a - p0b).norm(), 1.0e-10);
        const auto p1b = delfem2::Sample_BsplineCurve<dfm2::CVec2d, 3>(t0 - eps, poly);
        const auto dp = delfem2::Sample_BsplineDerivative<dfm2::CVec2d, 3>(t0, poly);
        EXPECT_LT(((p0a - p1b) / eps - dp).norm(), 1.0e-4);
      }
    }
  }
}

TEST(bspline, quadratic_curve) {
  namespace dfm2 = delfem2;
  std::mt19937 rndeng(std::random_device{}());
  std::uniform_real_distribution<double> dist_01(0, 1);
  std::uniform_real_distribution<double> dist_m1p1(-1, 1);
  for (unsigned int itr = 0; itr < 1000; ++itr) {
    const dfm2::CVec2d cp[3] = {
      dfm2::CVec2d(dist_01(rndeng), dist_01(rndeng)),
      dfm2::CVec2d(dist_01(rndeng), dist_01(rndeng)),
      dfm2::CVec2d(dist_01(rndeng), dist_01(rndeng)) };
    const double coeff[3][3] = {
      {dist_m1p1(rndeng), dist_m1p1(rndeng), dist_m1p1(rndeng)},
      {dist_m1p1(rndeng), dist_m1p1(rndeng), dist_m1p1(rndeng)},
      {dist_m1p1(rndeng), dist_m1p1(rndeng), dist_m1p1(rndeng)} };
    double l0 = delfem2::Length_ParametricCurve_Quadratic(coeff,cp);
    std::vector<dfm2::CVec2d> poly;
    const unsigned int N = 1000;
    for(unsigned int ip=0;ip<N+1;++ip){
      double t = static_cast<double>(ip) / static_cast<double>(N);
      const double w0 = coeff[0][0] + coeff[0][1] * t + coeff[0][2] * t * t;
      const double w1 = coeff[1][0] + coeff[1][1] * t + coeff[1][2] * t * t;
      const double w2 = coeff[2][0] + coeff[2][1] * t + coeff[2][2] * t * t;
      poly.push_back(cp[0] * w0 + cp[1] * w1 + cp[2] * w2);
    }
    double l1 = delfem2::Length_Polyline(poly);
    EXPECT_NEAR(l0,l1,1.0e-5);
  }
}

TEST(bspline, cubic_bezier) {
  namespace dfm2 = delfem2;
  std::mt19937 rndeng(std::random_device{}());
  std::uniform_real_distribution<double> dist_01(0, 1);
  for (unsigned int itr = 0; itr < 1000; ++itr) {
    std::vector<dfm2::CVec2d> poly = {
      dfm2::CVec2d(dist_01(rndeng), dist_01(rndeng)),
      dfm2::CVec2d(dist_01(rndeng), dist_01(rndeng)),
      dfm2::CVec2d(dist_01(rndeng), dist_01(rndeng)),
      dfm2::CVec2d(dist_01(rndeng), dist_01(rndeng)) };
    int ndiv = 10;
    for(unsigned int i=0;i<ndiv+1;++i) {
      double t = double(i)/double(ndiv);
      const dfm2::CVec2d pt = delfem2::PointOnCubicBezierCurve(t, poly[0], poly[1], poly[2], poly[3]);
      const auto pta = delfem2::Sample_CubicBsplineCurve<dfm2::CVec2d>(t, poly);
      const auto ptb = delfem2::Sample_BsplineCurve<dfm2::CVec2d,4>(t, poly);
      EXPECT_LT( (pt-ptb).norm(), 1.0e-10 );
      EXPECT_LT( (pt-pta).norm(), 1.0e-10 );
    }
  }
}


TEST(bspline, cubic_general) {
  namespace dfm2 = delfem2;
  std::mt19937 rndeng(std::random_device{}());
  std::uniform_real_distribution<double> dist_01(0, 1);
  for (unsigned int itr = 0; itr < 1000; ++itr) {
    std::vector<dfm2::CVec2d> poly = {
      dfm2::CVec2d(dist_01(rndeng), dist_01(rndeng)),
      dfm2::CVec2d(dist_01(rndeng), dist_01(rndeng)),
      dfm2::CVec2d(dist_01(rndeng), dist_01(rndeng)) };
    for(int ip=0;ip<10;++ip) {
      poly.emplace_back(dist_01(rndeng), dist_01(rndeng));
      const auto t_end = static_cast<double>(poly.size()-3);
      const double eps = 1.0e-5;
      { // start point match
        const auto p0a = delfem2::Sample_CubicBsplineCurve<dfm2::CVec2d>(0, poly);
        const auto p0b = delfem2::Sample_BsplineCurve<dfm2::CVec2d, 4>(0, poly);
        assert((p0a - poly[0]).norm() < 1.0e-10);
        assert((p0b - poly[0]).norm() < 1.0e-10);
        const auto p1b = delfem2::Sample_BsplineCurve<dfm2::CVec2d,4>(0+eps, poly);
        const auto dp = delfem2::Sample_BsplineDerivative<dfm2::CVec2d,4>(0,poly);
        EXPECT_LT( ((p1b-p0b)/eps-dp).norm(), 1.0e-4 );
      }
      { // end point match
        const auto pea = delfem2::Sample_CubicBsplineCurve<dfm2::CVec2d>(t_end, poly);
        const auto peb = delfem2::Sample_BsplineCurve<dfm2::CVec2d, 4>(t_end, poly);
        assert((pea - poly[poly.size() - 1]).norm() < 1.0e-10);
        assert((peb - poly[poly.size() - 1]).norm() < 1.0e-10);
        const auto pec = delfem2::Sample_BsplineCurve<dfm2::CVec2d,4>(t_end-eps, poly);
        const auto dp = delfem2::Sample_BsplineDerivative<dfm2::CVec2d,4>(t_end,poly);
        EXPECT_LT( ((peb-pec)/eps-dp).norm(), 1.0e-4 );
      }
      { // derivative at middle
        double t0 = dist_01(rndeng)*t_end*0.8 + 0.1;
        const auto p0a = delfem2::Sample_BsplineCurve<dfm2::CVec2d,4>(t0-eps, poly);
        const auto p0b = delfem2::Sample_BsplineCurve<dfm2::CVec2d,4>(t0, poly);
        const auto pta = delfem2::Sample_CubicBsplineCurve<dfm2::CVec2d>(t0, poly);
        const auto t = delfem2::Sample_BsplineDerivative<dfm2::CVec2d,4>(t0,poly);
        EXPECT_LT( ((p0b-p0a)/eps-t).norm(), 1.0e-4 );
        EXPECT_LT((p0b - pta).norm(), 1.0e-10);
      }
    }
  }
}

