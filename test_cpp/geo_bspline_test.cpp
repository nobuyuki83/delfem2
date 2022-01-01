//
// Created by Nobuyuki Umetani on 2021/12/09.
//

#include <random>

#include "gtest/gtest.h"

#include "delfem2/vec2.h"
#include "delfem2/geo_curve_ndegree.h"
#include "delfem2/geo_curve_quadratic.h"
#include "delfem2/geo_curve_cubic.h"
#include "delfem2/geo_polyline.h"

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

// compare quadratic-specific function and general function for BSpline
TEST(bspline, quadratic_near) {
  namespace dfm2 = delfem2;
  std::mt19937 rndeng(std::random_device{}());
  std::uniform_real_distribution<double> dist_01(0, 1);
  for (unsigned int itr = 0; itr < 1000; ++itr) {
    std::vector<dfm2::CVec2d> cps = {
      dfm2::CVec2d(dist_01(rndeng), dist_01(rndeng)),
      dfm2::CVec2d(dist_01(rndeng), dist_01(rndeng))};
    for (int ip = 0; ip < 10; ++ip) {  //  increasing number of points
      cps.emplace_back(dist_01(rndeng), dist_01(rndeng));
      const dfm2::CVec2d scr(0.5, 0.5);
      double t = dfm2::Nearest_QuadraticBSplineCurve(cps, scr);
      const dfm2::CVec2d p0 = delfem2::Sample_QuadraticBsplineCurve(t, cps);
      std::vector<dfm2::CVec2d> polyline;
      {
        double t_max = cps.size() - 2;
        unsigned int N = 1000;
        for (unsigned int i = 0; i < N; ++i) {
          double t = t_max * static_cast<double>(i) / static_cast<double>(N-1);
          dfm2::CVec2d q0 = delfem2::Sample_QuadraticBsplineCurve(t,cps);
          polyline.push_back(q0);
        }
      }
      const double param = delfem2::Nearest_Polyline(polyline, scr);
      const dfm2::CVec2d p1 = delfem2::Sample_Polyline(polyline, param);
      EXPECT_NEAR((p0-scr).norm(), (p1-scr).norm(), 1.0e-4);
    }
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

