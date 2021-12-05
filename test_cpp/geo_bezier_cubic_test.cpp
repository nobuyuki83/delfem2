//
// Created by Nobuyuki Umetani on 2021/12/01.
//

#include <random>

#include "gtest/gtest.h"

#include "delfem2/geo_bezier_cubic.h"
#include "delfem2/geo_polyline2.h"

TEST(bezier_cubic, test0) {
  namespace dfm2 = delfem2;
  std::mt19937 rndeng(std::random_device{}());
  std::uniform_real_distribution<double> dist_01(0, 1);
  for (int itr = 0; itr < 1000; ++itr) {
    dfm2::CVec2d p0, p1, p2, p3;
    p0 = dfm2::CVec2d(dist_01(rndeng), dist_01(rndeng));
    p1 = dfm2::CVec2d(dist_01(rndeng), dist_01(rndeng));
    p2 = dfm2::CVec2d(dist_01(rndeng), dist_01(rndeng));
    p3 = dfm2::CVec2d(dist_01(rndeng), dist_01(rndeng));
    std::vector<dfm2::CVec2d> polyline;
    dfm2::Polyline_BezierCubic(polyline, 1000, p0, p1, p2, p3);
    {
      double v0 = dfm2::Length_Polyline<dfm2::CVec2d>(polyline);
      double v2 = dfm2::Length_CubicBezierCurve_Quadrature<dfm2::CVec2d>(p0, p1, p2, p3, 3);
      EXPECT_NEAR(v0, v2, 2.0e-1);
    }
    {
      double a0 = dfm2::Area_CubicBezierCurve2(p0, p1, p2, p3);
      double a1 = dfm2::Area_Polyline2(polyline);
      EXPECT_NEAR(a0, a1, 1.0e-5);
    }
    {  // nearest
      dfm2::CVec2d scr(0.5, 0.5);
      const auto[i0, l0] = dfm2::FindNearestPointInPolyline(polyline, scr);
      const dfm2::CVec2d v0 = dfm2::PositionInPolyline(polyline, i0, l0);
      const double t1 = dfm2::Nearest_CubicBezierCurve<dfm2::CVec2d>(
        scr,
        p0, p1, p2, p3, 100, 5);
      const dfm2::CVec2d v1 = dfm2::PointOnCubicBezierCurve(t1, p0, p1, p2, p3);
      EXPECT_LT((v0 - v1).norm(), 0.02);
    }
    {  // bb
      auto bb0 = dfm2::AABB_CubicBezierCurve<2>(p0, p1, p2, p3);
      auto bb1 = dfm2::AABB_Polyline<2>(polyline);
      EXPECT_NEAR(bb0[0], bb1[0], 1.0e-3);
      EXPECT_NEAR(bb0[1], bb1[1], 1.0e-3);
      EXPECT_NEAR(bb0[2], bb1[2], 1.0e-3);
      EXPECT_NEAR(bb0[3], bb1[3], 1.0e-3);
    }
  }
}
