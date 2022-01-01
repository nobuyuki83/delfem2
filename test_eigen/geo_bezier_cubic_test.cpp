//
// Created by Nobuyuki Umetani on 2021/11/30.
//

#include <Eigen/Core>
#include <random>

#include "gtest/gtest.h"
#include "delfem2/geo_curve_cubic.h"
#include "delfem2/geo_polyline.h"

TEST(cad, cubic_bezier0) {
  namespace dfm2 = delfem2;
  std::mt19937 rndeng(std::random_device{}());
  std::uniform_real_distribution<double> dist_01(0, 1);
  for (int itr = 0; itr < 1000; ++itr) {
    Eigen::Vector2d p0, p1, p2, p3;
    p0 = Eigen::Vector2d(dist_01(rndeng), dist_01(rndeng));
    p1 = Eigen::Vector2d(dist_01(rndeng), dist_01(rndeng));
    p2 = Eigen::Vector2d(dist_01(rndeng), dist_01(rndeng));
    p3 = Eigen::Vector2d(dist_01(rndeng), dist_01(rndeng));
    //
    std::vector<Eigen::Vector2d> polyline;
    dfm2::Polyline_BezierCubic(polyline, 1000, p0, p1, p2, p3);
    {
      double v0 = dfm2::Length_Polyline<Eigen::Vector2d>(polyline);
      double v2 = dfm2::Length_CubicBezierCurve_Quadrature<Eigen::Vector2d>(p0, p1, p2, p3, 3);
      EXPECT_NEAR(v0, v2, 0.2);
    }
    {
      double a0 = dfm2::Area_CubicBezierCurve2(p0, p1, p2, p3);
      double a1 = dfm2::Area_Polyline2(polyline);
      EXPECT_NEAR(a0, a1, 1.0e-5);
    }
    {  // nearest
      Eigen::Vector2d scr(0.5, 0.5);
      const double param = dfm2::Nearest_Polyline(polyline, scr);
      const Eigen::Vector2d v0 = dfm2::Sample_Polyline(polyline, param);
      const double t1 = dfm2::Nearest_CubicBezierCurve<Eigen::Vector2d>(
          scr,
          p0, p1, p2, p3, 100, 3);
      const Eigen::Vector2d v1 = dfm2::PointOnCubicBezierCurve(t1, p0, p1, p2, p3);
      EXPECT_NEAR((v0 - scr).norm(), (v1 - scr).norm(), 0.002);
    }
    {  // bb
      auto bb0 = dfm2::AABB_CubicBezierCurve(p0, p1, p2, p3);
      auto bb1 = dfm2::AABB_Polyline(polyline);
      EXPECT_NEAR(bb0[0], bb1[0], 1.0e-3);
      EXPECT_NEAR(bb0[1], bb1[1], 1.0e-3);
      EXPECT_NEAR(bb0[2], bb1[2], 1.0e-3);
      EXPECT_NEAR(bb0[3], bb1[3], 1.0e-3);
    }
  }
}
