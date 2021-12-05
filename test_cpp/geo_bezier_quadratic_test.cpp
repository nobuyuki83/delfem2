//
// Created by Nobuyuki Umetani on 2021/12/01.
//

#include <string>
#include <random>

#include "gtest/gtest.h"

#include "delfem2/geo_bezier_cubic.h"
#include "delfem2/geo_bezier_quadratic.h"
#include "delfem2/geo_polyline2.h"

TEST(bezier_quadratic, test0) {
  namespace dfm2 = delfem2;
  std::mt19937 rndeng(std::random_device{}());
  std::uniform_real_distribution<double> dist_01(0, 1);
  for(unsigned int itr=0;itr<1000;++itr){
    dfm2::CVec2d p0, p1, p2;
    p0 = dfm2::CVec2d(dist_01(rndeng), dist_01(rndeng));
    p1 = dfm2::CVec2d(dist_01(rndeng), dist_01(rndeng));
    p2 = dfm2::CVec2d(dist_01(rndeng), dist_01(rndeng));
    std::vector<dfm2::CVec2d> polyline;
    dfm2::Polyline_BezierQuadratic(polyline, 100, p0, p1, p2);
    {  // length
      double v1 = dfm2::Length_QuadraticBezierCurve_Analytic<dfm2::CVec2d>(
        p0, p1, p2);
      double v0 = dfm2::Length_Polyline<dfm2::CVec2d>(polyline);
      double v2 = dfm2::Length_QuadraticBezierCurve_Quadrature<dfm2::CVec2d>(
        p0,p1,p2,3);
      EXPECT_NEAR(v1,v0,2.0e-4);
      EXPECT_NEAR(v1,v2,0.2);
    }
    {  // area
      double a0 = dfm2::Area_QuadraticBezierCurve2(p0, p1, p2);
      double a1 = dfm2::Area_Polyline2(polyline);
      EXPECT_NEAR(a1,a0,2.0e-4);
    }
    {  // nearest
      dfm2::CVec2d scr(0.5, 0.5);
      const auto [i0,l0] = dfm2::FindNearestPointInPolyline(polyline, scr);
      const dfm2::CVec2d v0 = dfm2::PositionInPolyline(polyline, i0, l0);
      const double t1 = dfm2::Nearest_QuadraticBezierCurve<dfm2::CVec2d>(
        scr,
        p0,p1,p2,100,5);
      const dfm2::CVec2d v1 = dfm2::PointOnQuadraticBezierCurve(t1,p0,p1,p2);
      EXPECT_LT( (v0-v1).norm(), 0.02 );
    }
    {  // bb
      auto bb0 = dfm2::AABB_QuadraticBezierCurve<2>(p0,p1,p2);
      auto bb1 = dfm2::AABB_Polyline<2>(polyline);
      EXPECT_NEAR(bb0[0], bb1[0], 1.0e-3);
      EXPECT_NEAR(bb0[1], bb1[1], 1.0e-3);
      EXPECT_NEAR(bb0[2], bb1[2], 1.0e-3);
      EXPECT_NEAR(bb0[3], bb1[3], 1.0e-3);
    }
  }
}

