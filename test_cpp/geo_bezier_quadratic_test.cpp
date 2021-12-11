//
// Created by Nobuyuki Umetani on 2021/12/01.
//

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
    dfm2::Polyline_BezierQuadratic(polyline, 1000, p0, p1, p2);
    {  // length
      double v1 = dfm2::Length_QuadraticBezierCurve_Analytic<dfm2::CVec2d>(
        p0, p1, p2);
      double v0 = dfm2::Length_Polyline<dfm2::CVec2d>(polyline);
      double v2 = dfm2::Length_QuadraticBezierCurve_Quadrature<dfm2::CVec2d>(
        p0,p1,p2,3);
      double v3 = dfm2::Length_QuadraticBezierCurve_QuadratureSubdiv<dfm2::CVec2d,11>(
        p0,p1,p2,1.0e-7,12);
      EXPECT_NEAR(v1,v0,2.0e-4);
      EXPECT_NEAR(v1,v2,0.2);
      EXPECT_NEAR(v1,v3,1.0e-3);
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
    {
      std::vector<dfm2::CVec2d> poly0;
      {
        std::array<dfm2::CVec2d,6> r = dfm2::Split_QuadraticBezierCurve(p0,p1,p2, 0.3);
        delfem2::Polyline_BezierQuadratic(
          poly0, 4,
          r[0], r[1], r[2]);
        std::vector<dfm2::CVec2d> poly1;
        delfem2::Polyline_BezierQuadratic(
          poly1, 8,
          r[3], r[4], r[5]);
        poly0.resize(poly0.size() - 1);
        poly0.insert(poly0.end(), poly1.begin(), poly1.end());
      }
      std::vector<dfm2::CVec2d> poly2;
      delfem2::Polyline_BezierQuadratic(
        poly2, 11,
        p0, p1, p2);
      EXPECT_EQ(poly0.size(),poly2.size());
      for(unsigned int ip=0;ip<poly0.size();++ip){
        EXPECT_LT((poly0[ip]-poly2[ip]).norm(),1.0e-10);
      }
    }
  }
}

