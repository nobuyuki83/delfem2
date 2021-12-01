//
// Created by Nobuyuki Umetani on 2021/11/30.
//

#include <Eigen/Core>
#include <random>

#include "gtest/gtest.h"
#include "delfem2/geo_polyline2.h"
#include "delfem2/geo_bezier_quadratic.h"

TEST(cad,quadratic_bezier0) {
  namespace dfm2 = delfem2;
  std::mt19937 rndeng(std::random_device{}());
  std::uniform_real_distribution<double> dist_01(0,1);
  for(unsigned int itr=0;itr<1000;++itr){
    Eigen::Vector2d p0, p1, p2;
    if( itr == 0 ) {
      p0 = Eigen::Vector2d(0, 0);
      p1 = Eigen::Vector2d(1, 1);
      p2 = Eigen::Vector2d(2, 0);
    }
    else if( itr == 1 ) {
      p0 = Eigen::Vector2d(0, 0);
      p1 = Eigen::Vector2d(1, 1);
      p2 = Eigen::Vector2d(2, 2);
    }
    else {
      p0 = Eigen::Vector2d(dist_01(rndeng), dist_01(rndeng));
      p1 = Eigen::Vector2d(dist_01(rndeng), dist_01(rndeng));
      p2 = Eigen::Vector2d(dist_01(rndeng), dist_01(rndeng));
    }
    std::vector<Eigen::Vector2d> polyline;
    dfm2::Polyline_BezierQuadratic(polyline, 1000, p0, p1, p2);
    {  // length
      double v1 = dfm2::Length_QuadraticBezierCurve_Analytic<Eigen::Vector2d>(
        p0, p1, p2);
      double v0 = dfm2::Length_Polyline<Eigen::Vector2d>(polyline);
      double v2 = dfm2::Length_QuadraticBezierCurve_Quadrature<Eigen::Vector2d>(
        p0,p1,p2,3);
      EXPECT_NEAR(v1,v0,1.0e-4);
      EXPECT_NEAR(v1,v2,5.0e-2);
    }
    {  // area
      double a0 = dfm2::Area_QuadraticBezierCurve2(p0, p1, p2);
      double a1 = dfm2::Area_Polyline2(polyline);
      EXPECT_NEAR(a1,a0,1.0e-4);
    }
    {  // nearest
      Eigen::Vector2d scr(0.5, 0.5);
      const auto [i0,l0] = dfm2::FindNearestPointInPolyline(polyline, scr);
      const Eigen::Vector2d v0 = dfm2::PositionInPolyline(polyline, i0, l0);
      const double t1 = dfm2::Nearest_QuadraticBezierCurve<Eigen::Vector2d>(
        scr,
        p0,p1,p2,100,3);
      const Eigen::Vector2d v1 = dfm2::PointOnQuadraticBezierCurve(t1,p0,p1,p2);
      EXPECT_LT( (v0-v1).norm(), 1.0e-2 );
    }
    {  // bb
      auto bb0 = dfm2::AABB_QuadraticBezierCurve<2>(p0,p1,p2);
      auto bb1 = dfm2::AABB_Polyline<2>(polyline);
      EXPECT_NEAR(bb0[0],bb1[0],1.0e-3);
      EXPECT_NEAR(bb0[1],bb1[1],1.0e-3);
      EXPECT_NEAR(bb0[2],bb1[2],1.0e-3);
      EXPECT_NEAR(bb0[3],bb1[3],1.0e-3);
    }
  }
}