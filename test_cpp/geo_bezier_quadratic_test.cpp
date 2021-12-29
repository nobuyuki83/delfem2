//
// Created by Nobuyuki Umetani on 2021/12/01.
//

#include <random>

#include "gtest/gtest.h"

#include "delfem2/geo_curve_quadratic.h"
#include "delfem2/geo_polyline.h"

TEST(bezier_quadratic, polyline2) {
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
      double v3 = dfm2::Length_QuadraticBezierCurve_QuadratureSubdiv<dfm2::CVec2d>(
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
      //
      const double t1 = dfm2::Nearest_QuadraticBezierCurve<dfm2::CVec2d>(
        scr,
        p0,p1,p2,30,3);
      const dfm2::CVec2d v1 = dfm2::PointOnQuadraticBezierCurve(t1,p0,p1,p2);
      EXPECT_NEAR( (v0-scr).norm(), (v1-scr).norm(),  0.02 );
      //
      const double t2 = dfm2::Nearest_QuadraticBezierCurve_Sturm<dfm2::CVec2d>(
        scr,
        p0,p1,p2,16);
      const dfm2::CVec2d v2 = dfm2::PointOnQuadraticBezierCurve(t2,p0,p1,p2);
      EXPECT_NEAR( (v0-scr).norm(), (v2-scr).norm(),  1.0e-6 );
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


