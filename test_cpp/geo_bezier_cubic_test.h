//
// Created by Nobuyuki Umetani on 2022/01/15.
//

#ifndef GEO_BEZIER_CUBIC_TEST_H_
#define GEO_BEZIER_CUBIC_TEST_H_

#include <random>

#include "gtest/gtest.h"

#include "delfem2/geo_curve_cubic.h"
#include "delfem2/geo_polyline.h"

template <typename VEC>
void TestBezierCubicUsingPolyline(unsigned int nitr)
{
  namespace dfm2 = delfem2;
  std::mt19937 rndeng(std::random_device{}());
  std::uniform_real_distribution<double> dist_01(0, 1);
  for (int itr = 0; itr < nitr; ++itr) {
    VEC p0, p1, p2, p3;
    p0 = VEC(dist_01(rndeng), dist_01(rndeng));
    p1 = VEC(dist_01(rndeng), dist_01(rndeng));
    p2 = VEC(dist_01(rndeng), dist_01(rndeng));
    p3 = VEC(dist_01(rndeng), dist_01(rndeng));
    std::vector<VEC> polyline;
    dfm2::Polyline_BezierCubic(polyline, 1000, p0, p1, p2, p3);
    {
      double v0 = dfm2::Length_Polyline<VEC>(polyline);
      double v2 = dfm2::Length_CubicBezierCurve_Quadrature<VEC>(p0, p1, p2, p3, 5);
      EXPECT_NEAR(v0, v2, 0.07);
      double v3 = dfm2::Length_CubicBezierCurve_QuadratureSubdivision<VEC>(
          p0, p1, p2, p3, 1.0e-3, 12);
      EXPECT_NEAR(v0, v3, 0.008);
    }
    {
      double a0 = dfm2::Area_CubicBezierCurve2(p0, p1, p2, p3);
      double a1 = dfm2::Area_Polyline2(polyline);
      EXPECT_NEAR(a0, a1, 1.0e-5);
    }
    {  // nearest
      const VEC scr(0.5, 0.5);
      const double param = dfm2::Nearest_Polyline(polyline, scr);
      const VEC v0 = dfm2::Sample_Polyline(polyline, param);
      const double t1 = dfm2::Nearest_CubicBezierCurve<VEC>(
          scr,
          p0, p1, p2, p3, 50, 3);
      const VEC v1 = dfm2::PointOnCubicBezierCurve(t1, p0, p1, p2, p3);
      EXPECT_NEAR((v0 - scr).norm(), (v1 - scr).norm(), 1.6e-2);
      //
      const double t2 = dfm2::Nearest_CubicBezierCurve_Strum<VEC>(
          scr,
          p0, p1, p2, p3, 16);
      const VEC v2 = dfm2::PointOnCubicBezierCurve(t2, p0, p1, p2, p3);
      EXPECT_NEAR((v0 - scr).norm(), (v2 - scr).norm(), 2.0e-6);
    }
    {  // bb
      auto bb0 = dfm2::AABB_CubicBezierCurve(p0, p1, p2, p3);
      auto bb1 = dfm2::AABB_Polyline(polyline);
      EXPECT_NEAR(bb0[0], bb1[0], 1.0e-3);
      EXPECT_NEAR(bb0[1], bb1[1], 1.0e-3);
      EXPECT_NEAR(bb0[2], bb1[2], 1.0e-3);
      EXPECT_NEAR(bb0[3], bb1[3], 1.0e-3);
    }
    {  // bezier split
      std::vector<VEC> poly0;
      {
        std::array<VEC, 8> r = dfm2::Split_CubicBezierCurve(p0, p1, p2, p3, 0.3);
        delfem2::Polyline_BezierCubic(
            poly0, 4,
            r[0], r[1], r[2], r[3]);
        std::vector<VEC> poly1;
        delfem2::Polyline_BezierCubic(
            poly1, 8,
            r[4], r[5], r[6], r[7]);
        poly0.resize(poly0.size() - 1);
        poly0.insert(poly0.end(), poly1.begin(), poly1.end());
      }
      std::vector<VEC> poly2;
      delfem2::Polyline_BezierCubic(
          poly2, 11,
          p0, p1, p2, p3);
      EXPECT_EQ(poly0.size(), poly2.size());
      for (unsigned int ip = 0; ip < poly0.size(); ++ip) {
        EXPECT_LT((poly0[ip] - poly2[ip]).norm(), 1.0e-6);
      }
    }
  }
}

#endif //GEO_BEZIER_CUBIC_TEST_H_
