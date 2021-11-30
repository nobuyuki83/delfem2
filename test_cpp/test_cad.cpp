/*
 * Copyright (c) 2019 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */

#include <string>

#include "gtest/gtest.h"

#include "delfem2/cad2.h"
#include "delfem2/cad2_io_svg.h"
#include "delfem2/geo_curve_cubic_bezier.h"
#include "delfem2/geo_curve_quadratic_bezier.h"
#include "delfem2/geo_polyline2.h"

namespace dfm2 = delfem2;

TEST(cad,read_svg) {
  {
    std::string path_svg = std::string(PATH_INPUT_DIR) + "/shape0.svg";
    dfm2::CCad2D cad;
    dfm2::ReadSVG_Cad2D(cad, path_svg, 1.0);
  }
  {
    std::string path_svg = std::string(PATH_INPUT_DIR) + "/shape1.svg";
    dfm2::CCad2D cad;
    dfm2::ReadSVG_Cad2D(cad, path_svg, 1.0);
  }
  {
    std::string path_svg = std::string(PATH_INPUT_DIR) + "/shape2.svg";
    dfm2::CCad2D cad;
    dfm2::ReadSVG_Cad2D(cad, path_svg, 1.0);
  }
  {
    std::string path_svg = std::string(PATH_INPUT_DIR) + "/tshirt.svg";
    dfm2::CCad2D cad;
    dfm2::ReadSVG_Cad2D(cad, path_svg, 1.0);
  }
  {
    std::string path_svg = std::string(PATH_INPUT_DIR) + "/ltshirt.svg";
    dfm2::CCad2D cad;
    dfm2::ReadSVG_Cad2D(cad, path_svg, 1.0);
  }
  {
    std::string path_svg = std::string(PATH_INPUT_DIR) + "/raglan.svg";
    dfm2::CCad2D cad;
    dfm2::ReadSVG_Cad2D(cad, path_svg, 1.0);
  }
  {
    std::string path_svg = std::string(PATH_INPUT_DIR) + "/raglan2.svg";
    dfm2::CCad2D cad;
    dfm2::ReadSVG_Cad2D(cad, path_svg, 1.0);
  }
}


TEST(cad,quadratic_bezier0) {
  for(unsigned int itr=0;itr<2;++itr){
    dfm2::CVec2d p0(0, 0), p1(1, 1), p2(2, 0);
    if( itr == 1 ) {
      p0 = dfm2::CVec2d(0, 0);
      p1 = dfm2::CVec2d(1, 1);
      p2 = dfm2::CVec2d(2, 2);
    }
    std::vector<dfm2::CVec2d> polyline;
    dfm2::Polyline_BezierQuadratic(polyline, 100, p0, p1, p2);
    {  // length
      double v1 = dfm2::Length_QuadraticBezierCurve_Analytic<dfm2::CVec2d>(
        p0, p1, p2);
      double v0 = dfm2::Length_Polyline<dfm2::CVec2d>(polyline);
      double v2 = dfm2::Length_QuadraticBezierCurve_Quadrature<dfm2::CVec2d>(
        p0,p1,p2,3);
      EXPECT_NEAR(v1,v0,2.0e-4);
      EXPECT_NEAR(v1,v2,2.0e-4);
    }
    {  // area
      double a0 = dfm2::Area_QuadraticBezierCurve2(p0, p1, p2);
      double a1 = dfm2::Area_Polyline2(polyline);
      EXPECT_NEAR(a1,a0,2.0e-4);
    }
    {  // nearest
      dfm2::CVec2d scr(0.8, 1.2);
      const auto [i0,l0] = dfm2::FindNearestPointInPolyline(polyline, scr);
      const dfm2::CVec2d v0 = dfm2::PositionInPolyline(polyline, i0, l0);
      const double t1 = dfm2::Nearest_QuadraticBezierCurve<dfm2::CVec2d>(
        scr,
        p0,p1,p2,10,5);
      const dfm2::CVec2d v1 = dfm2::PointOnQuadraticBezierCurve(t1,p0,p1,p2);
      EXPECT_LT( (v0-v1).norm(), 1.0e-2 );
    }
    {  // bb
      auto bb = dfm2::AABB_QuadraticBezierCurve<2>(p0,p1,p2);
    }
  }
}


TEST(cad,cubic_bezier0) {
  dfm2::CVec2d p0(0, 1), p1(1, 2), p2(2, 1), p3(3,2);
  std::vector<dfm2::CVec2d> polyline;
  dfm2::Polyline_BezierCubic(polyline, 100, p0, p1, p2, p3);
  {
    double v0 = dfm2::Length_Polyline<dfm2::CVec2d>(polyline);
    double v2 = dfm2::Length_CubicBezierCurve_Quadrature<dfm2::CVec2d>(p0,p1,p2,p3,3);
    EXPECT_NEAR(v0,v2,1.0e-3);
  }
  {
    double a0 = dfm2::Area_CubicBezierCurve2(p0, p1, p2, p3);
    double a1 = dfm2::Area_Polyline2(polyline);
    EXPECT_NEAR(a0, a1, 1.0e-5);
  }
  {  // nearest
    dfm2::CVec2d scr(0.8, 1.2);
    const auto [i0,l0] = dfm2::FindNearestPointInPolyline(polyline, scr);
    const dfm2::CVec2d v0 = dfm2::PositionInPolyline(polyline, i0, l0);
    const double t1 = dfm2::Nearest_CubicBezierCurve<dfm2::CVec2d>(
      scr,
      p0,p1,p2,p3,10,5);
    const dfm2::CVec2d v1 = dfm2::PointOnCubicBezierCurve(t1,p0,p1,p2,p3);
    EXPECT_LT( (v0-v1).norm(), 1.0e-3 );
  }
  {  // bb
      auto bb = dfm2::AABB_CubicBezierCurve<2>(p0,p1,p2,p3);
  }
}
