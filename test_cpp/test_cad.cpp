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
#include "delfem2/pgeo.h"
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
  /*
  if( iframe == nframe_interval*1 ){ path_svg = std::string(PATH_INPUT_DIR)+"/shape1.svg"; }
  if( iframe == nframe_interval*2 ){ path_svg = std::string(PATH_INPUT_DIR)+"/shape2.svg"; }
  if( iframe == nframe_interval*3 ){ path_svg = std::string(PATH_INPUT_DIR)+"/tshirt.svg"; }
  if( iframe == nframe_interval*4 ){ path_svg = std::string(PATH_INPUT_DIR)+"/ltshirt.svg"; }
  if( iframe == nframe_interval*5 ){ path_svg = std::string(PATH_INPUT_DIR)+"/raglan.svg"; }
  if( iframe == nframe_interval*6 ){ path_svg = std::string(PATH_INPUT_DIR)+"/raglan2.svg"; }
   */
}


TEST(cad,length_quadratic_bezier) {
  {
    dfm2::CVec2d p0(0, 0), p1(1, 1), p2(2, 0);
    std::vector<dfm2::CVec2d> aP;
    dfm2::Polyline_BezierQuadratic(aP, 100, p0, p1, p2);
    double v1 = dfm2::Length_QuadraticBezierCurve_Analytic<dfm2::CVec2d, double>(p0, p1, p2);
    double v0 = dfm2::LengthPolyline<dfm2::CVec2d,double>(aP);
    double v2 = dfm2::Length_QuadraticBezierCurve_Quadrature<dfm2::CVec2d, double>(p0,p1,p2,3);
    EXPECT_NEAR(v1,v0,1.0e-3);
    EXPECT_NEAR(v1,v2,1.0e-3);
  }
  {
    dfm2::CVec2d p0(0, 0), p1(1, 1), p2(2, 2);
    std::vector<dfm2::CVec2d> aP;
    dfm2::Polyline_BezierQuadratic(aP, 100, p0, p1, p2);
    double v1 = dfm2::Length_QuadraticBezierCurve_Analytic<dfm2::CVec2d, double>(p0, p1, p2);
    double v0 = dfm2::LengthPolyline<dfm2::CVec2d,double>(aP);
    double v2 = dfm2::Length_QuadraticBezierCurve_Quadrature<dfm2::CVec2d, double>(p0,p1,p2,3);
    EXPECT_NEAR(v1,v0,1.0e-5);
    EXPECT_NEAR(v1,v2,1.0e-5);
  }
}

TEST(cad,length_cubic_bezier) {
  {
    dfm2::CVec2d p0(0, 0), p1(1, 1), p2(2, 0), p3(3,1);
    std::vector<dfm2::CVec2d> aP;
    dfm2::Polyline_BezierCubic(aP, 100, p0, p1, p2,p3);
    double v0 = dfm2::LengthPolyline<dfm2::CVec2d,double>(aP);
    double v2 = dfm2::Length_CubicBezierCurve_Quadrature<dfm2::CVec2d, double>(p0,p1,p2,p3,3);
    EXPECT_NEAR(v0,v2,1.0e-3);
  }
}