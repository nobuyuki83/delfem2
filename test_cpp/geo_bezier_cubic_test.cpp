//
// Created by Nobuyuki Umetani on 2021/12/01.
//

#include "geo_bezier_cubic_test.h"

TEST(bezier_cubic, test0) {
  TestBezierCubicUsingPolyline<delfem2::CVec2d>(1000);
  TestBezierCubicUsingPolyline<delfem2::CVec2f>(1000);
}
