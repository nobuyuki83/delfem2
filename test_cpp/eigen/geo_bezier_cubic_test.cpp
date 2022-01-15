//
// Created by Nobuyuki Umetani on 2021/12/01.
//

#include <Eigen/Core>

#include "../geo_bezier_cubic_test.h"

TEST(bezier_cubic, test0) {
  TestBezierCubicUsingPolyline<Eigen::Vector2d>(1000);
  TestBezierCubicUsingPolyline<Eigen::Vector2f>(1000);
}
