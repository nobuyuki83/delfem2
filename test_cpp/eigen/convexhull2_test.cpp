//
// Created by Nobuyuki Umetani on 2022/01/12.
//

#include <vector>

#include <Eigen/Core>

#include "../convexhull2_test.h"


TEST(convexhull2, test0) {
  ConvexHull2_Test0<Eigen::Vector2d>(1000);
  ConvexHull2_Test0<Eigen::Vector2f>(1000);
}