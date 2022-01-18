//
// Created by Nobuyuki Umetani on 2022/01/18.
//

#include <Eigen/Core>

#include "gtest/gtest.h"

#include "../geo_gjksat2_test.h"

TEST(gjksat2, test0) {
  TestGjkSat2Test0<Eigen::Vector2d>(1000);
}