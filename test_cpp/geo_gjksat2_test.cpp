//
// Created by Nobuyuki Umetani on 2022/01/18.
//

#include <random>

#include "gtest/gtest.h"

#include "delfem2/vec2.h"
#include "geo_gjksat2_test.h"

TEST(gjksat2, test0) {
  TestGjkSat2Test0<delfem2::CVec2d>(1000);
}