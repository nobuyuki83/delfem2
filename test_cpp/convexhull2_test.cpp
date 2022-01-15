//
// Created by Nobuyuki Umetani on 2022/01/12.
//

#include <vector>

#include "delfem2/vec2.h"
#include "convexhull2_test.h"


TEST(convexhull2, test0) {
  ConvexHull2_Test0<delfem2::CVec2d>(1000);
  ConvexHull2_Test0<delfem2::CVec2f>(1000);
}