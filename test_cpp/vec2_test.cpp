

#include "gtest/gtest.h" // need to be defiend in the beginning

#include "delfem2/geo_polygon2.h"
#include "delfem2/geo_polyline2.h"
#include "delfem2/vec2.h"
#include "delfem2/quat.h"

namespace dfm2 = delfem2;

TEST(vec2, second_moment_of_area) {
  std::random_device randomDevice;
  std::mt19937 rdeng(randomDevice());
  std::uniform_real_distribution<double> dist01(0.0, 1.0);
  for (int itr = 0; itr < 100; itr++) {
    const double r0 = dist01(rdeng);
    const double r1 = dist01(rdeng);
    const double r2 = dist01(rdeng);
    const double r3 = dist01(rdeng);
    const double r4 = dist01(rdeng);
    const double a = 10 * r0;
    const double b = a * (3 * r1 + 1);
    std::vector<dfm2::CVec2d> aVec2;
    {
      aVec2.emplace_back(-a * 0.5, -b * 0.5);
      aVec2.emplace_back(+a * 0.5, -b * 0.5);
      aVec2.emplace_back(+a * 0.5, +b * 0.5);
      aVec2.emplace_back(-a * 0.5, +b * 0.5);
      double theta0 = r4 * 3.1415 * 2.0;
      dfm2::Rotate(aVec2, theta0);
      dfm2::Translate(aVec2, r2 * 10 - 5, r3 * 10 - 5);
    }
    dfm2::CVec2d cg, pa1, pa2;
    double area, I1, I2;
    dfm2::SecondMomentOfArea_Polygon(cg, area, pa1, I1, pa2, I2,
                                     aVec2);
    EXPECT_NEAR(area, a * b, 1.0e-10);
    EXPECT_NEAR(pa1.dot(pa2), 0.0, 1.0e-10);
    EXPECT_TRUE(I1 >= I2);
//    EXPECT_NEAR(pa1.x*pa1.x,  1.0,          1.0e-10);
//    EXPECT_NEAR(pa1.y,        0.0,          1.0e-10);
    EXPECT_NEAR(I1, a * b * b * b / 12.0, 1.0e-10);
///
//    EXPECT_NEAR(pa2.x,        0.0,          1.0e-10);
//    EXPECT_NEAR(pa2.y*pa2.y,  1.0,          1.0e-10);
    EXPECT_NEAR(I2, b * a * a * a / 12.0, 1.0e-10);
  }
}