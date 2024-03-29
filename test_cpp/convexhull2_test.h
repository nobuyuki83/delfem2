//
// Created by Nobuyuki Umetani on 2022/01/15.
//

#ifndef CONVEXHULL2_TEST_H_
#define CONVEXHULL2_TEST_H_

#include <random>
#include <set>

#include <gtest/gtest.h>

#include "delfem2/geo_convhull2.h"
#include "delfem2/geo_polyline.h"
#include "delfem2/vec2.h"

#ifndef M_PI
#  define M_PI 3.14159265359
#endif

template<typename VEC, typename SCALAR = typename VEC::Scalar>
void ConvexHull2_Test0(unsigned int nitr) {

  const unsigned int N = 10;
  std::mt19937 rndeng(std::random_device{}());
  std::uniform_real_distribution<SCALAR> dist_m2p2(-2, 2);

  for (int ii = 0; ii < nitr; ii++) {
    std::vector<VEC> points(N);
    for (int i = 0; i < N; i++) {
      points[i] = {dist_m2p2(rndeng), dist_m2p2(rndeng)};
    }

    std::vector<unsigned int> polygon_pointidx;
    delfem2::ConvexHull2<VEC>(polygon_pointidx, points);
    ASSERT_GT(polygon_pointidx.size(), 2);

    // checks whether all angles are less than Pi
    const unsigned int num_points_polygon = polygon_pointidx.size();
    for (unsigned int ip1 = 0; ip1 < num_points_polygon; ip1++) {
      const unsigned int ip0 = (ip1 + num_points_polygon - 1) % num_points_polygon;
      const unsigned int ip2 = (ip1 + 1) % num_points_polygon;
      EXPECT_NE(ip0, ip1);
      EXPECT_NE(ip1, ip2);
      const VEC p1p2 = points[polygon_pointidx[ip1]] - points[polygon_pointidx[ip0]];
      const VEC p2p3 = points[polygon_pointidx[ip2]] - points[polygon_pointidx[ip1]];
      const SCALAR sin_val = delfem2::Cross2(p1p2, p2p3);
      const SCALAR v0 = p1p2.norm() * p2p3.norm();
      EXPECT_GE(sin_val, 1.0e-10 * v0);
    }

    // checks winding number for internal points are 2 * Pi
    std::vector<VEC> polygon;
    std::set<unsigned int> boundary_point_idx;
    for (unsigned int i = 0; i < num_points_polygon; i++) {
      polygon.push_back(points[polygon_pointidx[i]]);
      boundary_point_idx.emplace(polygon_pointidx[i]);
    }
    for (unsigned int i = 0; i < points.size(); i++) {
      if (boundary_point_idx.count(i) != 0) { continue; }
      const SCALAR wn0 = delfem2::WindingNumber_Polyline2<VEC>(polygon, points[i]);
      EXPECT_NEAR(wn0, M_PI * 2, 0.0001);
    }
  }
};

#endif //CONVEXHULL2_TEST_H_
