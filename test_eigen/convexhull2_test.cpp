//
// Created by Nobuyuki Umetani on 2022/01/12.
//

#include <gtest/gtest.h>

#include <random>
#include <set>
#include <vector>

#include <Eigen/Core>

#include "delfem2/geo_convhull2.h"
#include "delfem2/geo_polyline.h"

TEST(convexhull2, test0) {
  std::mt19937 rndeng(std::random_device{}());
  std::uniform_real_distribution<> dist_m2p2(-2, 2);

  using VEC = Eigen::Vector2d;
  using SCALAR = typename VEC::Scalar;

  for(int ii = 0; ii < 100; ii++) {

    std::vector<VEC> points;
    for(int i = 0; i < 10; i++) {
      points.emplace_back(dist_m2p2(rndeng), dist_m2p2(rndeng));
    }

    std::vector<unsigned int> polygon_pointidx;
    delfem2::ConvexHull2<VEC>(polygon_pointidx, points);
    ASSERT_GT(polygon_pointidx.size(), 2);

    // checks whether all angles are less than Pi
    const unsigned int num_points_polygon = polygon_pointidx.size();
    for(unsigned int i = 0; i < num_points_polygon; i++) {
      const unsigned int last = (i + num_points_polygon - 1) % num_points_polygon;
      const unsigned int next = (i + 1) % num_points_polygon;
      const VEC p1p2 = points[i] - points[last];
      const VEC p2p3 = points[next] - points[i];
      const SCALAR cos_val = p1p2.dot(p2p3) / p1p2.norm() / p2p3.norm();
      EXPECT_GE(cos_val, -1);
    }

    // checks winding number for internal points are 2 * Pi
    std::vector<VEC> polygon;
    std::set<unsigned int> boundary_point_idx;
    for(unsigned int i = 0; i < num_points_polygon; i++) {
      polygon.push_back(points[polygon_pointidx[i]]);
      boundary_point_idx.emplace(polygon_pointidx[i]);
    }
    for(unsigned int i = 0; i < points.size(); i++) {
      if(boundary_point_idx.count(i) != 0) { continue; }
      const SCALAR wn0 = delfem2::WindingNumber_Polyline2<VEC>(polygon, points[i]);
      EXPECT_NEAR(wn0, M_PI * 2, 0.0001);
    }
  }
}