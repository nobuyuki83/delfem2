/*
* Copyright (c) 2019 Nobuyuki Umetani
*
* This source code is licensed under the MIT license found in the
* LICENSE file in the root directory of this source tree.
*/

/**
 * @file convex hull computation in 3D space
 */

#ifndef DFM2_GEO_CONVHULL2_H
#define DFM2_GEO_CONVHULL2_H

namespace delfem2 {

/**
 * This function return the array of indices of the input 2D points on the convex hull
 * in the counterclockwise orientation.
 * @tparam VEC CVec2X or Eigen::Vector2x
 * @param point_idxs
 * @param points
 */
template<typename VEC>
void ConvexHull2(
    std::vector<unsigned int> &point_idxs,
    const std::vector<VEC> &points) {
  using SCALAR = typename VEC::Scalar;
  constexpr SCALAR EPSILON = 1.0e-8;

  unsigned int p0_idx;
  { // find the index with minimum y coordinate
    VEC p0;
    for (unsigned int i = 0; i < points.size(); ++i) {
      if (points[i][1] > p0[1]) { continue; }
      if (points[i][1] == p0[1] && points[i][0] > p0[0]) { continue; }
      p0_idx = i;
      p0 = points[i];
    }
  }

  std::vector<std::pair<unsigned int, SCALAR> > idxcos;
  { // compute and sort points by cosine value
    VEC x_axis = {1, 0};
    for (unsigned int i = 0; i < points.size(); i++) {
      if (i == p0_idx) { continue; }
      VEC dir = points[i] - points[p0_idx];
      idxcos.emplace_back(i, x_axis.dot(dir) / dir.norm());
    }
  }
  { // sort idxcos
    auto comp = [&points, &p0_idx](
        const std::pair<unsigned int, SCALAR> &a,
        const std::pair<unsigned int, SCALAR> &b) {
      if (std::abs(a.second - b.second) > EPSILON) {
        return a.second > b.second;
      } else {
        const SCALAR dista = (points[a.first] - points[p0_idx]).squaredNorm();
        const SCALAR distb = (points[b.first] - points[p0_idx]).squaredNorm();
        return dista > distb;
      }
    };
    std::sort(idxcos.begin(), idxcos.end(), comp);
  }

  // check for collinear points
  for (auto itr = ++idxcos.begin(); std::next(itr) != idxcos.end(); itr++) {
    if (std::abs((*itr).second - (*std::next(itr)).second) < EPSILON) {
      idxcos.erase(std::next(itr)); // only keep the furthest
    }
  }
  idxcos.emplace_back(p0_idx, 0);

  point_idxs.clear();
  point_idxs.push_back(p0_idx);
  point_idxs.push_back(idxcos.begin()->first);
  unsigned int stack_top = 1;
  for (auto itr = ++idxcos.begin(); itr != idxcos.end(); itr++) {
    unsigned int p3_idx = itr->first;
    while (true) {
      unsigned int p1_idx = point_idxs[stack_top - 1];
      unsigned int p2_idx = point_idxs[stack_top];
      const VEC p1p2 = points[p2_idx] - points[p1_idx];
      const VEC p1p3 = points[p3_idx] - points[p1_idx];
      if (p1p2[0] * p1p3[1] - p1p2[1] * p1p3[0] <= 0) { // right turn or collinear
        point_idxs.pop_back(); // pop top of the stack
        stack_top--;
      } else {
        break;
      }
    }
    point_idxs.push_back(p3_idx);
    stack_top++;
  }
}

}

#endif // DFM2_GEO_CONVHULL2_H
