/*
 * Copyright (c) 2019 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */

/**
 * @detail The order of dependency in delfem2:
 * aabb ->
 * line -> ray -> edge -> polyline ->
 * curve_quadratic -> curve_cubic -> curve_ndegree ->
 * plane -> tri -> quad
 */

#ifndef DFM2_GEO_POLYLINE2_H
#define DFM2_GEO_POLYLINE2_H

#include <cstdlib>
#include <cassert>
#include <iostream>
#include <cmath>
#include <vector>
#include <climits>
#include <algorithm>  // for std::clamp

#include "delfem2/vec2.h"
#include "delfem2/geo_curve_cubic.h"
#include "delfem2/dfm2_inline.h"

// -----------------------------------------------------

namespace delfem2 {

//! @brief translate all the points
template<typename T>
void Translate(
  std::vector<CVec2<T> > &aP,
  double dx,
  double dy);

template<typename T>
DFM2_INLINE void Rotate(
  std::vector<CVec2<T> > &aP,
  double dt);

template<typename T>
void Polyline_CubicBezierCurve(
  std::vector<CVec2<T> > &aP,
  int n,
  const std::vector<CVec2<T> > &aCP);

template<class VEC>
void Polyline_BezierCubic(
  std::vector<VEC> &aP,
  unsigned int n,
  const VEC &p1,
  const VEC &p2,
  const VEC &p3,
  const VEC &p4);

template<class VEC>
std::vector<VEC> Polyline_Resample_Polyline(
  const std::vector<VEC> &stroke0,
  double l);

/**
 *
 * @tparam VEC delfem2::CVecX or Eigen::VectorX
 * @param polyline
 * @param scr
 * @return
 */
template<typename VEC>
unsigned int FindNearestPointInPoints(
  const std::vector<VEC> &polyline,
  const VEC &scr) {
  unsigned int idx_point_min_dist = UINT_MAX;
  float dist_min = -1;
  for (unsigned int ip = 0; ip < polyline.size(); ++ip) {
    auto dist = (scr - polyline[ip]).norm();
    if (dist_min < 0 || dist < dist_min) {
      idx_point_min_dist = ip;
      dist_min = dist;
    }
  }
  return idx_point_min_dist;
}

/**
 *
 * @tparam VEC delfem2::CVecX or Eigen::VectorX
 * @param polyline
 * @param scr
 * @return
 */
template<class VEC>
typename VEC::Scalar Length_Polyline(
  const std::vector<VEC> &polyline) {
  if (polyline.size() < 2) {
    return 0;
  }
  using SCALAR = typename VEC::Scalar;
  SCALAR len = 0;
  for (unsigned int ip = 0; ip < polyline.size() - 1; ++ip) {
    len += (polyline[ip + 1] - polyline[ip]).norm();
  }
  return len;
}

template<typename VEC>
auto Area_Polyline2(
  const std::vector<VEC> &polyline) -> typename VEC::Scalar {
  if (polyline.size() < 2) {
    return 0;
  }
  using SCALAR = typename VEC::Scalar;
  SCALAR area = 0;
  VEC p0(0, 0);
  for (unsigned int ip = 0; ip < polyline.size() - 1; ++ip) {
    area += Area_Tri2(p0, polyline[ip], polyline[ip + 1]);
  }
  return area;
}

template<typename VEC>
auto AABB_Polyline(
  const std::vector<VEC> &polyline) -> std::array<typename VEC::Scalar, VEC::SizeAtCompileTime * 2> {
  using SCALAR = typename VEC::Scalar;
  constexpr int ndim = VEC::SizeAtCompileTime;
  std::array<SCALAR, ndim * 2> res;
  for (int idim = 0; idim < ndim; ++idim) {
    res[idim] = res[idim + ndim] = polyline[0][idim];
  }
  for (int ip = 1; ip < polyline.size(); ++ip) {
    for (int idim = 0; idim < ndim; ++idim) {
      res[idim] = (polyline[ip][idim] < res[idim]) ? polyline[ip][idim] : res[idim];
      res[idim + ndim] = (polyline[ip][idim] > res[idim + ndim]) ? polyline[ip][idim] : res[idim + ndim];
    }
  }
  return res;
}

template<typename VEC>
VEC PositionInPolyline(
  const std::vector<VEC> &polyline,
  unsigned int ie,
  typename VEC::Scalar ratio) {
  assert(ie < polyline.size() - 1);
  return (1 - ratio) * polyline[ie] + ratio * polyline[ie + 1];
}

template<typename VEC>
VEC NormalInPolyline(
  const std::vector<VEC> &polyline,
  unsigned int ie,
  [[maybe_unused]] float ratio) {
  assert(ie < polyline.size() - 1);
  VEC ut = (polyline[ie + 1] - polyline[ie]).normalized();
  return rotate90(ut);
}

/**
 *
 * @tparam VEC delfem2::CVecX or Eigen::VectorX
 * @param polyline
 * @param scr
 * @return
 */
template<typename VEC>
[[nodiscard]] auto FindNearestPointInPolyline(
  const std::vector<VEC> &polyline,
  const VEC &scr) -> std::pair<unsigned int, typename VEC::Scalar> {
  using SCALAR = typename VEC::Scalar;
  if (polyline.empty()) { return {UINT_MAX, 0}; }
  unsigned int ie_min;
  SCALAR ratio_min;
  SCALAR dist_min = -1;
  for (unsigned int ip = 0; ip < polyline.size() - 1; ++ip) {
    const VEC &es = polyline[ip + 1] - polyline[ip];
    const VEC &sc = polyline[ip] - scr;
    const SCALAR a = es.squaredNorm();
    const SCALAR b = es.dot(sc);
    const SCALAR ratio = std::clamp(-b / a, static_cast<SCALAR>(0), static_cast<SCALAR>(1));
    VEC p = (1 - ratio) * polyline[ip] + ratio * polyline[ip + 1];
    SCALAR dist = (p - scr).norm();
    if (dist_min < 0 || dist < dist_min) {
      dist_min = dist;
      ie_min = ip;
      ratio_min = ratio;
    }
  }
  return {ie_min, ratio_min};
}

/**
 *
 * @tparam VEC delfem2::CVecX or Eigen::VectorX
 * @param polyline
 * @param scr
 * @return
 */
template<typename VEC>
auto ArcLengthPointInPolyline(
  const std::vector<VEC> &polyline,
  const VEC &scr) -> typename VEC::Scalar {
  using SCALAR = typename VEC::Scalar;
  if (polyline.size() < 2) { return 0; }
  float dist_min = -1;
  VEC p_min;
  unsigned int ip_min = -1;
  for (unsigned int ip = 0; ip < polyline.size() - 1; ++ip) {
    VEC p_near = Nearest_Edge_Point(scr, polyline[ip], polyline[ip + 1]);
    float dist = (p_near - scr).norm();
    if (dist_min < 0 || dist < dist_min) {
      dist_min = dist;
      p_min = p_near;
      ip_min = ip;
    }
  }
  SCALAR alen = 0;
  for (unsigned int ip = 0; ip < ip_min; ++ip) {
    alen += (polyline[ip + 1] - polyline[ip]).norm();
  }
  alen += (p_min - polyline[ip_min]).norm();
  return alen;
}

} // namespace delfem2

#ifndef DFM2_STATIC_LIBRARY
#  include "delfem2/geo_polyline.cpp"
#endif

#endif // DFM2_GEO_POLYLINE2_H


