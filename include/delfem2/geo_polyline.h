/*
 * Copyright (c) 2019 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */

/**
 * @detail The order of dependency in delfem2: \n
 * aabb -> \n
 * line -> ray -> edge -> polyline -> \n
 * curve_quadratic -> curve_cubic -> curve_ndegree -> \n
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

template<class VEC>
std::vector<VEC> Polyline_Resample_Polyline(
    const std::vector<VEC> &stroke0,
    typename VEC::Scalar l);

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
VEC Sample_Polyline(
    const std::vector<VEC> &polyline,
    typename VEC::Scalar param) {
  using SCALAR = typename VEC::Scalar;
  if (polyline.empty()) { return VEC(0, 0); }
  if (param < 0) { return polyline[0]; }
  if (param > static_cast<SCALAR>(polyline.size())) {
    return polyline[polyline.size() - 1];
  }
  const unsigned int ie = int(param);
  const SCALAR ratio = param - static_cast<SCALAR>(ie);
  assert(ratio >= 0 && ratio <= 1);
  return (1 - ratio) * polyline[ie] + ratio * polyline[ie + 1];
}

template<typename VEC>
VEC Tangent_Polyline(
    const std::vector<VEC> &polyline,
    unsigned int ie,
    [[maybe_unused]] float ratio) {
  assert(ie < polyline.size() - 1);
  return (polyline[ie + 1] - polyline[ie]).normalized();
}

/**
 *
 * @tparam VEC delfem2::CVecX or Eigen::VectorX
 * @param polyline
 * @param scr
 * @return
 */
template<typename VEC>
[[nodiscard]] typename VEC::Scalar Nearest_Polyline(
    const std::vector<VEC> &polyline,
    const VEC &scr) {
  using SCALAR = typename VEC::Scalar;
  if (polyline.empty()) { return -1; }
  if (polyline.size() == 1) { return 0; }
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
  return static_cast<SCALAR>(ie_min) + ratio_min;
}

template<typename VEC>
void Smooth_Polyline(
    std::vector<VEC> &xys,
    unsigned int ivtx,
    typename VEC::Scalar damping) {
  const int ixy0 = static_cast<int>(ivtx);
  for (int ixy = ixy0 - 2; ixy < ixy0 + 2; ++ixy) {
    if (ixy - 1 < 0 || ixy + 1 >= static_cast<int>(xys.size())) { continue; }
    xys[ixy] = xys[ixy] * damping + (1 - damping) * (xys[ixy - 1] * 0.5 + xys[ixy + 1] * 0.5);
  }
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


