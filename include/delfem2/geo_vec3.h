/*
* Copyright (c) 2021 Nobuyuki Umetani
*
* This source code is licensed under the MIT license found in the
* LICENSE file in the root directory of this source tree.
*/

/**
 * geometry of 3d vector
 * @detail The order of dependency in delfem2 is "vec2 -> mat2 -> vec3 -> quaternion -> mat3 -> mat4",
 */

#ifndef DFM2_GEO_VEC3_H
#define DFM2_GEO_VEC3_H

#include <cassert>
#include <cmath>
#include <iostream>
#include <vector>
#include <array>

// should not depends on "delfem2/vec3.h"
#include "delfem2/dfm2_inline.h"
#include "delfem2/geo_meta_funcs.h"

namespace delfem2 {

/**
 *
 * @tparam VEC dfm2::CVec3, Eigen::Vector3, std::array<*,3>, *[3]
 * @tparam T
 * @return
 */
template<typename VEC, typename T = value_type<VEC>>
T Dot3(
    const VEC &a,
    const VEC &b) {
  return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}

// --------------------------------------

template<typename VEC0, typename VEC1, typename T = value_type<VEC0>>
T Distance3(
    const VEC0 &p0,
    const VEC1 &p1) {
  const T v0 = p1[0] - p0[0];
  const T v1 = p1[1] - p0[1];
  const T v2 = p1[2] - p0[2];
  return std::sqrt(v0 * v0 + v1 * v1 + v2 * v2);
}

// ---------------------------

template<typename VEC, typename T = value_type<VEC>>
 T Length3(
    const VEC &v) {
  return std::sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);
}

// ---------------------------

template<typename VEC, typename T = value_type<VEC>>
void Normalize3(VEC &&v) {
  const T len = std::sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);
  v[0] /= len;
  v[1] /= len;
  v[2] /= len;
}

// ---------------------------

template<typename VEC, typename T = value_type<VEC>>
T SquareDistance3(
    const VEC& p0,
    const VEC& p1) {
  const T v0 = p1[0] - p0[0];
  const T v1 = p1[1] - p0[1];
  const T v2 = p1[2] - p0[2];
  return v0 * v0 + v1 * v1 + v2 * v2;
}

// ---------------------------

template<typename T>
T SquareLength3(
    const T v[3]) {
  return v[0] * v[0] + v[1] * v[1] + v[2] * v[2];
}

// -------------------

template<typename VEC0, typename VEC1, typename VEC2>
void Cross(
    VEC0 &&lhs,
    const VEC1 &v1,
    const VEC2 &v2) {
  lhs[0] = v1[1] * v2[2] - v2[1] * v1[2];
  lhs[1] = v1[2] * v2[0] - v2[2] * v1[0];
  lhs[2] = v1[0] * v2[1] - v2[0] * v1[1];
}

// ------------------------

template<typename VEC0>
VEC0 Cross(
    const VEC0 &v1,
    const VEC0 &v2){
  return {
      v1[1] * v2[2] - v2[1] * v1[2],
      v1[2] * v2[0] - v2[2] * v1[0],
      v1[0] * v2[1] - v2[0] * v1[1]
  };
}

// ------------------------

template<typename VEC, typename T = value_type<VEC>>
T ScalarTripleProduct(
    const VEC &a,
    const VEC &b,
    const VEC &c) {
  const T v0 = a[0] * (b[1] * c[2] - b[2] * c[1]);
  const T v1 = a[1] * (b[2] * c[0] - b[0] * c[2]);
  const T v2 = a[2] * (b[0] * c[1] - b[1] * c[0]);
  return v0 + v1 + v2;
}

// -------------------------

template<typename VEC0, typename VEC1>
void FrameFromVectorZ(
    VEC0 &&vec_x,
    VEC0 &&vec_y,
    const VEC1 &vec_n);


// =====================================
// dependency to the std::vector

template<typename VEC>
std::ostream &operator<<(
    std::ostream &output,
    const std::vector<VEC> &aV);

// -------------------------

template<typename VEC>
std::istream &operator>>(
    std::istream &input,
    std::vector<VEC> &aV);

}

#ifndef DFM2_STATIC_LIBRARY
#  include "delfem2/geo_vec3.cpp"
#endif

#endif // DFM2_GEO_VEC3_H
