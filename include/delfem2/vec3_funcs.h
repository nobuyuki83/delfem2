/*
* Copyright (c) 2021 Nobuyuki Umetani
*
* This source code is licensed under the MIT license found in the
* LICENSE file in the root directory of this source tree.
*/

/**
 * funcstions related to 3d vector
 * @detail The order of dependency in delfem2 is "vec2 -> mat2 -> vec3 -> quaternion -> mat3 -> mat4",
 */

#ifndef DFM2_VEC3_FUNCS_H
#define DFM2_VEC3_FUNCS_H

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
template<typename VEC, typename T = vecn_value_t<VEC,3>>
T Dot3(
    const VEC &a,
    const VEC &b) {
  return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}

// --------------------------------------

template<typename VEC0, typename VEC1, typename T = vecn_value_t<VEC0,3>>
T Distance3(
    const VEC0 &p0,
    const VEC1 &p1) {
  const T v0 = p1[0] - p0[0];
  const T v1 = p1[1] - p0[1];
  const T v2 = p1[2] - p0[2];
  return std::sqrt(v0 * v0 + v1 * v1 + v2 * v2);
}

// ---------------------------

template<typename VEC, typename T = vecn_value_t<VEC,3>>
 T Length3(
    const VEC &v) {
  return std::sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);
}

// ---------------------------

template<typename VEC, typename T = vecn_value_t<VEC,3>>
void Normalize3(VEC &&v) {
  const T len = std::sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);
  v[0] /= len;
  v[1] /= len;
  v[2] /= len;
}

// ---------------------------

template<typename VEC, typename T = vecn_value_t<VEC,3>>
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

template<typename VEC, typename T = vecn_value_t<VEC,3>>
T ScalarTripleProduct(
    const VEC &a,
    const VEC &b,
    const VEC &c) {
  const T v0 = a[0] * (b[1] * c[2] - b[2] * c[1]);
  const T v1 = a[1] * (b[2] * c[0] - b[0] * c[2]);
  const T v2 = a[2] * (b[0] * c[1] - b[1] * c[0]);
  return v0 + v1 + v2;
}

// ----------------

template<typename REAL>
void AverageTwo3(
    REAL po[3],
    const REAL p0[3], const REAL p1[3])  {
  constexpr REAL half = static_cast<REAL>(0.5);
  po[0] = (p0[0] + p1[0]) * half;
  po[1] = (p0[1] + p1[1]) * half;
  po[2] = (p0[2] + p1[2]) * half;
}

// --------------------

template<typename REAL>
void AverageFour3(
    REAL po[3],
    const REAL p0[3], const REAL p1[3], const REAL p2[3], const REAL p3[3]) {
  constexpr REAL quarter(0.25);
  po[0] = (p0[0] + p1[0] + p2[0] + p3[0]) * quarter;
  po[1] = (p0[1] + p1[1] + p2[1] + p3[1]) * quarter;
  po[2] = (p0[2] + p1[2] + p2[2] + p3[2]) * quarter;
}

// ------------------

/**
 * @func add values for 3darray (vo += vi)
 * @tparam REAL float and double
 * @param vo (out)
 * @param vi (in)
 */
template<typename REAL>
DFM2_INLINE void Add3(
    REAL vo[3],
    const REAL vi[3]){
  vo[0] += vi[0];
  vo[1] += vi[1];
  vo[2] += vi[2];
}

// -------------------------

template<typename VEC0, typename VEC1>
void FrameFromVectorZ(
    VEC0 &&vec_x,
    VEC0 &&vec_y,
    const VEC1 &vec_n);


template<typename VEC, typename T = vecn_value_t<VEC,3>>
VEC RotateVec3WithAxisAngleVector(
    const VEC &vec0,
    const VEC &rot);

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
#  include "delfem2/vec3_funcs.cpp"
#endif

#endif // DFM2_VEC3_FUNCS_H
