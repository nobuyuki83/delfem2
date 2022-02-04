/*
* Copyright (c) 2019 Nobuyuki Umetani
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

#define NEARLY_ZERO 1.e-16

namespace delfem2 {

template <class T>
class has_definition {
  template <class U>
  static constexpr int check(typename U::Scalar*){ return 0; }

  template <class U>
  static constexpr int check(typename U::value_type*){
    static_assert( std::tuple_size<U>::value == 3, "the size of the std::array needs to be 3");
    return 1; }

  template <class U>
  static constexpr int check(...)  { return 2; }
 public:
  static constexpr int value = check<T>(nullptr);
};

template <int F, class T>
struct conditional
{
  using type = typename T::Scalar;
};

template <class T>
struct conditional<1, T>
{
  using type = typename T::value_type;
};

template <class T>
struct conditional<2, T>
{
  using type = typename std::remove_extent<T>::type;
};

template <class T>
using type_scalar = typename conditional<has_definition<T>::value,T>::type;

/**
 *
 * @tparam VEC dfm2::CVec3, Eigen::Vector3, std::array<*,3>, *[3]
 * @tparam T
 * @return
 */
template<typename VEC, typename T = type_scalar<VEC>>
DFM2_INLINE T Dot3(const VEC &a, const VEC &b) {
  return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}

template<typename T>
DFM2_INLINE T Dot3(const T *a, const T *b) {
  return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}

// --------------------------------------

template<typename VEC, typename T = type_scalar<VEC>>
DFM2_INLINE T Distance3(
  const VEC &p0,
  const VEC &p1) {
  return std::sqrt(
    (p1[0] - p0[0]) * (p1[0] - p0[0]) +
      (p1[1] - p0[1]) * (p1[1] - p0[1]) +
      (p1[2] - p0[2]) * (p1[2] - p0[2]));
}

template<typename T>
DFM2_INLINE T Distance3(
    const T* p0,
    const T* p1) {
  return std::sqrt(
      (p1[0] - p0[0]) * (p1[0] - p0[0]) +
          (p1[1] - p0[1]) * (p1[1] - p0[1]) +
          (p1[2] - p0[2]) * (p1[2] - p0[2]));
}

// ---------------------------

template<typename REAL>
DFM2_INLINE REAL Length3(const REAL v[3]) {
  return std::sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);
}

// ---------------------------

template<typename T>
DFM2_INLINE void Normalize3(T v[3]) {
  T len = std::sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);
  v[0] /= len;
  v[1] /= len;
  v[2] /= len;
}

// ---------------------------

template<typename T>
DFM2_INLINE T SquareDistance3(const T p0[3], const T p1[3]) {
  return
    (p1[0] - p0[0]) * (p1[0] - p0[0]) +
      (p1[1] - p0[1]) * (p1[1] - p0[1]) +
      (p1[2] - p0[2]) * (p1[2] - p0[2]);
}

// ---------------------------

template<typename T>
DFM2_INLINE T SquareLength3(const T v[3]) {
  return v[0] * v[0] + v[1] * v[1] + v[2] * v[2];
}


}

#endif // DFM2_VEC3_H
