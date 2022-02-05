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
DFM2_INLINE T Dot3(
    const VEC &a,
    const VEC &b) {
  return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}

// --------------------------------------

template<typename VEC0, typename VEC1, typename T = value_type<VEC0>>
DFM2_INLINE T Distance3(
    const VEC0 &p0,
    const VEC1 &p1) {
  const T v0 = p1[0] - p0[0];
  const T v1 = p1[1] - p0[1];
  const T v2 = p1[2] - p0[2];
  return std::sqrt(v0 * v0 + v1 * v1 + v2 * v2);
}

// ---------------------------

template<typename VEC, typename T = value_type<VEC>>
DFM2_INLINE T Length3(
    const VEC &v) {
  return std::sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);
}

// ---------------------------

template<typename VEC, typename T = value_type<VEC>>
DFM2_INLINE void Normalize3(VEC &&v) {
  const T len = std::sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);
  v[0] /= len;
  v[1] /= len;
  v[2] /= len;
}

// ---------------------------

template<typename VEC, typename T = value_type<VEC>>
DFM2_INLINE T SquareDistance3(
    const VEC& p0,
    const VEC& p1) {
  const T v0 = p1[0] - p0[0];
  const T v1 = p1[1] - p0[1];
  const T v2 = p1[2] - p0[2];
  return v0 * v0 + v1 * v1 + v2 * v2;
}

// ---------------------------

template<typename T>
DFM2_INLINE T SquareLength3(
    const T v[3]) {
  return v[0] * v[0] + v[1] * v[1] + v[2] * v[2];
}

// -------------------

template<typename VEC0, typename VEC1, typename VEC2>
DFM2_INLINE void Cross(
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
DFM2_INLINE T ScalarTripleProduct(
    const VEC &a,
    const VEC &b,
    const VEC &c) {
  const T v0 = a[0] * (b[1] * c[2] - b[2] * c[1]);
  const T v1 = a[1] * (b[2] * c[0] - b[0] * c[2]);
  const T v2 = a[2] * (b[0] * c[1] - b[1] * c[0]);
  return v0 + v1 + v2;
}

// -------------------------

template<typename VEC0, typename VEC1, typename T = value_type<VEC0>>
DFM2_INLINE void FrameFromVectorZ(
    VEC0 &&vec_x,
    VEC0 &&vec_y,
    const VEC1 &vec_n) {
  const double vec_s[3] = {0, 1, 0};
  Cross(vec_x, vec_s, vec_n);
  const double len = Length3(vec_x);
  if (len < 1.0e-10) {
    const double vec_t[3] = {1, 0, 0};
    Cross(vec_x, vec_t, vec_n);  // z????
    Cross(vec_y, vec_n, vec_x);  // x????
  } else {
    const double invlen = 1.0 / len;
    vec_x[0] *= invlen;
    vec_x[1] *= invlen;
    vec_x[2] *= invlen;
    Cross(vec_y, vec_n, vec_x);
  }
}


// =====================================
// dependency to the std::vector

template<typename VEC>
std::ostream &operator<<(
    std::ostream &output,
    const std::vector<VEC> &aV) {
  output << aV.size() << std::endl;
  for (unsigned int iv = 0; iv < aV.size(); ++iv) {
    output << "  " << iv << "-->" << aV[iv][0] << " " << aV[iv][1] << " " << aV[iv][2] << std::endl;
  }
  return output;
}

// -------------------------

template<typename VEC>
std::istream &operator>>(
    std::istream &input,
    std::vector<VEC> &aV) {
  int nV;
  input >> nV;
  aV.resize(nV);
  for (int iv = 0; iv < nV; iv++) {
    input >> aV[iv][0] >> aV[iv][1] >> aV[iv][2];
  }
  return input;
}

}

#ifdef DFM2_STATIC_LIBRARY
#include "delfem2/vec3.h"
template float delfem2::Dot3(const float (&)[3], const float (&)[3]);
template double delfem2::Dot3(const double (&)[3], const double (&)[3]);
template float delfem2::Dot3(const float (&)[], const float (&)[]);
template double delfem2::Dot3(const double (&)[], const double (&)[]);
template double delfem2::Dot3(const CVec3d&, const CVec3d&);
template float delfem2::Dot3(const CVec3f&, const CVec3f&);
template float delfem2::Dot3(const std::array<float,3>&, const std::array<float,3>&);
template double delfem2::Dot3(const std::array<double,3>&, const std::array<double,3>&);
//
template float delfem2::Distance3(const float (&)[3], const float (&)[3]);
template double delfem2::Distance3(const double (&)[3], const double (&)[3]);
template float delfem2::Distance3(const float (&)[], const float (&)[]);
template double delfem2::Distance3(const double (&)[], const double (&)[]);
template double delfem2::Distance3(const CVec3d&, const CVec3d&);
template float delfem2::Distance3(const CVec3f&, const CVec3f&);
template float delfem2::Distance3(const std::array<float,3>&, const std::array<float,3>&);
template double delfem2::Distance3(const std::array<double,3>&, const std::array<double,3>&);
//
template float delfem2::Length3(const float (&)[3]);
template double delfem2::Length3(const double (&)[3]);
template float delfem2::Length3(const float (&)[]);
template double delfem2::Length3(const double (&)[]);
template float delfem2::Length3(const CVec3f&);
template double delfem2::Length3(const CVec3d&);
template float delfem2::Length3(const std::array<float,3>&);
template double delfem2::Length3(const std::array<double,3>&);
//
template void delfem2::Normalize3(float (&)[3]);
template void delfem2::Normalize3(double (&)[3]);
template void delfem2::Normalize3(float (&)[]);
template void delfem2::Normalize3(double (&)[]);
template void delfem2::Normalize3(CVec3f&);
template void delfem2::Normalize3(CVec3d&);
template void delfem2::Normalize3(std::array<float,3>&);
template void delfem2::Normalize3(std::array<double,3>&);
//
template delfem2::CVec3d delfem2::Cross(const CVec3d&, const CVec3d&);
template delfem2::CVec3f delfem2::Cross(const CVec3f&, const CVec3f&);
template std::array<float,3> delfem2::Cross(const std::array<float,3>&, const std::array<float,3>&);
template std::array<double,3> delfem2::Cross(const std::array<double,3>&, const std::array<double,3>&);
//
template float delfem2::ScalarTripleProduct(const float (&)[3], const float (&)[3], const float (&)[3]);
template double delfem2::ScalarTripleProduct(const double (&)[3], const double (&)[3], const double (&)[3]);
template float delfem2::ScalarTripleProduct(const float (&)[], const float (&)[], const float (&)[]);
template double delfem2::ScalarTripleProduct(const double (&)[], const double (&)[], const double (&)[]);
template float delfem2::ScalarTripleProduct(const CVec3f&, const CVec3f&, const CVec3f&);
template float delfem2::ScalarTripleProduct(const CVec3d&, const CVec3d&, const CVec3d&);
template float delfem2::ScalarTripleProduct(const std::array<float,3>&, const std::array<float,3>&, const std::array<float,3>&);
template double delfem2::ScalarTripleProduct(const std::array<double,3>&, const std::array<double,3>&, const std::array<double,3>&);
//
template void delfem2::Cross(float (&)[3], const float (&)[3], const float (&)[3]);
template void delfem2::Cross(double (&)[3], const double (&)[3], const double (&)[3]);
template void delfem2::Cross(float (&)[], const float (&)[], const float (&)[]);
template void delfem2::Cross(double (&)[], const double (&)[], const double (&)[]);
#endif

#endif // DFM2_GEO_VEC3_H
