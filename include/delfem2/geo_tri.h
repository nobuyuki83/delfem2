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
 * plane < tri < quad
 */


#ifndef DFM2_GEO_TRI_H
#define DFM2_GEO_TRI_H

#include <cassert>
#include <cmath>
#include <iostream>
#include <vector>

#include "delfem2/vec3.h"
#include "delfem2/geo_vec3.h"
#include "delfem2/dfm2_inline.h"

#define NEARLY_ZERO 1.e-16

namespace delfem2 {

DFM2_INLINE void Nearest_Triangle3_Point3(
  double nearest_position[3],
  double &r0,
  double &r1,
  const double ps[3], // origin point
  const double q0[3],
  const double q1[3],
  const double q2[3]);


template<typename T>
CVec3<T> Nearest_Origin3_Tri3(
    T &r0,
    T &r1,
    const CVec3<T> &q0,
    const CVec3<T> &q1,
    const CVec3<T> &q2);

template<typename T>
bool IntersectRay_Tri3(
  T &r0,
  T &r1,
  const CVec3<T> &org,
  const CVec3<T> &dir,
  const CVec3<T> &p0,
  const CVec3<T> &p1,
  const CVec3<T> &p2,
  T eps);

template<typename T>
bool isIntersectTriPair(
    CVec3<T> &P0, CVec3<T> &P1,
    int itri, int jtri,
    const std::vector<unsigned int> &aTri,
    const std::vector<double> &aXYZ);

template<typename T>
CVec3<T> ProjectPointOnTriangle(
    const CVec3<T> &p0,
    const CVec3<T> &tri_p1,
    const CVec3<T> &tri_p2,
    const CVec3<T> &tri_p3);

template<typename T>
bool isPointInsideTriangle(
    const CVec3<T> &p0,
    const CVec3<T> &tri_p1,
    const CVec3<T> &tri_p2,
    const CVec3<T> &tri_p3);

template<typename T>
bool isRayIntersectingTriangle(
    const CVec3<T> &line0,
    const CVec3<T> &line1,
    const CVec3<T> &tri0,
    const CVec3<T> &tri1,
    const CVec3<T> &tri2,
    CVec3<T> &intersectionPoint);

template<typename T>
double DistanceFaceVertex(
    const CVec3<T> &p0, 
    const CVec3<T> &p1, 
    const CVec3<T> &p2,
    const CVec3<T> &p3,
    double &w0, 
    double &w1);

template<typename VEC, typename T = value_type<VEC>>
DFM2_INLINE T Area_Tri3(
  const VEC &v1,
  const VEC &v2,
  const VEC &v3) {
  const T x = (v2[1] - v1[1]) * (v3[2] - v1[2]) - (v3[1] - v1[1]) * (v2[2] - v1[2]);
  const T y = (v2[2] - v1[2]) * (v3[0] - v1[0]) - (v3[2] - v1[2]) * (v2[0] - v1[0]);
  const T z = (v2[0] - v1[0]) * (v3[1] - v1[1]) - (v3[0] - v1[0]) * (v2[1] - v1[1]);
  return std::sqrt(x * x + y * y + z * z)/2;
}

template<typename VEC, typename T = value_type<VEC>>
T Area_Tri3(
    const int iv1,
    const int iv2,
    const int iv3,
    const std::vector<VEC> &aPoint) {
  return Area_Tri3(aPoint[iv1], aPoint[iv2], aPoint[iv3]);
}

template<typename VEC, typename T = value_type<VEC>>
T SquareTriArea(
    const VEC &v1,
    const VEC &v2,
    const VEC &v3) {
  const T dtmp_x = (v2[1] - v1[1]) * (v3[2] - v1[2]) - (v2[2] - v1[2]) * (v3[1] - v1[1]);
  const T dtmp_y = (v2[2] - v1[2]) * (v3[0] - v1[0]) - (v2[0] - v1[0]) * (v3[2] - v1[2]);
  const T dtmp_z = (v2[0] - v1[0]) * (v3[1] - v1[1]) - (v2[1] - v1[1]) * (v3[0] - v1[0]);
  return (dtmp_x * dtmp_x + dtmp_y * dtmp_y + dtmp_z * dtmp_z) / 4;
}

template<typename VEC0, typename VEC1>
void Normal_Tri3(
    VEC0 &&vnorm,
    const VEC1 &v1,
    const VEC1 &v2,
    const VEC1 &v3) {
  vnorm[0] = (v2[1] - v1[1]) * (v3[2] - v1[2]) - (v2[2] - v1[2]) * (v3[1] - v1[1]);
  vnorm[1] = (v2[2] - v1[2]) * (v3[0] - v1[0]) - (v2[0] - v1[0]) * (v3[2] - v1[2]);
  vnorm[2] = (v2[0] - v1[0]) * (v3[1] - v1[1]) - (v2[1] - v1[1]) * (v3[0] - v1[0]);
}

template<typename VEC0, typename VEC1, typename T = value_type<VEC0>>
void UnitNormal_Tri3(
    VEC0 &&vnorm,
    const VEC1 &v1,
    const VEC1 &v2,
    const VEC1 &v3) {
  vnorm[0] = (v2[1] - v1[1]) * (v3[2] - v1[2]) - (v2[2] - v1[2]) * (v3[1] - v1[1]);
  vnorm[1] = (v2[2] - v1[2]) * (v3[0] - v1[0]) - (v2[0] - v1[0]) * (v3[2] - v1[2]);
  vnorm[2] = (v2[0] - v1[0]) * (v3[1] - v1[1]) - (v2[1] - v1[1]) * (v3[0] - v1[0]);
  const T dtmp1 = 1 / vnorm.norm();
  vnorm[0] *= dtmp1;
  vnorm[1] *= dtmp1;
  vnorm[2] *= dtmp1;
}


// -----------------------------------------------
// should moved to mshuni

template<typename T>
std::array<T, 3> Normal_Tri3(
    unsigned int itri,
    const std::vector<unsigned int> &aTri,
    const std::vector<T> &aXYZ) {
  const unsigned int i0 = aTri[itri * 3 + 0];
  const unsigned int i1 = aTri[itri * 3 + 1];
  const unsigned int i2 = aTri[itri * 3 + 2];
  const T p01[3] = {
      aXYZ[i1 * 3 + 0] - aXYZ[i0 * 3 + 0],
      aXYZ[i1 * 3 + 1] - aXYZ[i0 * 3 + 1],
      aXYZ[i1 * 3 + 2] - aXYZ[i0 * 3 + 2]
  };
  const T p02[3] = {
      aXYZ[i2 * 3 + 0] - aXYZ[i0 * 3 + 0],
      aXYZ[i2 * 3 + 1] - aXYZ[i0 * 3 + 1],
      aXYZ[i2 * 3 + 2] - aXYZ[i0 * 3 + 2]
  };
  std::array<T, 3> ret;
  Cross(ret.data(), p01, p02);
  return ret;
}

// --------------------

template<typename VEC0, typename VEC1, typename T = value_type<VEC0>>
void UnitNormalAreaTri3(
    VEC0 &&n,
    T &a,
    const VEC1 &v1, const VEC1 &v2, const VEC1 &v3) {
  Normal_Tri3(
      n,
      v1, v2, v3);
  a = Length3(n) / 2;
  const T invlen = 1 / (a * 2);
  n[0] *= invlen;
  n[1] *= invlen;
  n[2] *= invlen;
}

template<typename VEC>
VEC Normal_Tri3(
    const VEC &v1,
    const VEC &v2,
    const VEC &v3) {
  return {
      (v2[1] - v1[1]) * (v3[2] - v1[2]) - (v2[2] - v1[2]) * (v3[1] - v1[1]),
      (v2[2] - v1[2]) * (v3[0] - v1[0]) - (v2[0] - v1[0]) * (v3[2] - v1[2]),
      (v2[0] - v1[0]) * (v3[1] - v1[1]) - (v2[1] - v1[1]) * (v3[0] - v1[0]) };
}

template<typename VEC, typename T = value_type<VEC>>
VEC UnitNormal_Tri3(
    const VEC &v1,
    const VEC &v2,
    const VEC &v3) {
  VEC vnorm{
      (v2[1] - v1[1]) * (v3[2] - v1[2]) - (v2[2] - v1[2]) * (v3[1] - v1[1]),
      (v2[2] - v1[2]) * (v3[0] - v1[0]) - (v2[0] - v1[0]) * (v3[2] - v1[2]),
      (v2[0] - v1[0]) * (v3[1] - v1[1]) - (v2[1] - v1[1]) * (v3[0] - v1[0]) };
  const T dtmp1 = 1 / vnorm.norm();
  vnorm[0] *= dtmp1;
  vnorm[1] *= dtmp1;
  vnorm[2] *= dtmp1;
  return vnorm;
}

// ----------------------------------
// this can be moved to meshuni

template<typename T>
std::array<T, 3> CG_Tri3(
    unsigned int itri,
    const std::vector<unsigned int> &aTri,
    const std::vector<T> &aXYZ) {
  const unsigned int i0 = aTri[itri * 3 + 0];
  const unsigned int i1 = aTri[itri * 3 + 1];
  const unsigned int i2 = aTri[itri * 3 + 2];
  return {
      (aXYZ[i0 * 3 + 0] + aXYZ[i1 * 3 + 0] + aXYZ[i2 * 3 + 0]) / 3.0,
      (aXYZ[i0 * 3 + 1] + aXYZ[i1 * 3 + 1] + aXYZ[i2 * 3 + 1]) / 3.0,
      (aXYZ[i0 * 3 + 2] + aXYZ[i1 * 3 + 2] + aXYZ[i2 * 3 + 2]) / 3.0};
}

} // end namespace delfem2

#ifndef DFM2_STATIC_LIBRARY
#  include "delfem2/geo_tri.cpp"
#endif

#endif // DFM2_GEO_TRI_H
