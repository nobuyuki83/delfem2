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


#ifndef DFM2_GEO_QUAD_H
#define DFM2_GEO_QUAD_H

#include <cassert>
#include <cmath>
#include <iostream>
#include <vector>

#include "delfem2/vec3.h"
#include "delfem2/dfm2_inline.h"

#define NEARLY_ZERO 1.e-16

namespace delfem2 {

template<typename T>
T QuadBilinear(
  int iq,
  double r0,
  double r1,
  std::vector<int> &aQuad,
  std::vector<T> &aPoint) {
  int i0 = aQuad[iq * 4 + 0];
  int i1 = aQuad[iq * 4 + 1];
  int i2 = aQuad[iq * 4 + 2];
  int i3 = aQuad[iq * 4 + 3];
  return (1 - r0) * (1 - r1) * aPoint[i0]
    + r0 * (1 - r1) * aPoint[i1]
    + r0 * r1 * aPoint[i2]
    + (1 - r0) * r1 * aPoint[i3];
}

template<typename T>
DFM2_INLINE CVec3<T> Nearst_Origin3_Quad3(
  double &s0,
  double &s1,
  const CVec3<T> &q0,
  const CVec3<T> &q1,
  const CVec3<T> &q2,
  const CVec3<T> &q3);

template<typename T>
bool intersection_Point_Quad(
    CVec3<T> &psec,
    double &s0,
    double &s1,
    const CVec3<T> &src,
    const CVec3<T> &dir,
    const CVec3<T> &q0,
    const CVec3<T> &q1,
    const CVec3<T> &q2,
    const CVec3<T> &q3);

template<typename T>
void iteration_intersection_Line_Quad(
    double &t0, double &t1,
    const CVec3<T> &src,
    const CVec3<T> &u,
    const CVec3<T> &v,
    const CVec3<T> &q0,
    const CVec3<T> &q1,
    const CVec3<T> &q2,
    const CVec3<T> &q3);

} // end namespace delfem2

#ifndef DFM2_STATIC_LIBRARY
#  include "delfem2/geo_quad.cpp"
#endif

#endif // DFM2_GEO_NEAREST3_H
