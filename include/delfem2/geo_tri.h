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

} // end namespace delfem2

#ifndef DFM2_STATIC_LIBRARY
#  include "delfem2/geo_tri.cpp"
#endif

#endif // DFM2_GEO_TRI_H
