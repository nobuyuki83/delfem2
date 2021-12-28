/*
* Copyright (c) 2019 Nobuyuki Umetani
*
* This source code is licensed under the MIT license found in the
* LICENSE file in the root directory of this source tree.
*/

/**
 * @detail The order of dependency in delfem2:
 * line < ray < edge < polyline < quadratic < cubic < bspline << plane < tri < quad
 */

#ifndef DFM2_PLANE_H
#define DFM2_PLANE_H

#include <cassert>
#include <cmath>
#include <iostream>
#include <vector>

#include "delfem2/vec3.h"
#include "delfem2/dfm2_inline.h"

#define NEARLY_ZERO 1.e-16

namespace delfem2 {

    
template<typename T>
CVec3<T> Nearest_Plane3_Point3(
  const CVec3<T> &p, // point
  const CVec3<T> &o, // origin
  const CVec3<T> &n); // normal

template<typename T>
CVec3<T> Nearest_Origin3_PlaneTri3(
  T &r0,
  T &r1,
  const CVec3<T> &q0,
  const CVec3<T> &q1,
  const CVec3<T> &q2);


template<typename T>
bool Intersection_Plane3_Line3(
    CVec3<T> &p0,
    double &r0,
    double &r1,
    double &r2,
    double eps,
    const CVec3<T> &src,
    const CVec3<T> &dir,
    const CVec3<T> &q0,
    const CVec3<T> &q1,
    const CVec3<T> &q2);

template<typename T>
CVec3<T> Intersection_Plane3_Line3(
    const CVec3<T> &o, // one point on plane
    const CVec3<T> &n, // plane normal
    const CVec3<T> &s, // one point on line
    const CVec3<T> &d); // direction of line

// ----------------------------------------------------------------------------------

template<typename T>
bool IsInside_Orgin_BoundingBoxPoint4(
    const CVec3<T> &p0,
    const CVec3<T> &p1,
    const CVec3<T> &p2,
    const CVec3<T> &p3);

template<typename T>
bool IsInside_Orgin_BoundingBoxPoint5(
    const CVec3<T> &p0,
    const CVec3<T> &p1,
    const CVec3<T> &p2,
    const CVec3<T> &p3,
    const CVec3<T> &p4);

template<typename T>
bool IsInside_Orgin_BoundingBoxPoint6(
    const CVec3<T> &p0,
    const CVec3<T> &p1,
    const CVec3<T> &p2,
    const CVec3<T> &p3,
    const CVec3<T> &p4,
    const CVec3<T> &p5);

} // end namespace delfem2

#ifndef DFM2_STATIC_LIBRARY
#  include "delfem2/geo_plane.cpp"
#endif

#endif // DFM2_PLANE_H
