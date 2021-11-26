/*
* Copyright (c) 2019 Nobuyuki Umetani
*
* This source code is licensed under the MIT license found in the
* LICENSE file in the root directory of this source tree.
*/


#ifndef DFM2_GEOPROXIMITY3_V3_H
#define DFM2_GEOPROXIMITY3_V3_H

#include <cassert>
#include <cmath>
#include <iostream>
#include <vector>

#include "delfem2/vec3.h"
#include "delfem2/dfm2_inline.h"

#define NEARLY_ZERO 1.e-16

namespace delfem2 {

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
bool intersection_Plane_Line(
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

template<typename T>
CVec3<T> intersection_Plane_Line(
    const CVec3<T> &o, // one point on plane
    const CVec3<T> &n, // plane normal
    const CVec3<T> &s, // one point on line
    const CVec3<T> &d); // direction of line

// ----------------------------------------------------------------------------------

template<typename T>
double DistanceFaceVertex(
    const CVec3<T> &p0, const CVec3<T> &p1, const CVec3<T> &p2,
    const CVec3<T> &p3,
    double &w0, double &w1);

template<typename T>
double DistanceEdgeEdge(
    const CVec3<T> &p0, const CVec3<T> &p1,
    const CVec3<T> &q0, const CVec3<T> &q1,
    double &ratio_p, double &ratio_q);

template<typename T>
bool FindCoplanerInterp(
    double &r,
    const CVec3<T> &s0, const CVec3<T> &s1, const CVec3<T> &s2, const CVec3<T> &s3,
    const CVec3<T> &e0, const CVec3<T> &e1, const CVec3<T> &e2, const CVec3<T> &e3);

template<typename T>
bool IsContact_EE_Proximity(
    int ino0, int ino1, int jno0, int jno1,
    const CVec3<T> &p0, const CVec3<T> &p1, const CVec3<T> &q0, const CVec3<T> &q1,
    double delta);

template<typename T>
bool IsContact_FV_CCD2(
    int ino0, int ino1, int ino2, int ino3,
    const CVec3<T> &p0, const CVec3<T> &p1, const CVec3<T> &p2, const CVec3<T> &p3,
    const CVec3<T> &q0, const CVec3<T> &q1, const CVec3<T> &q2, const CVec3<T> &q3);

template<typename T>
bool isIntersectTriPair(
    CVec3<T> &P0, CVec3<T> &P1,
    int itri, int jtri,
    const std::vector<unsigned int> &aTri,
    const std::vector<double> &aXYZ);

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
bool isPointSameSide(
    const CVec3<T> &p0,
    const CVec3<T> &p1,
    const CVec3<T> &line_p0,
    const CVec3<T> &line_p1);

template<typename T>
bool isRayIntersectingTriangle(
    const CVec3<T> &line0,
    const CVec3<T> &line1,
    const CVec3<T> &tri0,
    const CVec3<T> &tri1,
    const CVec3<T> &tri2,
    CVec3<T> &intersectionPoint);

} // end namespace delfem2

#ifndef DFM2_STATIC_LIBRARY
#  include "delfem2/geoproximity3_v3.cpp"
#endif

#endif // DFM2_GEOPROXIMITY3_V3_H
