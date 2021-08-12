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

DFM2_INLINE void GetNearest_LineSegPoint3D(
    double pn[3],
    const double p[3], // point
    const double s[3], // source
    const double e[3]); // end

DFM2_INLINE void GetNearest_TrianglePoint3D(
    double pn[3],
    double& r0,
    double& r1,
    const double ps[3], // origin point
    const double q0[3],
    const double q1[3],
    const double q2[3]);

template <typename T>
bool IntersectRay_Tri3(
    T& r0,
    T& r1,
    const CVec3<T>& org,
    const CVec3<T>& dir,
    const CVec3<T>& p0,
    const CVec3<T>& p1,
    const CVec3<T>& p2,
    T eps);

// --------------------------------------------------------------

template <typename T>
CVec3<T> nearest_Line_Point(
    const CVec3<T>& p, // point
    const CVec3<T>& s, // source
    const CVec3<T>& d);

template <typename T>
CVec3<T> nearest_Line_Point(
    double& t,
    const CVec3<T>& p, // point
    const CVec3<T>& s, // source
    const CVec3<T>& d); // direction

template <typename T>
CVec3<T> nearest_LineSeg_Point(
    const CVec3<T>& p, // point
    const CVec3<T>& s, // source
    const CVec3<T>& e); // end
  
template <typename T>
CVec3<T> nearest_LineSeg_Point(
    double& t,
    const CVec3<T>& p, // point
    const CVec3<T>& s,
    const CVec3<T>& e); // direction

template <typename T>
CVec3<T> nearest_Plane_Point(
    const CVec3<T>& p, // point
    const CVec3<T>& o, // origin
    const CVec3<T>& n); // normal
  
template <typename T>
CVec3<T> Nearest_Origin_Tri(
    double& r0,
    double& r1,
    const CVec3<T>& q0,
    const CVec3<T>& q1,
    const CVec3<T>& q2);
  
template <typename T>
CVec3<T> nearst_Origin_Quad(
    double& s0, double& s1,
    const CVec3<T>& p,
    const CVec3<T>& q0,
    const CVec3<T>& q1,
    const CVec3<T>& q2,
    const CVec3<T>& q3);

// --------------------------------

template <typename T>
void nearest_LineSeg_Line(
    CVec3<T>& a,
    CVec3<T>& b,
    const CVec3<T>& ps,
    const CVec3<T>& pe,
    const CVec3<T>& pb_,
    const CVec3<T>& vb);

// --------------------------------

/**
 *  @param D (out) scaling factor
 *  @param Da (out) nearest point scaled by D on line A
 *  @param Db (out) nearest point scaled by D on line B
 */
template <typename T>
void nearest_Line_Line(
    T& D,
    CVec3<T>& Da,
    CVec3<T>& Db,
    const CVec3<T>& pa_,
    const CVec3<T>& va,
    const CVec3<T>& pb_,
    const CVec3<T>& vb);

// ------------------------------------

DFM2_INLINE double Nearest_LineSeg_LineSeg_CCD_Iteration(
    double p[3],
    const CVec3d& p0s,
    const CVec3d& p0e,
    const CVec3d& p1s,
    const CVec3d& p1e,
    const CVec3d& q0s,
    const CVec3d& q0e,
    const CVec3d& q1s,
    const CVec3d& q1e,
    unsigned int nitr );

// ------------------------------------
  
/**
 *  @param D (out) scaling factor
 *  @param Da (out) nearest point scaled by D on line A
 *  @param Db (out) nearest point scaled by D on line B
 *  @param Dta (out) parameter for nearest pont one line A. Da = D*pa_ + Dta*va
 *  @param Dtb (out) parameter for nearest pont one line B. Db = D*pb_ + Dtb*vb
 */
template <typename T>
void nearest_Line_Line(
    T& D, CVec3<T>& Da, CVec3<T>& Db,
    T& Dta, T& Dtb,
    const CVec3<T>& pa_, const CVec3<T>& va,
    const CVec3<T>& pb_, const CVec3<T>& vb);

// ------------------------------------
  
/**
 * @param p0 (out)  nearest point on line
 * @param q0 (out)  nearest point on circle
 */
template <typename T>
void Nearest_Line_Circle(
    CVec3<T>& p0,
    CVec3<T>& q0,
    const CVec3<T>& src,
    const CVec3<T>& dir,
    const CVec3<T>& org,
    const CVec3<T>& normal,
    T rad);

template <typename T>
CVec3<T> nearst_Origin_Quad(
    double& s0,
    double& s1,
    const CVec3<T>& q0,
    const CVec3<T>& q1,
    const CVec3<T>& q2,
    const CVec3<T>& q3);

template <typename T>
CVec3<T> nearest_Origin_LineSeg(
    const CVec3<T>& s, // start
    const CVec3<T>& e); // end

// r0==0 -> p0==org
// r0==1 -> p1==org
template <typename T>
CVec3<T> nearest_Origin_LineSeg(
    double& r0,
    const CVec3<T>& p0, // start
    const CVec3<T>& p1); // end

template <typename T>
CVec3<T> Nearest_Orgin_PlaneTri(
    double& r0,
    double& r1,
    const CVec3<T>& q0,
    const CVec3<T>& q1,
    const CVec3<T>& q2);
    
// -----------------------------------------

template <typename T>
bool intersection_Plane_Line(
    CVec3<T>& p0,
    double& r0,
    double& r1,
    double& r2,
    double eps,
    const CVec3<T>& src,
    const CVec3<T>& dir,
    const CVec3<T>& q0,
    const CVec3<T>& q1,
    const CVec3<T>& q2);
  
template <typename T>
bool intersection_Point_Quad(
    CVec3<T>& psec,
    double& s0,
    double& s1,
    const CVec3<T>& src,
    const CVec3<T>& dir,
    const CVec3<T>& q0,
    const CVec3<T>& q1,
    const CVec3<T>& q2,
    const CVec3<T>& q3);

template <typename T>
void iteration_intersection_Line_Quad(
    double& t0, double& t1,
    const CVec3<T>& src,
    const CVec3<T>& u,
    const CVec3<T>& v,
    const CVec3<T>& q0,
    const CVec3<T>& q1,
    const CVec3<T>& q2,
    const CVec3<T>& q3);
  
template <typename T>
CVec3<T> intersection_Plane_Line(
    const CVec3<T>& o, // one point on plane
    const CVec3<T>& n, // plane normal
    const CVec3<T>& s, // one point on line
    const CVec3<T>& d); // direction of line

// ----------------------------------------------------------------------------------

template <typename T>
double DistanceFaceVertex(
    const CVec3<T>& p0, const CVec3<T>& p1, const CVec3<T>& p2,
    const CVec3<T>& p3,
    double& w0, double& w1);

template <typename T>
double DistanceEdgeEdge(
    const CVec3<T>& p0, const CVec3<T>& p1,
    const CVec3<T>& q0, const CVec3<T>& q1,
    double& ratio_p, double& ratio_q);

template <typename T>
bool FindCoplanerInterp(
    double& r,
    const CVec3<T>& s0, const CVec3<T>& s1, const CVec3<T>& s2, const CVec3<T>& s3,
    const CVec3<T>& e0, const CVec3<T>& e1, const CVec3<T>& e2, const CVec3<T>& e3);

template <typename T>
bool IsContact_EE_Proximity(
    int ino0,        int ino1,        int jno0,        int jno1,
    const CVec3<T>& p0, const CVec3<T>& p1, const CVec3<T>& q0, const CVec3<T>& q1,
    const double delta);

template <typename T>
bool IsContact_FV_CCD2(
    int ino0,        int ino1,        int ino2,        int ino3,
    const CVec3<T>& p0, const CVec3<T>& p1, const CVec3<T>& p2, const CVec3<T>& p3,
    const CVec3<T>& q0, const CVec3<T>& q1, const CVec3<T>& q2, const CVec3<T>& q3);

template <typename T>
bool isIntersectTriPair(
    CVec3<T>& P0, CVec3<T>& P1,
    int itri, int jtri,
    const std::vector<unsigned int>& aTri,
    const std::vector<double>& aXYZ);

template <typename T>
bool IsInside_Orgin_BoundingBoxPoint4(
    const CVec3<T>& p0,
    const CVec3<T>& p1,
    const CVec3<T>& p2,
    const CVec3<T>& p3);

template <typename T>
bool IsInside_Orgin_BoundingBoxPoint5(
    const CVec3<T>& p0,
    const CVec3<T>& p1,
    const CVec3<T>& p2,
    const CVec3<T>& p3,
    const CVec3<T>& p4);

template <typename T>
bool IsInside_Orgin_BoundingBoxPoint6(
    const CVec3<T>& p0,
    const CVec3<T>& p1,
    const CVec3<T>& p2,
    const CVec3<T>& p3,
    const CVec3<T>& p4,
    const CVec3<T>& p5);


template <typename T>
CVec3<T> ProjectPointOnTriangle(
    const CVec3<T> &p0,
    const CVec3<T> &tri_p1,
    const CVec3<T> &tri_p2,
    const CVec3<T> &tri_p3);

template <typename T>
bool isPointInsideTriangle(
    const CVec3<T> &p0,
    const CVec3<T> &tri_p1,
    const CVec3<T> &tri_p2,
    const CVec3<T> &tri_p3);
  
template <typename T>
bool isPointSameSide(
    const CVec3<T> &p0,
    const CVec3<T> &p1,
    const CVec3<T> &line_p0,
    const CVec3<T> &line_p1);
  
template <typename T>
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


#endif // VEC3_H
