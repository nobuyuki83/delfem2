/*
 * Copyright (c) 2019 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */

#include "delfem2/geo_plane.h"

#include <cmath>
#include <stack>

#ifndef M_PI
#  define M_PI 3.14159265358979323846
#endif

// =====================================
// below: unexposed 

namespace delfem2::geo_plane {

DFM2_INLINE bool MyIsnan(double x) { return x != x; }

// evaluate cubic function
template<typename REAL>
DFM2_INLINE REAL EvaluateCubic(
  REAL x,
  REAL k0, REAL k1, REAL k2, REAL k3) // coefficient of cubic function
{
  return k0 + k1 * x + k2 * x * x + k3 * x * x * x;
}
#ifdef DFM2_STATIC_LIBRARY
template float EvaluateCubic(float r2, float k0, float k1, float k2, float k3);
template double EvaluateCubic(double r2, double k0, double k1, double k2, double k3);
#endif


// find root of cubic function using bisection method
DFM2_INLINE double FindRootCubic_Bisect(
  double r0, double r1,
  double v0, double v1,
  double k0, double k1, double k2, double k3) {
  assert(v0 * v1 <= 0);
  if (v0 * v1 == 0) {
    if (v0 == 0) { return r0; }
    else { return r1; }
  }
  for (unsigned int itr = 0; itr < 15; itr++) {
    const double r2 = 0.5 * (r0 + r1);
    const double v2 = EvaluateCubic(r2, k0, k1, k2, k3);
    if (v2 == 0) { return r2; }
    if (v0 * v2 < 0) {
      r1 = r2;
      v1 = v2;
    } else {
      r0 = r2;
      v0 = v2;
    }
  }
  return 0.5 * (r0 + r1);
}

template<typename REAL>
DFM2_INLINE void MyInverse_Mat3
  (REAL Ainv[9],
   const REAL A[9]) {
  const REAL det =
    +A[0] * A[4] * A[8] + A[3] * A[7] * A[2] + A[6] * A[1] * A[5]
      - A[0] * A[7] * A[5] - A[6] * A[4] * A[2] - A[3] * A[1] * A[8];
  const REAL inv_det = 1.0 / det;
  Ainv[0] = inv_det * (A[4] * A[8] - A[5] * A[7]);
  Ainv[1] = inv_det * (A[2] * A[7] - A[1] * A[8]);
  Ainv[2] = inv_det * (A[1] * A[5] - A[2] * A[4]);
  Ainv[3] = inv_det * (A[5] * A[6] - A[3] * A[8]);
  Ainv[4] = inv_det * (A[0] * A[8] - A[2] * A[6]);
  Ainv[5] = inv_det * (A[2] * A[3] - A[0] * A[5]);
  Ainv[6] = inv_det * (A[3] * A[7] - A[4] * A[6]);
  Ainv[7] = inv_det * (A[1] * A[6] - A[0] * A[7]);
  Ainv[8] = inv_det * (A[0] * A[4] - A[1] * A[3]);
}

template<typename T>
DFM2_INLINE void MyMatVec3
  (T y[3],
   const T m[9], const T x[3]) {
  y[0] = m[0] * x[0] + m[1] * x[1] + m[2] * x[2];
  y[1] = m[3] * x[0] + m[4] * x[1] + m[5] * x[2];
  y[2] = m[6] * x[0] + m[7] * x[1] + m[8] * x[2];
}

//! Volume of a tetrahedra
template<typename T>
T Volume_Tet(
  const CVec3<T> &v0,
  const CVec3<T> &v1,
  const CVec3<T> &v2,
  const CVec3<T> &v3) {
//  return delfem2::Volume_Tet3(v0.p, v1.p, v2.p, v3.p);
  T v = (v1.p[0] - v0.p[0]) * ((v2.p[1] - v0.p[1]) * (v3.p[2] - v0.p[2]) - (v3.p[1] - v0.p[1]) * (v2.p[2] - v0.p[2]))
    + (v1.p[1] - v0.p[1]) * ((v2.p[2] - v0.p[2]) * (v3.p[0] - v0.p[0]) - (v3.p[2] - v0.p[2]) * (v2.p[0] - v0.p[0]))
    + (v1.p[2] - v0.p[2]) * ((v2.p[0] - v0.p[0]) * (v3.p[1] - v0.p[1]) - (v3.p[0] - v0.p[0]) * (v2.p[1] - v0.p[1]));
  return v * static_cast<T>(1.0 / 6.0);
}

template<typename T>
T Volume_OrgTet(
  const CVec3<T> &v1,
  const CVec3<T> &v2,
  const CVec3<T> &v3) {
  double v =
    v1.p[0] * (v2.p[1] * v3.p[2] - v3.p[1] * v2.p[2])
      + v1.p[1] * (v2.p[2] * v3.p[0] - v3.p[2] * v2.p[0])
      + v1.p[2] * (v2.p[0] * v3.p[1] - v3.p[0] * v2.p[1]);
  return v * 0.16666666666666666666666666666667;
};

template<typename T>
T Volume_Tet3(
  const T v1[3],
  const T v2[3],
  const T v3[3],
  const T v4[3]) {
  return
    (
      +(v2[0] - v1[0]) * ((v3[1] - v1[1]) * (v4[2] - v1[2]) - (v4[1] - v1[1]) * (v3[2] - v1[2]))
        - (v2[1] - v1[1]) * ((v3[0] - v1[0]) * (v4[2] - v1[2]) - (v4[0] - v1[0]) * (v3[2] - v1[2]))
        + (v2[2] - v1[2]) * ((v3[0] - v1[0]) * (v4[1] - v1[1]) - (v4[0] - v1[0]) * (v3[1] - v1[1]))
    ) * 0.16666666666666666666666666666667;
}

}

// ---------------------------------------------------------------------------


template<typename T>
delfem2::CVec3<T> delfem2::Nearest_Plane3_Point3(
  const CVec3<T> &p, // point
  const CVec3<T> &o, // origin
  const CVec3<T> &n) // normal
{
  const CVec3<T> n0 = n.normalized();
  return p + ((o - p) * n0) * n0;
}

// -------------------------------------

template<typename T>
delfem2::CVec3<T> delfem2::Nearest_Origin3_PlaneTri3(
  T &r0,
  T &r1,
  const CVec3<T> &q0,
  const CVec3<T> &q1,
  const CVec3<T> &q2) {
  namespace lcl = delfem2::geo_plane;
  assert((q1 - q0).norm() > 1.0e-10);
  assert((q2 - q0).norm() > 1.0e-10);
  assert((q1 - q2).norm() > 1.0e-10);
  assert(((q1 - q0) ^ (q2 - q0)).norm() > 1.0e-10);
  const CVec3<T> n1 = ((q1 - q0) ^ (q2 - q0)).normalized();
  const T v0 = lcl::Volume_OrgTet(q1, q2, n1);
  const T v1 = lcl::Volume_OrgTet(q2, q0, n1);
  const T v2 = lcl::Volume_OrgTet(q0, q1, n1);
  assert(fabs(v0 + v1 + v2) > 1.0e-10);
  T vt_inv = 1 / (v0 + v1 + v2);
  T r2;
  r0 = v0 * vt_inv;
  r1 = v1 * vt_inv;
  r2 = v2 * vt_inv;
  return q0 * r0 + q1 * r1 + q2 * r2;
}
#ifdef DFM2_STATIC_LIBRARY
template delfem2::CVec3d delfem2::Nearest_Origin3_PlaneTri3(
  double &r0,
  double &r1,
  const CVec3d &q0,
  const CVec3d &q1,
  const CVec3d &q2);
#endif


// ------------------------------------------------------------------------------

template<typename T>
bool delfem2::IsInside_Orgin_BoundingBoxPoint6(
  const CVec3<T> &p0,
  const CVec3<T> &p1,
  const CVec3<T> &p2,
  const CVec3<T> &p3,
  const CVec3<T> &p4,
  const CVec3<T> &p5) {
  if (p0.p[0] > 0 && p1.p[0] > 0 && p2.p[0] > 0 && p3.p[0] > 0 && p4.p[0] > 0 && p5.p[0] > 0) { return false; }
  if (p0.p[0] < 0 && p1.p[0] < 0 && p2.p[0] < 0 && p3.p[0] < 0 && p4.p[0] < 0 && p5.p[0] < 0) { return false; }
  if (p0.p[1] > 0 && p1.p[1] > 0 && p2.p[1] > 0 && p3.p[1] > 0 && p4.p[1] > 0 && p5.p[1] > 0) { return false; }
  if (p0.p[1] < 0 && p1.p[1] < 0 && p2.p[1] < 0 && p3.p[1] < 0 && p4.p[1] < 0 && p5.p[1] < 0) { return false; }
  if (p0.p[2] > 0 && p1.p[2] > 0 && p2.p[2] > 0 && p3.p[2] > 0 && p4.p[2] > 0 && p5.p[2] > 0) { return false; }
  if (p0.p[2] < 0 && p1.p[2] < 0 && p2.p[2] < 0 && p3.p[2] < 0 && p4.p[2] < 0 && p5.p[2] < 0) { return false; }
  return true;
}

template<typename T>
bool delfem2::IsInside_Orgin_BoundingBoxPoint5
  (const CVec3<T> &p0,
   const CVec3<T> &p1,
   const CVec3<T> &p2,
   const CVec3<T> &p3,
   const CVec3<T> &p4) {
  if (p0.p[0] > 0 && p1.p[0] > 0 && p2.p[0] > 0 && p3.p[0] > 0 && p4.p[0] > 0) { return false; }
  if (p0.p[0] < 0 && p1.p[0] < 0 && p2.p[0] < 0 && p3.p[0] < 0 && p4.p[0] < 0) { return false; }
  if (p0.p[1] > 0 && p1.p[1] > 0 && p2.p[1] > 0 && p3.p[1] > 0 && p4.p[1] > 0) { return false; }
  if (p0.p[1] < 0 && p1.p[1] < 0 && p2.p[1] < 0 && p3.p[1] < 0 && p4.p[1] < 0) { return false; }
  if (p0.p[2] > 0 && p1.p[2] > 0 && p2.p[2] > 0 && p3.p[2] > 0 && p4.p[2] > 0) { return false; }
  if (p0.p[2] < 0 && p1.p[2] < 0 && p2.p[2] < 0 && p3.p[2] < 0 && p4.p[2] < 0) { return false; }
  return true;
}

template<typename T>
bool delfem2::IsInside_Orgin_BoundingBoxPoint4
  (const CVec3<T> &p0,
   const CVec3<T> &p1,
   const CVec3<T> &p2,
   const CVec3<T> &p3) {
  if (p0.p[0] > 0 && p1.p[0] > 0 && p2.p[0] > 0 && p3.p[0] > 0) { return false; }
  if (p0.p[0] < 0 && p1.p[0] < 0 && p2.p[0] < 0 && p3.p[0] < 0) { return false; }
  if (p0.p[1] > 0 && p1.p[1] > 0 && p2.p[1] > 0 && p3.p[1] > 0) { return false; }
  if (p0.p[1] < 0 && p1.p[1] < 0 && p2.p[1] < 0 && p3.p[1] < 0) { return false; }
  if (p0.p[2] > 0 && p1.p[2] > 0 && p2.p[2] > 0 && p3.p[2] > 0) { return false; }
  if (p0.p[2] < 0 && p1.p[2] < 0 && p2.p[2] < 0 && p3.p[2] < 0) { return false; }
  return true;
}

// ----------------------------------------------------------------------


// ----------------

template<typename T>
bool delfem2::intersection_Plane_Line(
  CVec3<T> &p0,
  double &r0,
  double &r1,
  double &r2,
  double eps,
  const CVec3<T> &src,
  const CVec3<T> &dir,
  const CVec3<T> &q0,
  const CVec3<T> &q1,
  const CVec3<T> &q2) {
  namespace lcl = delfem2::geo_plane;
  r0 = lcl::Volume_Tet(src, src + dir, q1, q2);
  r1 = lcl::Volume_Tet(src, src + dir, q2, q0);
  r2 = lcl::Volume_Tet(src, src + dir, q0, q1);
  double v012 = (r0 + r1 + r2);
  double v012_inv = 1.0 / v012;
  r0 *= v012_inv;
  r1 *= v012_inv;
  r2 *= v012_inv;
  p0 = r0 * q0 + r1 * q1 + r2 * q2;
  return r0 > eps && r1 > eps && r2 > eps;
}
#ifdef DFM2_STATIC_LIBRARY
template bool delfem2::intersection_Plane_Line(
  CVec3d &p0, double &r0, double &r1, double &r2,
  double eps,
  const CVec3d &src, const CVec3d &dir,
  const CVec3d &q0, const CVec3d &q1, const CVec3d &q2);
#endif

// -------------------

template<typename T>
delfem2::CVec3<T> delfem2::intersection_Plane_Line(
  const CVec3<T> &o, // one point on plane
  const CVec3<T> &n, // plane normal
  const CVec3<T> &s, // one point on line
  const CVec3<T> &d) // direction of line
{
  double t = ((o - s).dot(n)) / (d.dot(n));
  return s + t * d;
}
#ifdef DFM2_STATIC_LIBRARY
template delfem2::CVec3d delfem2::intersection_Plane_Line(
  const CVec3d &o, // one point on plane
  const CVec3d &n, // plane normal
  const CVec3d &s, // one point on line
  const CVec3d &d); // direction of line
#endif



// ------------------

// ----------------------------------

// ----------------------------------------------------------------------------


