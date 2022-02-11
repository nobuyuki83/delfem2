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

}  // delfem2::geo_plane

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
  assert(((q1 - q0).cross(q2 - q0)).norm() > 1.0e-10);
  const CVec3<T> n1 = ((q1 - q0).cross(q2 - q0)).normalized();
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

// ----------------------------------------------------------------------

template<typename T>
bool delfem2::Intersection_Plane3_Line3(
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
template bool delfem2::Intersection_Plane3_Line3(
  CVec3d &p0, double &r0, double &r1, double &r2,
  double eps,
  const CVec3d &src, const CVec3d &dir,
  const CVec3d &q0, const CVec3d &q1, const CVec3d &q2);
#endif

// -------------------

template<typename T>
delfem2::CVec3<T> delfem2::Intersection_Plane3_Line3(
  const CVec3<T> &o, // one point on plane
  const CVec3<T> &n, // plane normal
  const CVec3<T> &s, // one point on line
  const CVec3<T> &d) // direction of line
{
  double t = ((o - s).dot(n)) / (d.dot(n));
  return s + t * d;
}
#ifdef DFM2_STATIC_LIBRARY
template delfem2::CVec3d delfem2::Intersection_Plane3_Line3(
  const CVec3d &o, // one point on plane
  const CVec3d &n, // plane normal
  const CVec3d &s, // one point on line
  const CVec3d &d); // direction of line
#endif


#ifdef DFM2_STATIC_LIBRARY

#include "delfem2/vec3.h"

namespace delfem2 {
  using f0 = double [3];
  using d0 = double [3];
  using f1 = float*;
  using d1 = double*;
  using f2 = std::array<float,3>;
  using d2 = std::array<double,3>;
  using f3 = CVec3f;
  using d3 = CVec3d;

template d3 delfem2::Nearest_Origin3_PlaneTri3(double &, double &, const d3 &, const d3 &q1, const d3 &);
template f3 delfem2::Nearest_Origin3_PlaneTri3(float &, float &, const f3 &, const f3 &q1, const f3 &);

}
#endif