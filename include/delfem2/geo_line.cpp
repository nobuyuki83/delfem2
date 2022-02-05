/*
 * Copyright (c) 2019 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */

#include "delfem2/geo_line.h"

#include <cmath>
#include <stack>

#include "delfem2/geo_vec3.h"

#ifndef M_PI
#  define M_PI 3.14159265358979323846
#endif

// =====================================
// below: unexposed 

namespace delfem2::geo_nearest3 {

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

template<class VEC>
VEC delfem2::Nearest_Line3_Point3(
  const VEC &point, // point
  const VEC &line_src, // source
  const VEC &line_dir) // direction
{
  using SCALAR = decltype(point[0]);
  assert(line_dir.dot(line_dir) > 1.0e-20);
  const VEC ps = line_src - point;
  const SCALAR a = line_dir.dot(line_dir);
  const SCALAR b = line_dir.dot(line_src - point);
  const SCALAR t = -b / a;
  return line_src + t * line_dir;
}
#ifdef DFM2_STATIC_LIBRARY
template delfem2::CVec3d delfem2::Nearest_Line3_Point3<delfem2::CVec3d>(
  const CVec3d &p,
  const CVec3d &s,
  const CVec3d &d);
template delfem2::CVec3f delfem2::Nearest_Line3_Point3<delfem2::CVec3f>(
  const CVec3f &p,
  const CVec3f &s,
  const CVec3f &d);
#endif

// -------------------------------------------------------

template<typename T>
delfem2::CVec3<T> delfem2::Nearest_Line3_Point3(
  double &t,
  const CVec3<T> &p, // point
  const CVec3<T> &s, // source
  const CVec3<T> &d) // direction
{
  if (Dot(d, d) < 1.0e-20) {
    t = 0;
    return s;
  }
  const CVec3<T> ps = s - p;
  double a = Dot(d, d);
  double b = Dot(d, s - p);
  t = -b / a;
  return s + t * d;
}

// ---------------------------------------------

template<typename VEC, typename T>
void delfem2::Nearest_Line3_Line3(
  T &scale,
  VEC &scaled_neraest_a,
  VEC &scaled_nearest_b,
  const VEC &line_org_a,
  const VEC &line_dir_a,
  const VEC &line_org_b,
  const VEC &line_dir_b) {
  const T xaa = line_dir_a.dot(line_dir_a);
  const T xab = line_dir_b.dot(line_dir_a);
  const T xbb = line_dir_b.dot(line_dir_b);
  scale = (xaa * xbb - xab * xab);
  const T xac = line_dir_a.dot(line_org_b - line_org_a);
  const T xbc = line_dir_b.dot(line_org_b - line_org_a);
  const T da = xbb * xac - xab * xbc;
  const T db = xab * xac - xaa * xbc;
  scaled_neraest_a = scale * line_org_a + da * line_dir_a;
  scaled_nearest_b = scale * line_org_b + db * line_dir_b;
}
#ifdef DFM2_STATIC_LIBRARY
template void delfem2::Nearest_Line3_Line3(
  float &, CVec3f &, CVec3f &,
  const CVec3f &, const CVec3f &,
  const CVec3f &, const CVec3f &);
template void delfem2::Nearest_Line3_Line3(
  double &, CVec3d &, CVec3d &,
  const CVec3d &, const CVec3d &,
  const CVec3d &, const CVec3d &);
#endif

// ---------------------------------------------

template<class VEC, typename T>
void delfem2::Nearest_Line3_Line3(
  T &D,
  VEC &Da,
  VEC &Db,
  T &Dta,
  T &Dtb,
  const VEC &pa_,
  const VEC &va,
  const VEC &pb_,
  const VEC &vb) {
  T xaa = va.dot(va);
  T xab = vb.dot(va);
  T xbb = vb.dot(vb);
  D = (xaa * xbb - xab * xab);
  T xac = va.dot(pb_ - pa_);
  T xbc = vb.dot(pb_ - pa_);
  Dta = xbb * xac - xab * xbc;
  Dtb = xab * xac - xaa * xbc;
  Da = D * pa_ + Dta * va;
  Db = D * pb_ + Dtb * vb;
}
#ifdef DFM2_STATIC_LIBRARY
template void delfem2::Nearest_Line3_Line3(
  float &D, CVec3f &Da, CVec3f &Db,
  float &Dta, float &Dtb,
  const CVec3f &pa_, const CVec3f &va,
  const CVec3f &pb_, const CVec3f &vb);
template void delfem2::Nearest_Line3_Line3(
  double &D, CVec3d &Da, CVec3d &Db,
  double &Dta, double &Dtb,
  const CVec3d &pa_, const CVec3d &va,
  const CVec3d &pb_, const CVec3d &vb);
#endif

// ------------------------------------------

template<typename T>
void delfem2::Nearest_Line_Circle(
  CVec3<T> &p0,
  CVec3<T> &q0,
  const CVec3<T> &src,
  const CVec3<T> &dir,
  const CVec3<T> &org, // center of the circle
  const CVec3<T> &normal, // normal of the circle
  T rad) {
  const int nitr = 4;
  // ---------------------------------------
  CVec3<T> ex, ey;
  FrameFromVectorZ(ex, ey, normal);
  T u0;
  {
    if (fabs(dir.dot(normal)) > fabs((org - src).dot(normal)) * 1.0e-4) {
      u0 = ((org - src).dot(normal)) / (dir.dot(normal));
    } else {
      u0 = (org - src).dot(dir) / (dir.dot(dir));
    }
  }
  for (int itr = 0; itr < nitr; ++itr) {
    p0 = src + u0 * dir;
    double t0 = atan2(ey.dot(p0 - org), ex.dot(p0 - org));
    q0 = (T) (rad * cos(t0)) * ex + (T) (rad * sin(t0)) * ey + org;
    u0 = (q0 - src).dot(dir) / (dir.dot(dir));
  }
}
#ifdef DFM2_STATIC_LIBRARY
template void delfem2::Nearest_Line_Circle(
  CVec3f &p0,
  CVec3f &q0,
  const CVec3f &src,
  const CVec3f &dir,
  const CVec3f &org, // center of the circle
  const CVec3f &normal, // normal of the circle
  float rad);
template void delfem2::Nearest_Line_Circle(
  CVec3d &p0,
  CVec3d &q0,
  const CVec3d &src,
  const CVec3d &dir,
  const CVec3d &org, // center of the circle
  const CVec3d &normal, // normal of the circle
  double rad);
#endif
