/*
 * Copyright (c) 2019 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */

#include "delfem2/geo_nearest3.h"

#include <cstdlib>
#include <cmath>
#include <stack>

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

// ------------------------------------------

DFM2_INLINE void delfem2::Nearest_LineSeg3_Point3(
  double nearest_position[3],
  const double point_position[3], // point
  const double lineseg_s[3], // source
  const double lineseg_e[3]) // end
{
  const double d[3] = {
    lineseg_e[0] - lineseg_s[0],
    lineseg_e[1] - lineseg_s[1],
    lineseg_e[2] - lineseg_s[2]};
  double t = 0.5;
  if (Dot3(d, d) > 1.0e-20) {
    const double ps[3] = {
      lineseg_s[0] - point_position[0],
      lineseg_s[1] - point_position[1],
      lineseg_s[2] - point_position[2]};
    double a = Dot3(d, d);
    double b = Dot3(d, ps);
    t = -b / a;
    if (t < 0) t = 0;
    if (t > 1) t = 1;
  }
  nearest_position[0] = lineseg_s[0] + t * d[0];
  nearest_position[1] = lineseg_s[1] + t * d[1];
  nearest_position[2] = lineseg_s[2] + t * d[2];
}

// ==============================

DFM2_INLINE void delfem2::Nearest_Triangle3_Point3(
  double nearest_position[3],
  double &r0,
  double &r1,
  const double ps[3], // origin point
  const double q0[3],
  const double q1[3],
  const double q2[3]) {
  namespace lcl = delfem2::geo_nearest3;
  double area, n012[3];
  UnitNormalAreaTri3(n012, area, q0, q1, q2);
  const double pe[3] = {ps[0] + n012[0], ps[1] + n012[1], ps[2] + n012[2]};
  const double v012 = lcl::Volume_Tet3(ps, q0, q1, q2);
  if (fabs(v012) > 1.0e-10) {
    const double sign = (v012 > 0) ? +1 : -1;
    const double v0 = lcl::Volume_Tet3(ps, q1, q2, pe) * sign;
    const double v1 = lcl::Volume_Tet3(ps, q2, q0, pe) * sign;
    const double v2 = lcl::Volume_Tet3(ps, q0, q1, pe) * sign;
    assert(fabs(v0 + v1 + v2) > 1.0e-10);
    double inv_v012 = 1.0 / (v0 + v1 + v2);
    r0 = v0 * inv_v012;
    r1 = v1 * inv_v012;
    const double r2 = (1.0 - r0 - r1);
    const double tol = 1.0e-4;
    if (r0 > -tol && r1 > -tol && r2 > -tol) {
      nearest_position[0] = q0[0] * r0 + q1[0] * r1 + q2[0] * r2;
      nearest_position[1] = q0[1] * r0 + q1[1] * r1 + q2[1] * r2;
      nearest_position[2] = q0[2] * r0 + q1[2] * r1 + q2[2] * r2;
      return;
    }
  }
  double r12[3];
  Nearest_LineSeg3_Point3(r12, ps, q1, q2);
  double r20[3];
  Nearest_LineSeg3_Point3(r20, ps, q2, q0);
  double r01[3];
  Nearest_LineSeg3_Point3(r01, ps, q0, q1);
  const double d12 = Distance3(r12, ps);
  const double d20 = Distance3(r20, ps);
  const double d01 = Distance3(r01, ps);
  if (d12 < d20) {
    if (d12 < d01) { // 12 is the smallest
      nearest_position[0] = r12[0];
      nearest_position[1] = r12[1];
      nearest_position[2] = r12[2];
      r0 = 0;
      r1 = Distance3(nearest_position, q2) / Distance3(q1, q2);
      return;
    }
  } else {
    if (d20 < d01) { // d20 is the smallest
      nearest_position[0] = r20[0];
      nearest_position[1] = r20[1];
      nearest_position[2] = r20[2];
      r0 = Distance3(nearest_position, q2) / Distance3(q0, q2);
      r1 = 0;
      return;
    }
  }
  nearest_position[0] = r01[0];
  nearest_position[1] = r01[1];
  nearest_position[2] = r01[2];
  r0 = Distance3(nearest_position, q1) / Distance3(q0, q1);
  r1 = 1 - r0;
}

// ---------------------------------------------------------------------------

template<class VEC, typename T>
VEC delfem2::Nearest_Line3_Point3(
  const VEC &point, // point
  const VEC &line_src, // source
  const VEC &line_dir) // direction
{
  assert(line_dir.dot(line_dir) > 1.0e-20);
  const VEC ps = line_src - point;
  const T a = line_dir.dot(line_dir);
  const T b = line_dir.dot(line_src - point);
  const T t = -b / a;
  return line_src + t * line_dir;
}
#ifdef DFM2_STATIC_LIBRARY
template delfem2::CVec3d delfem2::Nearest_Line3_Point3<delfem2::CVec3d, double>(
  const CVec3d &p,
  const CVec3d &s,
  const CVec3d &d);
template delfem2::CVec3f delfem2::Nearest_Line3_Point3<delfem2::CVec3f, float>(
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

// ---------------------------------------

template<typename T>
delfem2::CVec3<T> delfem2::Nearest_Origin3_LineSeg3(
  const CVec3<T> &s, // start
  const CVec3<T> &e) // end
{
  CVec3<T> d = e - s;
  double a = Dot(d, d);
  if (a < 1.0e-20) { return (s + e) * 0.5; }
  double b = Dot(d, s);
  double t = -b / a;
  if (t < 0) t = 0;
  if (t > 1) t = 1;
  return s + t * d;
}

// ----------------------------------------

// r0==0 -> p0==org
// r0==1 -> p1==org
template<typename T>
delfem2::CVec3<T> delfem2::Nearest_Origin3_LineSeg3(
  double &r0,
  const CVec3<T> &p0, // start
  const CVec3<T> &p1) // end
{
  CVec3<T> d = p1 - p0;
  double a = Dot(d, d);
  if (a < 1.0e-20) {
    r0 = 0.5;
    return (p0 + p1) * 0.5;
  }
  double b = Dot(d, p0);
  r0 = -b / a;
  if (r0 < 0) r0 = 0;
  if (r0 > 1) r0 = 1;
  return (1.0 - r0) * p0 + r0 * p1;
}

// ---------------------------------------

template<typename T>
delfem2::CVec3<T> delfem2::Nearest_LineSeg3_Point3(
  const CVec3<T> &p, // point
  const CVec3<T> &s, // start
  const CVec3<T> &e) // end
{
  CVec3<T> d = e - s;
  if (Dot(d, d) < 1.0e-20) {
    return (s + e) * static_cast<T>(0.5);
  }
  const CVec3<T> ps = s - p;
  T a = Dot(d, d);
  T b = Dot(d, s - p);
  T t = -b / a;
  if (t < 0) t = 0;
  if (t > 1) t = 1;
  return s + t * d;
}
#ifdef DFM2_STATIC_LIBRARY
template delfem2::CVec3d delfem2::Nearest_LineSeg3_Point3(
  const CVec3d &p,
  const CVec3d &s,
  const CVec3d &e);
template delfem2::CVec3f delfem2::Nearest_LineSeg3_Point3(
  const CVec3f &p,
  const CVec3f &s,
  const CVec3f &e);
#endif

// --------------------------------------

template<typename T>
delfem2::CVec3<T> delfem2::Nearest_LineSeg3_Point3(
  T &t,
  const CVec3<T> &p, // point
  const CVec3<T> &s, // source
  const CVec3<T> &e) // end
{
  CVec3<T> d = e - s;
  if (Dot(d, d) < 1.0e-20) {
    t = static_cast<T>(0.5);
    return (1 - t) * s + t * e;
  }
  const CVec3<T> ps = s - p;
  T a = Dot(d, d);
  T b = Dot(d, s - p);
  t = -b / a;
  if (t < 0) t = 0;
  if (t > 1) t = 1;
  return s + t * d;
}
#ifdef DFM2_STATIC_LIBRARY
template delfem2::CVec3d delfem2::Nearest_LineSeg3_Point3(
  double &t,
  const CVec3d &p, // point
  const CVec3d &s, // source
  const CVec3d &e); // end
template delfem2::CVec3f delfem2::Nearest_LineSeg3_Point3(
  float &t,
  const CVec3f &p, // point
  const CVec3f &s, // source
  const CVec3f &e); // end
#endif

// ---------------------------------------------

template<class VEC, typename T>
void delfem2::Nearest_LineSeg3_Line3(
  VEC &nearest_lineseg,
  VEC &nearest_line,
  const VEC &lineseg_start,
  const VEC &lineseg_end,
  const VEC &line_origin,
  const VEC &line_direction) {
  T D0, Dta0, Dtb0;
  VEC Da0, Db0;
  Nearest_Line3_Line3(
    D0, Da0, Db0, Dta0, Dtb0,
    lineseg_start, lineseg_end - lineseg_start, line_origin, line_direction);
  if (abs(D0) < 1.0e-10) { // pararell
    nearest_lineseg = (lineseg_start + lineseg_end) * static_cast<T>(0.5);
    nearest_line = ::delfem2::Nearest_Line3_Point3<VEC, T>(
      nearest_lineseg, line_origin, line_direction);
    return;
  }
  const T ta = Dta0 / D0;
  if (ta > 0 && ta < 1) { // nearst point is inside the segment
    nearest_lineseg = Da0 / D0;
    nearest_line = Db0 / D0;
    return;
  }
  //
  const VEC p1 = Nearest_Line3_Point3<VEC, T>(lineseg_start, line_origin, line_direction);
  const VEC p2 = Nearest_Line3_Point3<VEC, T>(lineseg_end, line_origin, line_direction);
  const T Dist1 = (p1 - lineseg_start).norm();
  const T Dist2 = (p2 - lineseg_end).norm();
  if (Dist1 < Dist2) {
    nearest_lineseg = lineseg_start;
    nearest_line = p1;
    return;
  }
  nearest_lineseg = lineseg_end;
  nearest_line = p2;
}
#ifdef DFM2_STATIC_LIBRARY
template void delfem2::Nearest_LineSeg3_Line3<delfem2::CVec3f, float>(
  CVec3f &, CVec3f &,
  const CVec3f &, const CVec3f &,
  const CVec3f &, const CVec3f &);
template void delfem2::Nearest_LineSeg3_Line3<delfem2::CVec3d, double>(
  CVec3d &, CVec3d &,
  const CVec3d &, const CVec3d &,
  const CVec3d &, const CVec3d &);
#endif

// ---------------------------------------------

DFM2_INLINE double delfem2::Nearest_LineSeg_LineSeg_CCD_Iteration(
  double p[3],
  const CVec3d &p0s,
  const CVec3d &p0e,
  const CVec3d &p1s,
  const CVec3d &p1e,
  const CVec3d &q0s,
  const CVec3d &q0e,
  const CVec3d &q1s,
  const CVec3d &q1e,
  unsigned int nitr) {
  namespace lcl = delfem2::geo_nearest3;
  CVec3d v0;
  for (unsigned int itr = 0; itr < nitr; ++itr) {
    const double s0 = p[0], t0 = p[1], u0 = p[2];
    v0 =
      +((1 - s0) * (1 - u0)) * p0s + ((1 - s0) * u0) * p0e + (s0 * (1 - u0)) * p1s + (s0 * u0) * p1e
      - ((1 - t0) * (1 - u0)) * q0s - ((1 - t0) * u0) * q0e - (t0 * (1 - u0)) * q1s - (t0 * u0) * q1e;
//    std::cout << "   " << itr << " " << v0.Length() << "  " << p[0] << " " << p[1] << " " << p[2] << std::endl;
    const CVec3d ds = -(1 - u0) * p0s - u0 * p0e + (1 - u0) * p1s + u0 * p1e;
    const CVec3d dt = +(1 - u0) * q0s + u0 * q0e - (1 - u0) * q1s - u0 * q1e;
    const CVec3d du =
      -(1 - s0) * p0s + (1 - s0) * p0e - s0 * p1s + s0 * p1e
        + (1 - t0) * q0s - (1 - t0) * q0e + t0 * q1s - t0 * q1e;
    const CVec3d dsu = +p0s - p0e - p1s + p1e;
    const CVec3d dtu = -q0s + q0e + q1s - q1e;
    double R[3] = {v0.dot(ds), v0.dot(dt), v0.dot(du)};
    double A[9] = {
      ds.dot(ds), ds.dot(dt), ds.dot(du) + v0.dot(dsu),
      dt.dot(ds), dt.dot(dt), dt.dot(du) * v0.dot(dtu),
      du.dot(ds) + v0.dot(dsu), du.dot(dt) + v0.dot(dtu), du.dot(du)};
    {
      double eps = (A[0] + A[4] + A[8]) * 1.0e-10 + 1.0e-20;
      A[0] += eps;
      A[4] += eps;
      A[8] += eps;
    }
    double Ainv[9];
    lcl::MyInverse_Mat3(Ainv, A);
    double D[3];
    lcl::MyMatVec3(D, Ainv, R);
    p[0] -= D[0];
    p[1] -= D[1];
    p[2] -= D[2];
    if (p[0] < 0) { p[0] = 0.0; } else if (p[0] > 1) { p[0] = 1.0; }
    if (p[1] < 0) { p[1] = 0.0; } else if (p[1] > 1) { p[1] = 1.0; }
    if (p[2] < 0) { p[2] = 0.0; } else if (p[2] > 1) { p[2] = 1.0; }
  }
  return v0.norm();
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
delfem2::CVec3<T> delfem2::Nearest_Orgin3_PlaneTri3(
  T &r0,
  T &r1,
  const CVec3<T> &q0,
  const CVec3<T> &q1,
  const CVec3<T> &q2) {
  namespace lcl = delfem2::geo_nearest3;
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

// -----------------------------

template<typename T>
delfem2::CVec3<T> delfem2::Nearest_Origin3_Tri3(
  T &r0,
  T &r1,
  const CVec3<T> &q0,
  const CVec3<T> &q1,
  const CVec3<T> &q2) {

  if (((q1 - q0) ^ (q2 - q0)).norm() > 1.0e-10) {
    CVec3<T> p012 = Nearest_Orgin3_PlaneTri3(r0, r1, q0, q1, q2);
    if (r0 > 0 && r1 > 0 && (1 - r0 - r1) > 0) { return p012; }
  }
  CVec3<T> p_min = q0;
  T d_min = q0.norm();
  r0 = 1;
  r1 = 0;
  {
    T s2;
    CVec3<T> p12 = Nearest_Origin3_LineSeg3(s2, q1, q2);
    const T d12 = p12.norm();
    if (d12 < d_min) {
      d_min = d12;
      p_min = p12;
      r1 = 1 - s2;
      r0 = 0;
    }
  }
  {
    T s0;
    CVec3<T> p20 = Nearest_Origin3_LineSeg3(s0, q2, q0);
    const T d20 = p20.norm();
    if (d20 < d_min) {
      d_min = d20;
      p_min = p20;
      r1 = 0;
      r0 = s0;
    }
  }
  {
    T s1;
    CVec3<T> p01 = Nearest_Origin3_LineSeg3(s1, q0, q1);
    const T d01 = p01.norm();
    if (d01 < d_min) {
      d_min = d01;
      p_min = p01;
      r0 = 1 - s1;
      r1 = s1;
    }
  }
  return p_min;
}
#ifdef DFM2_STATIC_LIBRARY
template delfem2::CVec3d delfem2::Nearest_Origin3_Tri3(
  double &r0, double &r1,
  const CVec3d &q0, const CVec3d &q1, const CVec3d &q2);
#endif

// -------------------------------------------

template<typename T>
delfem2::CVec3<T> delfem2::Nearst_Origin3_Quad3(
  double &s0,
  double &s1,
  const CVec3<T> &q0,
  const CVec3<T> &q1,
  const CVec3<T> &q2,
  const CVec3<T> &q3) {
  double dist_min = -1;
  CVec3<T> q_min;
  for (int ip = 0; ip < 5; ++ip) {
    double t0 = 0, t1 = 0;
    if (ip == 0) {
      t0 = 0.0;
      t1 = 0.0;
    } else if (ip == 1) {
      t0 = 1.0;
      t1 = 0.0;
    } else if (ip == 2) {
      t0 = 1.0;
      t1 = 1.0;
    } else if (ip == 3) {
      t0 = 0.0;
      t1 = 1.0;
    } else if (ip == 4) {
      t0 = 0.5;
      t1 = 0.5;
    }
    CVec3<T> q;
    for (int itr = 0; itr < 4; ++itr) {
      CVec3<T> pq = (1 - t0) * (1 - t1) * q0 + t0 * (1 - t1) * q1 + t0 * t1 * q2 + (1 - t0) * t1 * q3;
      CVec3<T> dqt0 = -(1 - t1) * q0 + (1 - t1) * q1 + t1 * q2 - t1 * q3;
      CVec3<T> dqt1 = -(1 - t0) * q0 - t0 * q1 + t0 * q2 + (1 - t0) * q3;
      CVec3<T> ddqt0t1 = q0 - q1 + q2 - q3;
      double f0 = -dqt0 * pq;
      double f1 = -dqt1 * pq;
      double A00 = dqt0 * dqt0;
      double A11 = dqt1 * dqt1;
      double A01 = dqt1 * dqt0 + ddqt0t1 * pq;
      double det = A00 * A11 - A01 * A01;
      double detinv = 1.0 / det;
      double B00 = +A11 * detinv;
      double B11 = +A00 * detinv;
      double B01 = -A01 * detinv;
      double d0 = B00 * f0 + B01 * f1;
      double d1 = B01 * f0 + B11 * f1;
      t0 += d0;
      t1 += d1;
    }
    double tol = 1.0e-4;
    if (t0 > -tol && t0 < 1.0 + tol && t1 > -tol && t1 < 1.0 + tol) {
      double d0 = q.norm();
      if (dist_min < 0 || d0 < dist_min) {
        dist_min = d0;
        s0 = t0;
        s1 = t1;
        q_min = q;
      }
    }
  }
  if (dist_min > 0) { return q_min; }
  //
  const CVec3<T> q01 = Nearest_Origin3_LineSeg3(q0, q1);
  const double d01 = q01.norm();
  if (dist_min < 0 || d01 < dist_min) {
    dist_min = d01;
    s0 = Distance(q01, q0) / Distance(q0, q1);
    s1 = 0.0;
    q_min = q01;
  }
  //
  CVec3<T> q12 = Nearest_Origin3_LineSeg3(q1, q2);
  const double d12 = q12.norm();
  if (dist_min < 0 || d12 < dist_min) {
    dist_min = d12;
    s0 = 1.0;
    s1 = Distance(q12, q1) / Distance(q1, q2);
    q_min = q12;
  }
  //
  CVec3<T> q23 = Nearest_Origin3_LineSeg3(q2, q3);
  const double d23 = q23.norm();
  if (dist_min < 0 || d23 < dist_min) {
    dist_min = d23;
    s0 = Distance(q23, q3) / Distance(q2, q3);
    s1 = 1.0;
    q_min = q23;
  }
  //
  CVec3<T> q30 = Nearest_Origin3_LineSeg3(q3, q0);
  const double d30 = q30.norm();
  if (dist_min < 0 || d30 < dist_min) {
    dist_min = d30;
    s0 = 0.0;
    s1 = Distance(q30, q0) / Distance(q3, q0);
    q_min = q30;
  }
  return q_min;
}

/*
CVector3 nearest_Origin_Tet
(double& r0, double& r1, double& r2,
 const CVector3& q0,
 const CVector3& q1,
 const CVector3& q2,
 const CVector3& q3)
{
  CVector3 p_min = q0;
  {
    bool res = barycentricCoord_Origin_Tet(r0, r1, r2, q0, q1, q2, q3);
    p_min = r0*q0 + r1*q1 + r2*q2 + (1-r0-r1-r2)*q3;
    if( r0>0 && r1>0 && r2>0 && (1-r0-r1-r2)>0 ){ return p_min; }
  }
  ////////////////////////
  double r3;
  { // face123
    r0 = 0;
    p_min = nearest_Orgin_PlaneTri(r1,r2,r3, q1,q2,q3);
    if( r1>0 && r2>0 && r3>0 ){ return p_min; }
  }
  { // face230
    r1 = 0;
    p_min = nearest_Orgin_PlaneTri(r2,r3,r0, q2,q3,q0);
    if( r2>0 && r3>0 && r0>0 ){ return p_min; }
  }
  { // face301
    r2 = 0;
    p_min = nearest_Orgin_PlaneTri(r3,r0,r1, q3,q0,q1);
    if( r3>0 && r0>0 && r1>0 ){ return p_min; }
  }
  { // face012
    r3 = 0;
    p_min = nearest_Orgin_PlaneTri(r0,r1,r2, q0,q1,q2);
    if( r0>0 && r1>0 && r2>0 ){ return p_min; }
  }
  ////////////////////////
  double d_min = q0.Length();
  double s0,s1,s2,s3;
  { // edge01
    CVector3 p01 = nearest_Origin_LineSeg(s0,s1,q0,q1);
    double d01 = p01.Length();
    if( d01<d_min ){ d_min=d01; p_min=p01; r0=s0; r1=s1; r2=0; r3=0; }
  }
  { // edge02
    CVector3 p02 = nearest_Origin_LineSeg(s0,s2,q0,q2);
    double d02 = p02.Length();
    if( d02<d_min ){ d_min=d02; p_min=p02; r0=s0; r1=0; r2=s2; r3=0; }
  }
  { // edge03
    CVector3 p03 = nearest_Origin_LineSeg(s0,s3,q0,q3);
    double d03 = p03.Length();
    if( d03<d_min ){ d_min=d03; p_min=p03; r0=s0; r1=0; r2=0; r3=s3; }
  }
  { // edge12
    CVector3 p12 = nearest_Origin_LineSeg(s1,s2,q1,q2);
    double d12 = p12.Length();
    if( d12<d_min ){ d_min=d12; p_min=p12; r0=0; r1=s1; r2=s2; r3=0; }
  }
  { // edge13
    CVector3 p13 = nearest_Origin_LineSeg(s1,s3,q1,q3);
    double d13 = p13.Length();
    if( d13<d_min ){ d_min=d13; p_min=p13; r0=0; r1=s1; r2=0; r3=s3; }
  }
  { // edge23
    CVector3 p23 = nearest_Origin_LineSeg(s2,s3,q2,q3);
    double d23 = p23.Length();
    if( d23<d_min ){ d_min=d23; p_min=p23; r0=0; r1=0; r2=s2; r3=s3; }
  }
  return p_min;
}
 */

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
  GetVertical2Vector(normal, ex, ey);
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
