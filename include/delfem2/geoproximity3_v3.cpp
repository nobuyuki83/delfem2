/*
 * Copyright (c) 2019 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */

#include "delfem2/geoproximity3_v3.h"
#include "delfem2/geo_nearest3.h"

#include <cstdlib>
#include <cmath>
#include <stack>

#ifndef M_PI
#  define M_PI 3.14159265358979323846
#endif

// =====================================
// below: unexposed 

namespace delfem2::proximity3 {

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

template<typename REAL>
bool delfem2::IntersectRay_Tri3(
  REAL &r0,
  REAL &r1,
  const CVec3<REAL> &org,
  const CVec3<REAL> &dir,
  const CVec3<REAL> &p0,
  const CVec3<REAL> &p1,
  const CVec3<REAL> &p2,
  REAL eps) {
  namespace lcl = delfem2::proximity3;
  const REAL v0 = lcl::Volume_Tet(p1, p2, org, org + dir);
  const REAL v1 = lcl::Volume_Tet(p2, p0, org, org + dir);
  const REAL v2 = lcl::Volume_Tet(p0, p1, org, org + dir);
  const REAL vt = v0 + v1 + v2;
  r0 = v0 / vt;
  r1 = v1 / vt;
  const REAL r2 = v2 / vt;
  return (r0 >= -eps && r1 >= -eps && r2 >= -eps);
}
#ifdef DFM2_STATIC_LIBRARY
template bool delfem2::IntersectRay_Tri3(
  double &r0, double &r1,
  const CVec3d &org, const CVec3d &dir,
  const CVec3d &p0, const CVec3d &p1, const CVec3d &p2,
  double eps);
template bool delfem2::IntersectRay_Tri3(
  float &r0, float &r1,
  const CVec3f &org, const CVec3f &dir,
  const CVec3f &p0, const CVec3f &p1, const CVec3f &p2,
  float eps);
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

// distance VF
template<typename T>
double delfem2::DistanceFaceVertex
  (const CVec3<T> &p0, const CVec3<T> &p1, const CVec3<T> &p2,
   const CVec3<T> &p3,
   double &w0, double &w1) {
  CVec3<T> v20 = p0 - p2;
  CVec3<T> v21 = p1 - p2;
  double t0 = Dot(v20, v20);
  double t1 = Dot(v21, v21);
  double t2 = Dot(v20, v21);
  double t3 = Dot(v20, p3 - p2);
  double t4 = Dot(v21, p3 - p2);
  double det = t0 * t1 - t2 * t2;
  double invdet = 1.0 / det;
  w0 = (+t1 * t3 - t2 * t4) * invdet;
  w1 = (-t2 * t3 + t0 * t4) * invdet;
  const double w2 = 1 - w0 - w1;
  CVec3<T> pw = w0 * p0 + w1 * p1 + w2 * p2;
  return (pw - p3).norm();
}
#ifdef DFM2_STATIC_LIBRARY
template double delfem2::DistanceFaceVertex(
  const CVec3d &p0, const CVec3d &p1,
  const CVec3d &p2, const CVec3d &p3,
  double &w0, double &w1);
#endif

// ---------------------------------------------------------

//　distance EE
template<typename T>
double delfem2::DistanceEdgeEdge(
  const CVec3<T> &p0,
  const CVec3<T> &p1,
  const CVec3<T> &q0,
  const CVec3<T> &q1,
  double &ratio_p,
  double &ratio_q) {
  const CVec3<T> &vp = p1 - p0;
  const CVec3<T> &vq = q1 - q0;
  if (Cross(vp, vq).norm() < 1.0e-10) { // handling parallel edge
    CVec3<T> pq0 = p0 - q0;
    CVec3<T> nvp = vp;
    nvp.normalize();
    CVec3<T> vert = pq0 - Dot(pq0, nvp) * nvp;
    double dist = vert.norm();
    double lp0 = Dot(p0, nvp);
    double lp1 = Dot(p1, nvp);
    double lq0 = Dot(q0, nvp);
    double lq1 = Dot(q1, nvp);
    double p_min = (lp0 < lp1) ? lp0 : lp1;
    double p_max = (lp0 > lp1) ? lp0 : lp1;
    double q_min = (lq0 < lq1) ? lq0 : lq1;
    double q_max = (lq0 > lq1) ? lq0 : lq1;
    double lm;
    if (p_max < q_min) { lm = (p_max + q_min) * 0.5; }
    else if (q_max < p_min) { lm = (q_max + p_min) * 0.5; }
    else if (p_max < q_max) { lm = (p_max + q_min) * 0.5; }
    else { lm = (q_max + p_min) * 0.5; }
    ratio_p = (lm - lp0) / (lp1 - lp0);
    ratio_q = (lm - lq0) / (lq1 - lq0);
    return dist;
  }
  double t0 = Dot(vp, vp);
  double t1 = Dot(vq, vq);
  double t2 = Dot(vp, vq);
  double t3 = Dot(vp, q0 - p0);
  double t4 = Dot(vq, q0 - p0);
  double det = t0 * t1 - t2 * t2;
  double invdet = 1.0 / det;
  ratio_p = (+t1 * t3 - t2 * t4) * invdet;
  ratio_q = (+t2 * t3 - t0 * t4) * invdet;
  CVec3<T> pc = p0 + ratio_p * vp;
  CVec3<T> qc = q0 + ratio_q * vq;
  return (pc - qc).norm();
}
#ifdef DFM2_STATIC_LIBRARY
template double delfem2::DistanceEdgeEdge(
  const CVec3d &p0, const CVec3d &p1,
  const CVec3d &q0, const CVec3d &q1,
  double &ratio_p, double &ratio_q);
#endif

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
  namespace lcl = delfem2::proximity3;
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

template<typename T>
void delfem2::iteration_intersection_Line_Quad
  (double &t0, double &t1,
   const CVec3<T> &src,
   const CVec3<T> &u,
   const CVec3<T> &v,
   const CVec3<T> &q0,
   const CVec3<T> &q1,
   const CVec3<T> &q2,
   const CVec3<T> &q3) {
  CVec3<T> q = (1 - t0) * (1 - t1) * q0 + t0 * (1 - t1) * q1 + t0 * t1 * q2 + (1 - t0) * t1 * q3;
  CVec3<T> pq = q - src;
  CVec3<T> dqt0 = -(1 - t1) * q0 + (1 - t1) * q1 + t1 * q2 - t1 * q3;
  CVec3<T> dqt1 = -(1 - t0) * q0 - t0 * q1 + t0 * q2 + (1 - t0) * q3;
  CVec3<T> ddqt0t1 = q0 - q1 + q2 - q3;
  double f0 = -u * pq;
  double f1 = -v * pq;
  double A00 = u * dqt0;
  double A01 = u * dqt1;
  double A10 = v * dqt0;
  double A11 = v * dqt1;
  double det = A00 * A11 - A01 * A10;
  double detinv = 1.0 / det;
  double B00 = +A11 * detinv;
  double B01 = -A01 * detinv;
  double B10 = -A10 * detinv;
  double B11 = +A00 * detinv;
  double d0 = B00 * f0 + B01 * f1;
  double d1 = B10 * f0 + B11 * f1;
  t0 += d0;
  t1 += d1;
}

template<typename T>
bool delfem2::intersection_Point_Quad(
  CVec3<T> &psec, double &s0, double &s1,
  const CVec3<T> &src, const CVec3<T> &dir,
  const CVec3<T> &q0, const CVec3<T> &q1, const CVec3<T> &q2, const CVec3<T> &q3) {
  CVec3<T> u, v;
  GetVertical2Vector(dir, u, v);
  //
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
    for (int itr = 0; itr < 4; ++itr) {
      iteration_intersection_Line_Quad(t0, t1,
                                       src, u, v,
                                       q0, q1, q2, q3);
    }
    CVec3<T> q = (1 - t0) * (1 - t1) * q0 + t0 * (1 - t1) * q1 + t0 * t1 * q2 + (1 - t0) * t1 * q3;
    double tol = 1.0e-4;
    //    std::cout << t0 << " " << t1 << std::endl;
    if (t0 > -tol && t0 < 1.0 + tol && t1 > -tol && t1 < 1.0 + tol) {
      double d0 = (q - src).Length();
      if (dist_min < 0 || d0 < dist_min) {
        dist_min = d0;
        s0 = t0;
        s1 = t1;
        q_min = q;
      }
    }
  }
  //  std::cout << dist_min << std::endl;
  if (dist_min > 0) {
    psec = q_min;
    return true;
  }
  return false;
}

// EEの距離が所定の距離以下にあるかどうか
template<typename T>
bool delfem2::IsContact_EE_Proximity(
  int ino0,
  int ino1,
  int jno0,
  int jno1,
  const CVec3<T> &p0,
  const CVec3<T> &p1,
  const CVec3<T> &q0,
  const CVec3<T> &q1,
  const double delta) {
  if (ino0 == jno0 || ino0 == jno1 || ino1 == jno0 || ino1 == jno1) return false;
  if (q0.p[0] + delta < p0.p[0] && q0.p[0] + delta < p1.p[0] && q1.p[0] + delta < p0.p[0]
    && q1.p[0] + delta < p1.p[0])
    return false;
  if (q0.p[0] - delta > p0.p[0] && q0.p[0] - delta > p1.p[0] && q1.p[0] - delta > p0.p[0]
    && q1.p[0] - delta > p1.p[0])
    return false;
  if (q0.p[1] + delta < p0.p[1] && q0.p[1] + delta < p1.p[1] && q1.p[1] + delta < p0.p[1]
    && q1.p[1] + delta < p1.p[1])
    return false;
  if (q0.p[1] - delta > p0.p[1] && q0.p[1] - delta > p1.p[1] && q1.p[1] - delta > p0.p[1]
    && q1.p[1] - delta > p1.p[1])
    return false;
  if (q0.p[2] + delta < p0.p[2] && q0.p[2] + delta < p1.p[2] && q1.p[2] + delta < p0.p[2]
    && q1.p[2] + delta < p1.p[2])
    return false;
  if (q0.p[2] - delta > p0.p[2] && q0.p[2] - delta > p1.p[2] && q1.p[2] - delta > p0.p[2]
    && q1.p[2] - delta > p1.p[2])
    return false;
  double ratio_p, ratio_q;
  double dist = DistanceEdgeEdge(p0, p1, q0, q1, ratio_p, ratio_q);
  if (dist > delta) return false;
  if (ratio_p < 0) return false;
  if (ratio_p > 1) return false;
  if (ratio_q < 0) return false;
  if (ratio_q > 1) return false;
  const CVec3<T> &pm = (1 - ratio_p) * p0 + ratio_p * p1;
  const CVec3<T> &qm = (1 - ratio_q) * q0 + ratio_q * q1;
  return (pm - qm).norm() <= delta;
}
#ifdef DFM2_STATIC_LIBRARY
template bool delfem2::IsContact_EE_Proximity(
  int ino0, int ino1,
  int jno0, int jno1,
  const CVec3d &p0, const CVec3d &p1,
  const CVec3d &q0, const CVec3d &q1,
  const double delta);
#endif

// ------------------

// compute time where four points gets coplaner
template<typename T>
bool delfem2::FindCoplanerInterp(
  double &r,
  const CVec3<T> &s0,
  const CVec3<T> &s1,
  const CVec3<T> &s2,
  const CVec3<T> &s3,
  const CVec3<T> &e0,
  const CVec3<T> &e1,
  const CVec3<T> &e2,
  const CVec3<T> &e3) {
  namespace lcl = delfem2::proximity3;
  const CVec3<T> x1 = s1 - s0;
  const CVec3<T> x2 = s2 - s0;
  const CVec3<T> x3 = s3 - s0;
  const CVec3<T> v1 = e1 - e0 - x1;
  const CVec3<T> v2 = e2 - e0 - x2;
  const CVec3<T> v3 = e3 - e0 - x3;
  // compute coefficient for cubic function
  const T k0 = ScalarTripleProduct(x3, x1, x2);
  const T k1 = ScalarTripleProduct(v3, x1, x2) + ScalarTripleProduct(x3, v1, x2) + ScalarTripleProduct(x3, x1, v2);
  const T k2 = ScalarTripleProduct(v3, v1, x2) + ScalarTripleProduct(v3, x1, v2) + ScalarTripleProduct(x3, v1, v2);
  const T k3 = ScalarTripleProduct(v3, v1, v2);
  // cubic funciton is f(x) = k0 + k1*x + k2*x^2 + k3*x^3
  const T r0 = +0.0;
  const T r1 = +1.0;
  const T f0 = lcl::EvaluateCubic(r0, k0, k1, k2, k3);
  const T f1 = lcl::EvaluateCubic(r1, k0, k1, k2, k3);
  if (f0 * f1 <= 0) {
    r = lcl::FindRootCubic_Bisect(r0, r1, f0, f1, k0, k1, k2, k3);
    return true;
  }
  if (fabs(k3) > 1.0e-30) { // cubic function
    const double det = k2 * k2 - 3 * k1 * k3; // if det > 0, the cubic function takes extreme value
    if (det < 0) { return false; } // monotonus function
    //
    const double r3 = (-k2 - sqrt(det)) / (3 * k3); // smaller extreme value
    if (r3 > 0 && r3 < 1) {
      const double f3 = lcl::EvaluateCubic(r3, k0, k1, k2, k3);
      if (f3 == 0) {
        r = r3;
        return true;
      }
      if (f0 * f3 < 0) {
        r = lcl::FindRootCubic_Bisect(r0, r3, f0, f3, k0, k1, k2, k3);
        return true;
      }
    }
    const double r4 = (-k2 + sqrt(det)) / (3 * k3); // larger extreme value
    if (r4 > 0 && r4 < 1) {
      const double f4 = lcl::EvaluateCubic(r4, k0, k1, k2, k3);
      if (f4 == 0) {
        r = r4;
        return true;
      }
      if (f0 * f4 < 0) {
        r = lcl::FindRootCubic_Bisect(r0, r4, f0, f4, k0, k1, k2, k3);
        return true;
      }
    }
    return false;
  }
  //
  if (fabs(k2) > 1.0e-30) { // quadric function
    const double r2 = -k1 / (2 * k2); // extreme valuse
    if (r2 > 0 && r2 < 1) {
      const double f2 = lcl::EvaluateCubic(r2, k0, k1, k2, k3);
      if (f0 * f2 < 0) {
        r = lcl::FindRootCubic_Bisect(r0, r2, f0, f2, k0, k1, k2, k3);
        return true;
      }
    }
    return false;
  }
  return false;
}

// CCDのFVで接触する要素を検出
template<typename T>
bool delfem2::IsContact_FV_CCD2(
  [[maybe_unused]] int ino0,
  [[maybe_unused]] int ino1,
  [[maybe_unused]] int ino2,
  [[maybe_unused]] int ino3,
  const CVec3<T> &p0,
  const CVec3<T> &p1,
  const CVec3<T> &p2,
  const CVec3<T> &p3,
  const CVec3<T> &q0,
  const CVec3<T> &q1,
  const CVec3<T> &q2,
  const CVec3<T> &q3) {
  { // CSAT
    CVec3<T> n = Cross(p1 - p0, p2 - p0);
    double t0 = Dot(p0 - p3, n);
    double t1 = Dot(q0 - q3, n);
    double t2 = Dot(q1 - q3, n);
    double t3 = Dot(q2 - q3, n);
    if (t0 * t1 > 0 && t0 * t2 > 0 && t0 * t3 > 0) { return false; }
  }
  double r0, r1;
  double dist = DistanceFaceVertex(p0, p1, p2, p3, r0, r1);
  {
    double vn0 = (p0 - q0).norm();
    double vn1 = (p1 - q1).norm();
    double vn2 = (p2 - q2).norm();
    double vn3 = (p3 - q3).norm();
    double vnt = (vn0 > vn1) ? vn0 : vn1;
    vnt = (vn2 > vnt) ? vn2 : vnt;
    double max_app = (vnt + vn3);
    const double r2 = 1 - r0 - r1;
    if (dist > max_app) return false;
    if (r0 < 0 || r0 > 1 || r1 < 0 || r1 > 1 || r2 < 0 || r2 > 1) {
      double dist01 = (Nearest_LineSeg3_Point3(p3, p0, p1) - p3).norm();
      double dist12 = (Nearest_LineSeg3_Point3(p3, p1, p2) - p3).norm();
      double dist20 = (Nearest_LineSeg3_Point3(p3, p2, p0) - p3).norm();
      if (dist01 > max_app && dist12 > max_app && dist20 > max_app) { return false; }
    }
  }
  double t;
  {
    bool res = FindCoplanerInterp(t,
                                  p0, p1, p2, p3, q0, q1, q2, q3);
    if (!res) return false;
    assert(t >= 0 && t <= 1);
  }
  CVec3<T> p0m = (1 - t) * p0 + t * q0;
  CVec3<T> p1m = (1 - t) * p1 + t * q1;
  CVec3<T> p2m = (1 - t) * p2 + t * q2;
  CVec3<T> p3m = (1 - t) * p3 + t * q3;
  double w0, w1;
  DistanceFaceVertex(p0m, p1m, p2m, p3m, w0, w1);
  double w2 = 1 - w0 - w1;
  if (w0 < 0 || w0 > 1) return false;
  if (w1 < 0 || w1 > 1) return false;
  if (w2 < 0 || w2 > 1) return false;
  return true;
}
#ifdef DFM2_STATIC_LIBRARY
template bool delfem2::IsContact_FV_CCD2(
  int ino0,
  int ino1,
  int ino2,
  int ino3,
  const CVec3d &p0,
  const CVec3d &p1,
  const CVec3d &p2,
  const CVec3d &p3,
  const CVec3d &q0,
  const CVec3d &q1,
  const CVec3d &q2,
  const CVec3d &q3);
#endif

// ----------------------------------

template<typename T>
bool delfem2::isIntersectTriPair(
  CVec3<T> &P0,
  CVec3<T> &P1,
  int itri,
  int jtri,
  const std::vector<unsigned int> &aTri,
  const std::vector<double> &aXYZ) {
  const unsigned int i0 = aTri[itri * 3 + 0];
  const unsigned int i1 = aTri[itri * 3 + 1];
  const unsigned int i2 = aTri[itri * 3 + 2];
  const unsigned int j0 = aTri[jtri * 3 + 0];
  const unsigned int j1 = aTri[jtri * 3 + 1];
  const unsigned int j2 = aTri[jtri * 3 + 2];
  if (i0 == j0 || i0 == j1 || i0 == j2) return false;
  if (i1 == j0 || i1 == j1 || i1 == j2) return false;
  if (i2 == j0 || i2 == j1 || i2 == j2) return false;
  const CVec3<T> p0(aXYZ[i0 * 3 + 0], aXYZ[i0 * 3 + 1], aXYZ[i0 * 3 + 2]);
  const CVec3<T> p1(aXYZ[i1 * 3 + 0], aXYZ[i1 * 3 + 1], aXYZ[i1 * 3 + 2]);
  const CVec3<T> p2(aXYZ[i2 * 3 + 0], aXYZ[i2 * 3 + 1], aXYZ[i2 * 3 + 2]);
  const CVec3<T> q0(aXYZ[j0 * 3 + 0], aXYZ[j0 * 3 + 1], aXYZ[j0 * 3 + 2]);
  const CVec3<T> q1(aXYZ[j1 * 3 + 0], aXYZ[j1 * 3 + 1], aXYZ[j1 * 3 + 2]);
  const CVec3<T> q2(aXYZ[j2 * 3 + 0], aXYZ[j2 * 3 + 1], aXYZ[j2 * 3 + 2]);
  const CVec3<T> np = Normal(p0, p1, p2);
  const CVec3<T> nq = Normal(q0, q1, q2);
  double dp0 = (p0 - q0).dot(nq);
  double dp1 = (p1 - q0).dot(nq);
  double dp2 = (p2 - q0).dot(nq);
  double dq0 = (q0 - p0).dot(np);
  double dq1 = (q1 - p0).dot(np);
  double dq2 = (q2 - p0).dot(np);
  if (((dp0 > 0) == (dp1 > 0)) && ((dp1 > 0) == (dp2 > 0))) return false;
  if (((dq0 > 0) == (dq1 > 0)) && ((dq1 > 0) == (dq2 > 0))) return false;
  const CVec3<T> p01 = (1.0 / (dp0 - dp1)) * (dp0 * p1 - dp1 * p0);
  const CVec3<T> p12 = (1.0 / (dp1 - dp2)) * (dp1 * p2 - dp2 * p1);
  const CVec3<T> p20 = (1.0 / (dp2 - dp0)) * (dp2 * p0 - dp0 * p2);
  const CVec3<T> q01 = (1.0 / (dq0 - dq1)) * (dq0 * q1 - dq1 * q0);
  const CVec3<T> q12 = (1.0 / (dq1 - dq2)) * (dq1 * q2 - dq2 * q1);
  const CVec3<T> q20 = (1.0 / (dq2 - dq0)) * (dq2 * q0 - dq0 * q2);
  const CVec3<T> vz = Cross(np, nq);
  CVec3<T> ps, pe;
  if (dp0 * dp1 > 0) {
    ps = p20;
    pe = p12;
  } else if (dp1 * dp2 > 0) {
    ps = p01;
    pe = p20;
  } else {
    ps = p12;
    pe = p01;
  }
  if (ps.dot(vz) > pe.dot(vz)) {
    CVec3<T> pt = ps;
    ps = pe;
    pe = pt;
  }
  double zps = ps.dot(vz);
  double zpe = pe.dot(vz);
  assert(zps <= zpe);
  //
  CVec3<T> qs, qe;
  if (dq0 * dq1 > 0) {
    qs = q20;
    qe = q12;
  } else if (dq1 * dq2 > 0) {
    qs = q01;
    qe = q20;
  } else {
    qs = q12;
    qe = q01;
  }
  if (qs.dot(vz) > qe.dot(vz)) {
    CVec3<T> qt = qs;
    qs = qe;
    qe = qt;
  }
  double zqs = qs.dot(vz);
  double zqe = qe.dot(vz);
  assert(zqs <= zqe);
  //
  if (zps > zqe || zqs > zpe) return false;
  CVec3<T> P[4];
  int icnt = 0;
  if (zps > zqs && zps < zqe) {
    P[icnt] = ps;
    icnt++;
  }
  if (zpe > zqs && zpe < zqe) {
    P[icnt] = pe;
    icnt++;
  }
  if (zqs > zps && zqs < zpe) {
    P[icnt] = qs;
    icnt++;
  }
  if (zqe > zps && zqe < zpe) {
    P[icnt] = qe;
    icnt++;
  }
  if (icnt != 2) return false;
  P0 = P[0];
  P1 = P[1];
  return true;
}
#ifdef DFM2_STATIC_LIBRARY
template bool delfem2::isIntersectTriPair(
  CVec3d &P0, CVec3d &P1,
  int itri, int jtri,
  const std::vector<unsigned int> &aTri,
  const std::vector<double> &aXYZ);
#endif

// ----------------------------------------------------------------------------

template<typename T>
delfem2::CVec3<T> delfem2::ProjectPointOnTriangle
  (const CVec3<T> &p0,
   const CVec3<T> &tri_p1, const CVec3<T> &tri_p2, const CVec3<T> &tri_p3) {
  CVec3<T> normal = Cross(tri_p2 - tri_p1, tri_p3 - tri_p1);
  double cosAlpha = Dot(p0 - tri_p1, normal) / (Length(p0 - tri_p1) * Length(normal));
  double lenP0ProjectedP0 = Length(tri_p1 - p0) * cosAlpha;
  CVec3<T> p0ProjectedP0 = -1 * lenP0ProjectedP0 * normal / Length(normal);

  return p0 + p0ProjectedP0;
}

// ----------------------

template<typename T>
bool delfem2::isRayIntersectingTriangle
  (const CVec3<T> &line0, const CVec3<T> &line1,
   const CVec3<T> &tri0, const CVec3<T> &tri1, const CVec3<T> &tri2,
   CVec3<T> &intersectionPoint) {
  CVec3<T> normal = Cross(tri1 - tri0, tri2 - tri0);

  // The ray is parallel to the triangle plane
  if (Dot(normal, line1 - line0) == 0) {
    return false;
  }

  double r = Dot(normal, tri0 - line0) / Dot(normal, line1 - line0);

  // The ray does not intersect the triangle plane
  if (r < 0) {
    return false;
  }

  // Find the intersection point
  intersectionPoint = line0 + r * (line1 - line0);

  if (!isPointInsideTriangle(intersectionPoint,
                             tri0, tri1, tri2)) {
    return false;
  }

  return true;
}
#ifdef DFM2_STATIC_LIBRARY
template bool delfem2::isRayIntersectingTriangle(
  const CVec3d &line0, const CVec3d &line1,
  const CVec3d &tri0, const CVec3d &tri1, const CVec3d &tri2,
  CVec3d &intersectionPoint);
#endif

// ----------------------------------------

template<typename T>
bool delfem2::isPointInsideTriangle
  (const CVec3<T> &p0,
   const CVec3<T> &tri_p1, const CVec3<T> &tri_p2, const CVec3<T> &tri_p3) {
  if (isPointSameSide(p0, tri_p1, tri_p2, tri_p3)
    && isPointSameSide(p0, tri_p2, tri_p1, tri_p3)
    && isPointSameSide(p0, tri_p3, tri_p1, tri_p2)) {
    return true;
  } else {
    return false;
  }
}

template<typename T>
bool delfem2::isPointSameSide
  (const CVec3<T> &p0, const CVec3<T> &p1,
   const CVec3<T> &line_p0, const CVec3<T> &line_p1) {
  CVec3<T> crossProd1 = Cross(line_p1 - line_p0, p0 - line_p0);
  CVec3<T> crossProd2 = Cross(line_p1 - line_p0, p1 - line_p0);

  if (Dot(crossProd1, crossProd2) >= 0) {
    return true;
  } else {
    return false;
  }
}
