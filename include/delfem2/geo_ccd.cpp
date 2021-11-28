/*
 * Copyright (c) 2019 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */

#include "delfem2/geo_ccd.h"

#include <cmath>
#include <stack>

#include "delfem2/geo_edge.h"
#include "delfem2/geo_tri.h"

#ifndef M_PI
#  define M_PI 3.14159265358979323846
#endif


// =====================================
// below: unexposed 

namespace delfem2::geo_ccd {

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

}  // delfem2::geo_ccd

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
  namespace lcl = delfem2::geo_ccd;
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
  namespace lcl = delfem2::geo_ccd;
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
      double dist01 = (Nearest_Edge3_Point3(p3, p0, p1) - p3).norm();
      double dist12 = (Nearest_Edge3_Point3(p3, p1, p2) - p3).norm();
      double dist20 = (Nearest_Edge3_Point3(p3, p2, p0) - p3).norm();
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
