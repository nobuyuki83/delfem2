/*
 * Copyright (c) 2019 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */

#include "delfem2/geo_edge.h"

#include <cstdlib>
#include <cmath>

#include "delfem2/geo_line.h"

#ifndef M_PI
#  define M_PI 3.14159265358979323846
#endif

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

// ---------------------------------------------------------------------------

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
#ifdef DFM2_STATIC_LIBRARY
template delfem2::CVec3d delfem2::Nearest_Origin3_LineSeg3(
  const CVec3d &s,
  const CVec3d &e);
#endif

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
#ifdef DFM2_STATIC_LIBRARY
template delfem2::CVec3d delfem2::Nearest_Origin3_LineSeg3(
  double &r0,
  const CVec3d &p0, // start
  const CVec3d &p1); // end
#endif

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


// ===================


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

// =============

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