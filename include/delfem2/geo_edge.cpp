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
#include "delfem2/vec3_funcs.h"

#ifndef M_PI
#  define M_PI 3.14159265358979323846
#endif

// ------------------

// r0==0 -> p0==org
// r0==1 -> p1==org
template<class VEC>
VEC delfem2::Nearest_Origin_Edge(
  typename VEC::Scalar &r0,
  const VEC &p0, // start
  const VEC &p1) // end
{
  using SCALAR = typename VEC::Scalar;
  VEC d = p1 - p0;
  SCALAR a = d.dot(d);
  if (a < 1.0e-20) {
    r0 = static_cast<SCALAR>(0.5);
    return (p0 + p1) * r0;
  }
  SCALAR b = d.dot(p0);
  r0 = -b / a;
  if (r0 < 0) { r0 = 0; }
  if (r0 > 1) { r0 = 1; }
  return (1 - r0) * p0 + r0 * p1;
}
#ifdef DFM2_STATIC_LIBRARY
template delfem2::CVec3d delfem2::Nearest_Origin_Edge(
  double &r0,
  const CVec3d &p0, // start
  const CVec3d &p1); // end
#endif

// ------------------------------------

template<class VEC>
VEC delfem2::Nearest_Origin_Edge(
  const VEC &s, // start
  const VEC &e) // end
{
  using SCALAR = typename VEC::Scalar;
  SCALAR t;
  return Nearest_Origin_Edge(t, s, e);
}
#ifdef DFM2_STATIC_LIBRARY
template delfem2::CVec3d delfem2::Nearest_Origin_Edge(
  const CVec3d &s,
  const CVec3d &e);
template delfem2::CVec2d delfem2::Nearest_Origin_Edge(
  const CVec2d &s,
  const CVec2d &e);
#endif

// ---------------------------------------

template<class VEC>
VEC delfem2::Nearest_Edge_Point(
  typename VEC::Scalar &t,
  const VEC &p, // point
  const VEC &s, // source
  const VEC &e) // end
{
  return Nearest_Origin_Edge(t, s - p, e - p) + p;
}
#ifdef DFM2_STATIC_LIBRARY
template delfem2::CVec3d delfem2::Nearest_Edge_Point(
  double &t,
  const CVec3d &p, // point
  const CVec3d &s, // source
  const CVec3d &e); // end
template delfem2::CVec3f delfem2::Nearest_Edge_Point(
  float &t,
  const CVec3f &p, // point
  const CVec3f &s, // source
  const CVec3f &e); // end
#endif

// ------------------------------

template<class VEC>
VEC delfem2::Nearest_Edge_Point(
  const VEC &po_c,
  const VEC &po_s,
  const VEC &po_e) {
  using SCALAR = typename VEC::Scalar;
  SCALAR t;
  return Nearest_Edge_Point(t, po_c, po_s, po_e);
}
#ifdef DFM2_STATIC_LIBRARY
template delfem2::CVec2d delfem2::Nearest_Edge_Point(
  const CVec2d &po_c, const CVec2d &po_s, const CVec2d &po_e);
template delfem2::CVec3d delfem2::Nearest_Edge_Point(
  const CVec3d &po_c, const CVec3d &po_s, const CVec3d &po_e);
#endif

// --------------------------------------

template<class VEC>
typename VEC::Scalar delfem2::Distance_Edge_Point(
  const VEC &po_c,
  const VEC &po_s,
  const VEC &po_e) {
  return Nearest_Origin_Edge(po_s - po_c, po_e - po_c).norm();
}
#ifdef DFM2_STATIC_LIBRARY
template double delfem2::Distance_Edge_Point(
  const CVec2d &po_c, const CVec2d &po_s, const CVec2d &po_e);
template double delfem2::Distance_Edge_Point(
  const CVec3d &po_c, const CVec3d &po_s, const CVec3d &po_e);
#endif

// above: 2D and 3D
// =========================================================
// below: 2D

template<typename VEC>
bool delfem2::IsIntersect_Edge2_Edge2(
  const VEC &po_s0,
  const VEC &po_e0,
  const VEC &po_s1,
  const VEC &po_e1) {
  using SCALAR = typename VEC::Scalar;
  {
    const SCALAR min0x = (po_s0.p[0] < po_e0.p[0]) ? po_s0.p[0] : po_e0.p[0];
    const SCALAR max0x = (po_s0.p[0] > po_e0.p[0]) ? po_s0.p[0] : po_e0.p[0];
    const SCALAR max1x = (po_s1.p[0] > po_e1.p[0]) ? po_s1.p[0] : po_e1.p[0];
    const SCALAR min1x = (po_s1.p[0] < po_e1.p[0]) ? po_s1.p[0] : po_e1.p[0];
    const SCALAR min0y = (po_s0.p[1] < po_e0.p[1]) ? po_s0.p[1] : po_e0.p[1];
    const SCALAR max0y = (po_s0.p[1] > po_e0.p[1]) ? po_s0.p[1] : po_e0.p[1];
    const SCALAR max1y = (po_s1.p[1] > po_e1.p[1]) ? po_s1.p[1] : po_e1.p[1];
    const SCALAR min1y = (po_s1.p[1] < po_e1.p[1]) ? po_s1.p[1] : po_e1.p[1];
    const SCALAR len = ((max0x - min0x) + (max0y - min0y) + (max1x - min1x) + (max1y - min1y)) * 0.0001;
    //    std::cout << len << std::endl;
    if (max1x + len < min0x) return false;
    if (max0x + len < min1x) return false;
    if (max1y + len < min0y) return false;
    if (max0y + len < min1y) return false;
  }
  const SCALAR area1 = Area_Tri2(po_s0, po_e0, po_s1);
  const SCALAR area2 = Area_Tri2(po_s0, po_e0, po_e1);
  const SCALAR area3 = Area_Tri2(po_s1, po_e1, po_s0);
  const SCALAR area4 = Area_Tri2(po_s1, po_e1, po_e0);
  //  std::cout << area1 << " " << area2 << " " << area3 << " " << area4 << std::endl;
  const double a12 = area1 * area2;
  if (a12 > 0) return false;
  const double a34 = area3 * area4;
  if (a34 > 0) return false;
  return true;
}
#ifdef DFM2_STATIC_LIBRARY
template bool delfem2::IsIntersect_Edge2_Edge2(
  const delfem2::CVec2d &po_s0,
  const delfem2::CVec2d &po_e0,
  const delfem2::CVec2d &po_s1,
  const delfem2::CVec2d &po_e1);
template bool delfem2::IsIntersect_Edge2_Edge2(
  const delfem2::CVec2f &po_s0,
  const delfem2::CVec2f &po_e0,
  const delfem2::CVec2f &po_s1,
  const delfem2::CVec2f &po_e1);
#endif

template<typename T>
double delfem2::Distance_Edge2_Edge2(
  const CVec2<T> &po_s0,
  const CVec2<T> &po_e0,
  const CVec2<T> &po_s1,
  const CVec2<T> &po_e1) {
  if (IsCross_LineSeg_LineSeg(po_s0, po_e0, po_s1, po_e1)) return -1;
  const double ds1 = Distance_Edge_Point(po_s0, po_s1, po_e1);
  const double de1 = Distance_Edge_Point(po_e0, po_s1, po_e1);
  const double ds0 = Distance_Edge_Point(po_s1, po_s0, po_e0);
  const double de0 = Distance_Edge_Point(po_e1, po_s0, po_e0);
  double min_dist = ds1;
  min_dist = (de1 < min_dist) ? de1 : min_dist;
  min_dist = (ds0 < min_dist) ? ds0 : min_dist;
  min_dist = (de0 < min_dist) ? de0 : min_dist;
  return min_dist;
}

template<typename VEC>
bool IsIntersect_AABB2_Edge2(
  const VEC &pmin,
  const VEC &pmax,
  const VEC &p0,
  const VEC &p1) {
  typedef typename VEC::Scalar SCALAR;

  assert(pmin[0]<=pmax[0]);
  assert(pmin[1]<=pmax[1]);

  {
    SCALAR xmin, xmax;
    // x axis
    if (p0[0] > p1[0]) {
      xmax = p0[0];
      xmin = p1[0];
    } else {
      xmax = p1[0];
      xmin = p0[0];
    }
    if (xmin > pmax[0]) { return false; }
    if (xmax < pmin[0]) { return false; }
  }

  { // y axis
    SCALAR ymin, ymax;
    if (p0[1] > p1[1]) {
      ymax = p0[1];
      ymin = p1[1];
    } else {
      ymax = p1[1];
      ymin = p0[1];
    }
    if (ymin > pmax[1]) { return false; }
    if (ymax < pmin[1]) { return false; }
  }

  // normal axis
  const VEC direction = p1 - p0;
  const VEC normal = {direction[1], -direction[0]};

  bool positive = false, negative = false;
  if (normal.dot(pmin - p0) > 0) {
    positive = true;
  } else {
    negative = true;
  }

  if (normal.dot(pmax - p0) > 0) {
    positive = true;
  } else {
    negative = true;
  }

  if (normal.dot(VEC(pmin[0], pmax[1]) - p0) > 0) {
    positive = true;
  } else {
    negative = true;
  }

  if (normal.dot(VEC(pmax[0], pmin[1]) - p0) > 0) {
    positive = true;
  } else {
    negative = true;
  }

  if (positive && negative) { return true; }

  return false;
}
#ifdef DFM2_STATIC_LIBRARY
template bool IsIntersect_AABB2_Edge2(
  const delfem2::CVec2d &pmin,
  const delfem2::CVec2d &pmax,
  const delfem2::CVec2d &p0,
  const delfem2::CVec2d &p1);
template bool IsIntersect_AABB2_Edge2(
  const delfem2::CVec2f &pmin,
  const delfem2::CVec2f &pmax,
  const delfem2::CVec2f &p0,
  const delfem2::CVec2f &p1);
#endif


// above: 2D
// =========================
// below: 3D

DFM2_INLINE void delfem2::Nearest_Edge3_Point3(
  double nearest_position[3],
  const double point_pos[3], // point
  const double edge_pos0[3], // source
  const double edge_pos1[3]) // end
{
  const double d[3] = {
    edge_pos1[0] - edge_pos0[0],
    edge_pos1[1] - edge_pos0[1],
    edge_pos1[2] - edge_pos0[2]};
  double t = 0.5;
  if (Dot3(d, d) > 1.0e-20) {
    const double ps[3] = {
      edge_pos0[0] - point_pos[0],
      edge_pos0[1] - point_pos[1],
      edge_pos0[2] - point_pos[2]};
    double a = Dot3(d, d);
    double b = Dot3(d, ps);
    t = -b / a;
    if (t < 0) t = 0;
    if (t > 1) t = 1;
  }
  nearest_position[0] = edge_pos0[0] + t * d[0];
  nearest_position[1] = edge_pos0[1] + t * d[1];
  nearest_position[2] = edge_pos0[2] + t * d[2];
}

// ---------------------------------------------------------------------------

template<class VEC>
void delfem2::Nearest_Edge3_Line3(
  VEC &nearest_lineseg,
  VEC &nearest_line,
  const VEC &lineseg_start,
  const VEC &lineseg_end,
  const VEC &line_origin,
  const VEC &line_direction) {
  using T = decltype(lineseg_start[0]);
  T D0, Dta0, Dtb0;
  VEC Da0, Db0;
  Nearest_Line3_Line3(
    D0, Da0, Db0, Dta0, Dtb0,
    lineseg_start, lineseg_end - lineseg_start, line_origin, line_direction);
  if (abs(D0) < 1.0e-10) { // pararell
    nearest_lineseg = (lineseg_start + lineseg_end) * static_cast<T>(0.5);
    nearest_line = ::delfem2::Nearest_Line3_Point3<VEC>(
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
  const VEC p1 = Nearest_Line3_Point3<VEC>(lineseg_start, line_origin, line_direction);
  const VEC p2 = Nearest_Line3_Point3<VEC>(lineseg_end, line_origin, line_direction);
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
template void delfem2::Nearest_Edge3_Line3<delfem2::CVec3f>(
  CVec3f &, CVec3f &,
  const CVec3f &, const CVec3f &,
  const CVec3f &, const CVec3f &);
template void delfem2::Nearest_Edge3_Line3<delfem2::CVec3d>(
  CVec3d &, CVec3d &,
  const CVec3d &, const CVec3d &,
  const CVec3d &, const CVec3d &);
#endif

// -----------------------------

//　distance EE
template<typename T>
double delfem2::Distance_Edge3_Edge3(
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
    CVec3<T> vert = pq0 - pq0.dot(nvp) * nvp;
    double dist = vert.norm();
    double lp0 = p0.dot(nvp);
    double lp1 = p1.dot(nvp);
    double lq0 = q0.dot(nvp);
    double lq1 = q1.dot(nvp);
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
  double t0 = vp.dot(vp);
  double t1 = vq.dot(vq);
  double t2 = vp.dot(vq);
  double t3 = vp.dot(q0 - p0);
  double t4 = vq.dot(q0 - p0);
  double det = t0 * t1 - t2 * t2;
  double invdet = 1.0 / det;
  ratio_p = (+t1 * t3 - t2 * t4) * invdet;
  ratio_q = (+t2 * t3 - t0 * t4) * invdet;
  CVec3<T> pc = p0 + ratio_p * vp;
  CVec3<T> qc = q0 + ratio_q * vq;
  return (pc - qc).norm();
}
#ifdef DFM2_STATIC_LIBRARY
template double delfem2::Distance_Edge3_Edge3(
  const CVec3d &p0, const CVec3d &p1,
  const CVec3d &q0, const CVec3d &q1,
  double &ratio_p, double &ratio_q);
#endif

// ------------------------

template<typename T>
bool delfem2::IsContact_Edge3_Edge3_Proximity(
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
  double dist = Distance_Edge3_Edge3(p0, p1, q0, q1, ratio_p, ratio_q);
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
template bool delfem2::IsContact_Edge3_Edge3_Proximity(
  int ino0, int ino1,
  int jno0, int jno1,
  const CVec3d &p0, const CVec3d &p1,
  const CVec3d &q0, const CVec3d &q1,
  const double delta);
#endif
