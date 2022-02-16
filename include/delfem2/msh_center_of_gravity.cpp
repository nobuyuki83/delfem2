/*
 * Copyright (c) 2019 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */

#include "delfem2/msh_center_of_gravity.h"

#include <cassert>
#include <cmath>
#include <vector>
#include <functional>
#include <array>

#ifndef M_PI
#  define M_PI 3.14159265358979323846
#endif

// ------------------------------------------------

namespace delfem2::cnter_of_gvavity {

//! @details we have "float" and "double" versions Length3 because of sqrtf and sqrt
DFM2_INLINE double Length3(const double p[3]) {
  return sqrt(p[0] * p[0] + p[1] * p[1] + p[2] * p[2]);
}

//! @details we have "float" and "double" versions Length3 because of sqrtf and sqrt
DFM2_INLINE float Length3(const float p[3]) {
  return sqrtf(p[0] * p[0] + p[1] * p[1] + p[2] * p[2]);
}

DFM2_INLINE void Cross3D(double r[3], const double v1[3], const double v2[3]) {
  r[0] = v1[1] * v2[2] - v2[1] * v1[2];
  r[1] = v1[2] * v2[0] - v2[2] * v1[0];
  r[2] = v1[0] * v2[1] - v2[0] * v1[1];
}

DFM2_INLINE double TriArea2D(const double p0[], const double p1[], const double p2[]) {
  return 0.5 * ((p1[0] - p0[0]) * (p2[1] - p0[1]) - (p2[0] - p0[0]) * (p1[1] - p0[1]));
}

template<typename T>
DFM2_INLINE T TriArea3D(
    const T v1[3],
    const T v2[3],
    const T v3[3]) {
  T n[3];
  n[0] = (v2[1] - v1[1]) * (v3[2] - v1[2]) - (v3[1] - v1[1]) * (v2[2] - v1[2]);
  n[1] = (v2[2] - v1[2]) * (v3[0] - v1[0]) - (v3[2] - v1[2]) * (v2[0] - v1[0]);
  n[2] = (v2[0] - v1[0]) * (v3[1] - v1[1]) - (v3[0] - v1[0]) * (v2[1] - v1[1]);
  return std::sqrt(n[0] * n[0] + n[1] * n[1] + n[2] * n[2]) / 2;
}

template<typename T>
DFM2_INLINE void UnitNormalAreaTri3(
    T n[3],
    T &a,
    const T v1[3],
    const T v2[3],
    const T v3[3]) {
  n[0] = (v2[1] - v1[1]) * (v3[2] - v1[2]) - (v3[1] - v1[1]) * (v2[2] - v1[2]);
  n[1] = (v2[2] - v1[2]) * (v3[0] - v1[0]) - (v3[2] - v1[2]) * (v2[0] - v1[0]);
  n[2] = (v2[0] - v1[0]) * (v3[1] - v1[1]) - (v3[0] - v1[0]) * (v2[1] - v1[1]);
  a = std::sqrt(n[0] * n[0] + n[1] * n[1] + n[2] * n[2]) / 2;
  const T invlen = 1 / (a * 2);
  n[0] *= invlen;
  n[1] *= invlen;
  n[2] *= invlen;
}

template<typename T>
std::array<T,3> Normal_Tri3(
    const T v1[3],
    const T v2[3],
    const T v3[3]) {
  return {
      (v2[1] - v1[1]) * (v3[2] - v1[2]) - (v3[1] - v1[1]) * (v2[2] - v1[2]),
      (v2[2] - v1[2]) * (v3[0] - v1[0]) - (v3[2] - v1[2]) * (v2[0] - v1[0]),
      (v2[0] - v1[0]) * (v3[1] - v1[1]) - (v3[0] - v1[0]) * (v2[1] - v1[1]) };
}

template<typename T>
DFM2_INLINE void MatVec3(T y[3],
                         const T m[9], const T x[3]) {
  y[0] = m[0] * x[0] + m[1] * x[1] + m[2] * x[2];
  y[1] = m[3] * x[0] + m[4] * x[1] + m[5] * x[2];
  y[2] = m[6] * x[0] + m[7] * x[1] + m[8] * x[2];
}

//! @details we have "float" and "double" versions Distance3 because of sqrtf and sqrt
DFM2_INLINE double Distance3(const double p0[3], const double p1[3]) {
  return sqrt(
      (p1[0] - p0[0]) * (p1[0] - p0[0]) +
      (p1[1] - p0[1]) * (p1[1] - p0[1]) +
      (p1[2] - p0[2]) * (p1[2] - p0[2]));
}

//! @details we have "float" and "double" versions Distance3 because of sqrtf and sqrt
DFM2_INLINE float Distance3(const float p0[3], const float p1[3]) {
  return sqrtf(
      (p1[0] - p0[0]) * (p1[0] - p0[0]) +
      (p1[1] - p0[1]) * (p1[1] - p0[1]) +
      (p1[2] - p0[2]) * (p1[2] - p0[2]));
}

DFM2_INLINE double Distance2D(const double p0[3], const double p1[3]) {
  return sqrt((p1[0] - p0[0]) * (p1[0] - p0[0]) + (p1[1] - p0[1]) * (p1[1] - p0[1]));
}

DFM2_INLINE double Dot3(const double p0[3], const double p1[3]) {
  return p0[0] * p1[0] + p0[1] * p1[1] + p0[2] * p1[2];
}

DFM2_INLINE double largest(double x0, double x1, double x2) {
  double wmax = x0;
  wmax = (x1 > wmax) ? x1 : wmax;
  wmax = (x2 > wmax) ? x2 : wmax;
  return wmax;
}

template<typename T>
DFM2_INLINE  T TetVolume3D(
    const T v1[3],
    const T v2[3],
    const T v3[3],
    const T v4[3]) {
  return
      ((v2[0] - v1[0]) * ((v3[1] - v1[1]) * (v4[2] - v1[2]) - (v4[1] - v1[1]) * (v3[2] - v1[2]))
          - (v2[1] - v1[1]) * ((v3[0] - v1[0]) * (v4[2] - v1[2]) - (v4[0] - v1[0]) * (v3[2] - v1[2]))
          + (v2[2] - v1[2]) * ((v3[0] - v1[0]) * (v4[1] - v1[1]) - (v4[0] - v1[0]) * (v3[1] - v1[1]))
      ) * static_cast<T>(1.0 / 6.0);
}

DFM2_INLINE void Mat3_Bryant
    (double m[9],
     double rx, double ry, double rz) {
  m[0] = cos(rz) * cos(ry);
  m[1] = cos(rz) * sin(ry) * sin(rx) - sin(rz) * cos(rx);
  m[2] = cos(rz) * sin(ry) * cos(rx) + sin(rz) * sin(rx);
  m[3] = sin(rz) * cos(ry);
  m[4] = sin(rz) * sin(ry) * sin(rx) + cos(rz) * cos(rx);
  m[5] = sin(rz) * sin(ry) * cos(rx) - cos(rz) * sin(rx);
  m[6] = -sin(ry);
  m[7] = cos(ry) * sin(rx);
  m[8] = cos(ry) * cos(rx);
}

DFM2_INLINE void Mat3_Bryant
    (float m[9],
     float rx, float ry, float rz) {
  m[0] = cosf(rz) * cosf(ry);
  m[1] = cosf(rz) * sinf(ry) * sinf(rx) - sinf(rz) * cosf(rx);
  m[2] = cosf(rz) * sinf(ry) * cosf(rx) + sinf(rz) * sinf(rx);
  m[3] = sinf(rz) * cosf(ry);
  m[4] = sinf(rz) * sinf(ry) * sinf(rx) + cosf(rz) * cosf(rx);
  m[5] = sinf(rz) * sinf(ry) * cosf(rx) - cosf(rz) * sinf(rx);
  m[6] = -sinf(ry);
  m[7] = cosf(ry) * sinf(rx);
  m[8] = cosf(ry) * cosf(rx);
}

template<typename T>
void CenterWidth_MinMaxXYZ
    (T &cx, T &cy, T &cz,
     T &wx, T &wy, T &wz,
        //
     T x_min, T x_max,
     T y_min, T y_max,
     T z_min, T z_max) {
  cx = (x_min + x_max) * 0.5;
  cy = (y_min + y_max) * 0.5;
  cz = (z_min + z_max) * 0.5;
  wx = x_max - x_min;
  wy = y_max - y_min;
  wz = z_max - z_min;
}

template<typename T>
void updateMinMaxXYZ(
    T &x_min, T &x_max,
    T &y_min, T &y_max,
    T &z_min, T &z_max,
    T x, T y, T z) {
  if (x_min > x_max) {
    x_min = x_max = x;
    y_min = y_max = y;
    z_min = z_max = z;
    return;
  }
  x_min = (x_min < x) ? x_min : x;
  x_max = (x_max > x) ? x_max : x;
  y_min = (y_min < y) ? y_min : y;
  y_max = (y_max > y) ? y_max : y;
  z_min = (z_min < z) ? z_min : z;
  z_max = (z_max > z) ? z_max : z;
}
#ifdef DFM2_STATIC_LIBRARY
template void updateMinMaxXYZ(
    float &x_min, float &x_max,
    float &y_min, float &y_max,
    float &z_min, float &z_max,
    float X, float Y, float Z);
template void updateMinMaxXYZ(
    double &x_min, double &x_max,
    double &y_min, double &y_max,
    double &z_min, double &z_max,
    double X, double Y, double Z);
#endif

}

// static function above
// ==============================================
// exposed function below


template<typename T>
void delfem2::CG_Point3(
    T *center_of_gravity,
    const std::vector<T> &vtx_xyz) {
  center_of_gravity[0] = center_of_gravity[1] = center_of_gravity[2] = 0;
  const size_t nXYZ = vtx_xyz.size() / 3;
  for (unsigned int ixyz = 0; ixyz < nXYZ; ixyz++) {
    center_of_gravity[0] += vtx_xyz[ixyz * 3 + 0];
    center_of_gravity[1] += vtx_xyz[ixyz * 3 + 1];
    center_of_gravity[2] += vtx_xyz[ixyz * 3 + 2];
  }
  center_of_gravity[0] /= nXYZ;
  center_of_gravity[1] /= nXYZ;
  center_of_gravity[2] /= nXYZ;
}
#ifdef DFM2_STATIC_LIBRARY
template void delfem2::CG_Point3(float *, const std::vector<float> &);
template void delfem2::CG_Point3(double *, const std::vector<double> &);
#endif


// ---------------------

template<typename T>
T delfem2::CentsMaxRad_MeshTri3(
    std::vector<T> &tri_centerxyz,
    const std::vector<T> &vtx_xyz,
    const std::vector<unsigned int> &tri_vtx) {
  namespace lcl = delfem2::cnter_of_gvavity;
  T max_rad0 = -1;
  const size_t nTri = tri_vtx.size() / 3;
  tri_centerxyz.resize(nTri * 3);
  for (unsigned int itri = 0; itri < nTri; ++itri) {
    const unsigned int i0 = tri_vtx[itri * 3 + 0];
    const unsigned int i1 = tri_vtx[itri * 3 + 1];
    const unsigned int i2 = tri_vtx[itri * 3 + 2];
    const T pc[3] = {
        (vtx_xyz[i0 * 3 + 0] + vtx_xyz[i1 * 3 + 0] + vtx_xyz[i2 * 3 + 0]) / 3,
        (vtx_xyz[i0 * 3 + 1] + vtx_xyz[i1 * 3 + 1] + vtx_xyz[i2 * 3 + 1]) / 3,
        (vtx_xyz[i0 * 3 + 2] + vtx_xyz[i1 * 3 + 2] + vtx_xyz[i2 * 3 + 2]) / 3};
    tri_centerxyz[itri * 3 + 0] = pc[0];
    tri_centerxyz[itri * 3 + 1] = pc[1];
    tri_centerxyz[itri * 3 + 2] = pc[2];
    const T l0 = lcl::Distance3(pc, vtx_xyz.data() + i0 * 3);
    const T l1 = lcl::Distance3(pc, vtx_xyz.data() + i1 * 3);
    const T l2 = lcl::Distance3(pc, vtx_xyz.data() + i2 * 3);
    if (max_rad0 < 0 || l0 > max_rad0) { max_rad0 = l0; }
    if (max_rad0 < 0 || l1 > max_rad0) { max_rad0 = l1; }
    if (max_rad0 < 0 || l2 > max_rad0) { max_rad0 = l2; }
  }
  return max_rad0;
}
#ifdef DFM2_STATIC_LIBRARY
template float delfem2::CentsMaxRad_MeshTri3(
    std::vector<float> &tri_centerxyz,
    const std::vector<float> &vtx_xyz,
    const std::vector<unsigned int> &tri_vtx);
template double delfem2::CentsMaxRad_MeshTri3(
    std::vector<double> &tri_centerxyz,
    const std::vector<double> &vtx_xyz,
    const std::vector<unsigned int> &tri_vtx);
#endif

// ----------------------------------------------

template<typename T>
void delfem2::CG_MeshTri3_Solid(
    T center_of_gravity_xyz[3],
    const std::vector<T> &vtx_xyz,
    const std::vector<unsigned int> &tri_vtx) { // center of gravity
  namespace lcl = delfem2::cnter_of_gvavity;
  center_of_gravity_xyz[0] = 0;
  center_of_gravity_xyz[1] = 0;
  center_of_gravity_xyz[2] = 0;
  T tw = 0;
  const size_t nTri = tri_vtx.size() / 3;
  constexpr T quarter = static_cast<T>(0.25);
  for (unsigned int itri = 0; itri < nTri; itri++) {
    const unsigned int i1 = tri_vtx[itri * 3 + 0];
    const unsigned int i2 = tri_vtx[itri * 3 + 1];
    const unsigned int i3 = tri_vtx[itri * 3 + 2];
    const T q0[3] = {0, 0, 0};
    const T q1[3] = {vtx_xyz[i1 * 3 + 0], vtx_xyz[i1 * 3 + 1], vtx_xyz[i1 * 3 + 2]};
    const T q2[3] = {vtx_xyz[i2 * 3 + 0], vtx_xyz[i2 * 3 + 1], vtx_xyz[i2 * 3 + 2]};
    const T q3[3] = {vtx_xyz[i3 * 3 + 0], vtx_xyz[i3 * 3 + 1], vtx_xyz[i3 * 3 + 2]};
    T v = lcl::TetVolume3D(q0, q1, q2, q3);
    tw += v;
    center_of_gravity_xyz[0] += (q0[0] + q1[0] + q2[0] + q3[0]) * quarter * v;
    center_of_gravity_xyz[1] += (q0[1] + q1[1] + q2[1] + q3[1]) * quarter * v;
    center_of_gravity_xyz[2] += (q0[2] + q1[2] + q2[2] + q3[2]) * quarter * v;
  }
  center_of_gravity_xyz[0] /= tw;
  center_of_gravity_xyz[1] /= tw;
  center_of_gravity_xyz[2] /= tw;
}
#ifdef DFM2_STATIC_LIBRARY
template void delfem2::CG_MeshTri3_Solid(
    float cg[3],
    const std::vector<float> &aXYZ,
    const std::vector<unsigned int> &aTri);
template void delfem2::CG_MeshTri3_Solid(
    double cg[3],
    const std::vector<double> &aXYZ,
    const std::vector<unsigned int> &aTri);
#endif

// ----------------------------------------

template<typename T>
void delfem2::CG_MeshTri3_Shell(
    T cg[3],
    const std::vector<T> &aXYZ,
    const std::vector<unsigned int> &aTri) {  // center of gravity
  namespace lcl = delfem2::cnter_of_gvavity;
  constexpr T onethird = static_cast<T>(1.0 / 3.0);
  cg[0] = cg[1] = cg[2] = 0;
  T tw = 0;
  const size_t nTri = aTri.size() / 3;
  for (unsigned int itri = 0; itri < nTri; itri++) {
    const unsigned int i1 = aTri[itri * 3 + 0];
    const unsigned int i2 = aTri[itri * 3 + 1];
    const unsigned int i3 = aTri[itri * 3 + 2];
    const T q1[3] = {aXYZ[i1 * 3 + 0], aXYZ[i1 * 3 + 1], aXYZ[i1 * 3 + 2]};
    const T q2[3] = {aXYZ[i2 * 3 + 0], aXYZ[i2 * 3 + 1], aXYZ[i2 * 3 + 2]};
    const T q3[3] = {aXYZ[i3 * 3 + 0], aXYZ[i3 * 3 + 1], aXYZ[i3 * 3 + 2]};
    T a = lcl::TriArea3D(q1, q2, q3);
    tw += a;
    cg[0] += (q1[0] + q2[0] + q3[0]) * onethird * a;
    cg[1] += (q1[1] + q2[1] + q3[1]) * onethird * a;
    cg[2] += (q1[2] + q2[2] + q3[2]) * onethird * a;
  }
  cg[0] /= tw;
  cg[1] /= tw;
  cg[2] /= tw;
}
#ifdef DFM2_STATIC_LIBRARY
template void delfem2::CG_MeshTri3_Shell(
    float cg[3],
    const std::vector<float> &aXYZ,
    const std::vector<unsigned int> &aTri);
template void delfem2::CG_MeshTri3_Shell(
    double cg[3],
    const std::vector<double> &aXYZ,
    const std::vector<unsigned int> &aTri);
#endif

// ------------------------------------------

template<typename T>
T delfem2::CG_TriMsh3Flg_Shell(
    T cg[3],
    const std::vector<T> &aXYZ,
    const std::vector<unsigned int> &aTri,
    int iflg,
    const std::vector<int> &aFlg) {
  namespace lcl = delfem2::cnter_of_gvavity;
  constexpr T onethird = static_cast<T>(1.0 / 3.0);
  cg[0] = cg[1] = cg[2] = 0;
  T tw = 0;
  const std::size_t nTri = aTri.size() / 3;
  for (std::size_t itri = 0; itri < nTri; itri++) {
    if (aFlg[itri] != iflg) continue;
    const unsigned int i1 = aTri[itri * 3 + 0];
    const unsigned int i2 = aTri[itri * 3 + 1];
    const unsigned int i3 = aTri[itri * 3 + 2];
    const T q1[3] = {aXYZ[i1 * 3 + 0], aXYZ[i1 * 3 + 1], aXYZ[i1 * 3 + 2]};
    const T q2[3] = {aXYZ[i2 * 3 + 0], aXYZ[i2 * 3 + 1], aXYZ[i2 * 3 + 2]};
    const T q3[3] = {aXYZ[i3 * 3 + 0], aXYZ[i3 * 3 + 1], aXYZ[i3 * 3 + 2]};
    T a = lcl::TriArea3D(q1, q2, q3);
    tw += a;
    cg[0] += (q1[0] + q2[0] + q3[0]) * onethird * a;
    cg[1] += (q1[1] + q2[1] + q3[1]) * onethird * a;
    cg[2] += (q1[2] + q2[2] + q3[2]) * onethird * a;
  }
  cg[0] /= tw;
  cg[1] /= tw;
  cg[2] /= tw;
  return tw;
}
#ifdef DFM2_STATIC_LIBRARY
template float delfem2::CG_TriMsh3Flg_Shell(
    float cg[3],
    const std::vector<float> &aXYZ,
    const std::vector<unsigned int> &aTri,
    int iflg,
    const std::vector<int> &aFlg);
template double delfem2::CG_TriMsh3Flg_Shell(
    double cg[3],
    const std::vector<double> &aXYZ,
    const std::vector<unsigned int> &aTri,
    int iflg,
    const std::vector<int> &aFlg);
#endif

// ------------------------------------------

void delfem2::CG_Tri(
    double &cgx,
    double &cgy,
    double &cgz,
    int itri,
    const std::vector<double> &aXYZ,
    const std::vector<int> &aTri) { // center of gravity
  assert(itri >= 0 && itri < (int) aTri.size() / 3);
  const int i1 = aTri[itri * 3 + 0];
  const int i2 = aTri[itri * 3 + 1];
  const int i3 = aTri[itri * 3 + 2];
  const double q1[3] = {aXYZ[i1 * 3 + 0], aXYZ[i1 * 3 + 1], aXYZ[i1 * 3 + 2]};
  const double q2[3] = {aXYZ[i2 * 3 + 0], aXYZ[i2 * 3 + 1], aXYZ[i2 * 3 + 2]};
  const double q3[3] = {aXYZ[i3 * 3 + 0], aXYZ[i3 * 3 + 1], aXYZ[i3 * 3 + 2]};
//  double a = TriArea3D(q1, q2, q3);
  cgx = (q1[0] + q2[0] + q3[0]) * 0.333333;
  cgy = (q1[1] + q2[1] + q3[1]) * 0.333333;
  cgz = (q1[2] + q2[2] + q3[2]) * 0.333333;
}

// -------------------------------

template<typename T>
void delfem2::CG_MeshTet3(
    T &v_tot,
    T cg[3],
    const std::vector<T> &aXYZC,
    const std::vector<unsigned int> &aTet) {
  namespace lcl = delfem2::cnter_of_gvavity;
  constexpr T quarter = static_cast<T>(1.0 / 4.0);
  v_tot = cg[0] = cg[1] = cg[2] = 0.0;
  const T *pXYZ = aXYZC.data();
  const std::size_t nTet = aTet.size() / 4;
  for (std::size_t it = 0; it < nTet; ++it) {
    const T *p0 = pXYZ + aTet[it * 4 + 0] * 3;
    const T *p1 = pXYZ + aTet[it * 4 + 1] * 3;
    const T *p2 = pXYZ + aTet[it * 4 + 2] * 3;
    const T *p3 = pXYZ + aTet[it * 4 + 3] * 3;
    const T v = lcl::TetVolume3D(p0, p1, p2, p3);
    v_tot += v;
    cg[0] += v * (p0[0] + p1[0] + p2[0] + p3[0]) * quarter;
    cg[1] += v * (p0[1] + p1[1] + p2[1] + p3[1]) * quarter;
    cg[2] += v * (p0[2] + p1[2] + p2[2] + p3[2]) * quarter;
  }
  cg[0] /= v_tot;
  cg[1] /= v_tot;
  cg[2] /= v_tot;
}
#ifdef DFM2_STATIC_LIBRARY
template void delfem2::CG_MeshTet3(
    float &v_tot, float cg[3],
    const std::vector<float> &aXYZ,
    const std::vector<unsigned int> &aTet);
template void delfem2::CG_MeshTet3(
    double &v_tot, double cg[3],
    const std::vector<double> &aXYZ,
    const std::vector<unsigned int> &aTet);
#endif

// -------------------------------

template<typename T>
void delfem2::CG_MeshLine3(
  T &v_tot,
  T cg[3],
  const std::vector<T> &vtx_xyz,
  const std::vector<unsigned int> &line_vtx) {
  namespace lcl = delfem2::cnter_of_gvavity;
  constexpr T half = static_cast<T>(1.0 / 2.0);
  v_tot = cg[0] = cg[1] = cg[2] = 0.0;
  const std::size_t nline = line_vtx.size() / 2;
  for (std::size_t il = 0; il < nline; ++il) {
    const T *p0 = vtx_xyz.data() + line_vtx[il * 2 + 0] * 3;
    const T *p1 = vtx_xyz.data() + line_vtx[il * 2 + 1] * 3;
    const T v = lcl::Distance3(p0, p1);
    v_tot += v;
    cg[0] += v * (p0[0] + p1[0]) * half;
    cg[1] += v * (p0[1] + p1[1]) * half;
    cg[2] += v * (p0[2] + p1[2]) * half;
  }
  cg[0] /= v_tot;
  cg[1] /= v_tot;
  cg[2] /= v_tot;
}
#ifdef DFM2_STATIC_LIBRARY
template void delfem2::CG_MeshLine3(
    float &v_tot, float cg[3],
    const std::vector<float> &aXYZ,
    const std::vector<unsigned int> &aTet);
template void delfem2::CG_MeshLine3(
    double &v_tot, double cg[3],
    const std::vector<double> &aXYZ,
    const std::vector<unsigned int> &aTet);
#endif
