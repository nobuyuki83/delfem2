/*
 * Copyright (c) 2019 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */

#include "delfem2/msh_boundingbox.h"

#include <cassert>
#include <cmath>
#include <vector>
#include <functional>
#include <array>

#ifndef M_PI
#  define M_PI 3.14159265358979323846
#endif

// ------------------------------------------------

namespace delfem2::msh_boundingbox {

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

// -----------------------------------------------------------------------------

void delfem2::GetCenterWidthGroup(
    double &cx, double &cy, double &cz,
    double &wx, double &wy, double &wz,
    // ----------
    const std::vector<double> &vtx_xyz,
    const std::vector<unsigned int> &elm_vtx,
    unsigned int nnoel,
    int igroup,
    const std::vector<int> &elm_groupidx) {
  const std::size_t nelem = elm_vtx.size() / nnoel;
  assert(elm_vtx.size() == nelem * nnoel);
  assert(elm_groupidx.size() == nelem);
  bool is_ini = true;
  double x_min = 0, x_max = 0, y_min = 0, y_max = 0, z_min = 0, z_max = 0;
  for (unsigned int ielem = 0; ielem < nelem; ielem++) {
    if (elm_groupidx[ielem] != igroup) { continue; }
    for (unsigned int inotri = 0; inotri < nnoel; inotri++) {
      const unsigned int ip = elm_vtx[ielem * 3 + inotri];
      if (is_ini) {
        x_min = x_max = vtx_xyz[ip * 3 + 0];
        y_min = y_max = vtx_xyz[ip * 3 + 1];
        z_min = z_max = vtx_xyz[ip * 3 + 2];
        is_ini = false;
        continue;
      }
      msh_boundingbox::updateMinMaxXYZ(
          x_min, x_max, y_min, y_max, z_min, z_max,
          vtx_xyz[ip * 3 + 0], vtx_xyz[ip * 3 + 1], vtx_xyz[ip * 3 + 2]);
    }
  }
  if (is_ini) {
    cx = cy = cz = 0;
    wx = wy = wz = 1;
    return;
  }
  msh_boundingbox::CenterWidth_MinMaxXYZ(
      cx, cy, cz, wx, wy, wz,
      x_min, x_max, y_min, y_max, z_min, z_max);
}

void delfem2::GetCenterWidthGroup(
    double &cx, double &cy, double &cz,
    double &wx, double &wy, double &wz,
    //
    const std::vector<double> &aXYZ,
    const std::vector<unsigned int> &aElemInd,
    const std::vector<unsigned int> &aElem,
    int igroup,
    const std::vector<int> &aIndGroup) {
  assert(!aElemInd.empty());
  const std::size_t nelem = aElemInd.size() - 1;
  assert(aIndGroup.size() == nelem);
  bool is_ini = true;
  double x_min = 0, x_max = 0, y_min = 0, y_max = 0, z_min = 0, z_max = 0;
  for (unsigned int ielem = 0; ielem < nelem; ielem++) {
    if (aIndGroup[ielem] != igroup) { continue; }
    for (unsigned int iip = aElemInd[ielem]; iip < aElemInd[ielem + 1]; iip++) {
      const unsigned int ip = aElem[iip];
      if (is_ini) {
        x_min = x_max = aXYZ[ip * 3 + 0];
        y_min = y_max = aXYZ[ip * 3 + 1];
        z_min = z_max = aXYZ[ip * 3 + 2];
        is_ini = false;
        continue;
      }
      msh_boundingbox::updateMinMaxXYZ(
          x_min, x_max, y_min, y_max, z_min, z_max,
          aXYZ[ip * 3 + 0], aXYZ[ip * 3 + 1], aXYZ[ip * 3 + 2]);
    }
  }
  if (is_ini) {
    cx = cy = cz = 0;
    wx = wy = wz = 1;
    return;
  }
  msh_boundingbox::CenterWidth_MinMaxXYZ(
      cx, cy, cz, wx, wy, wz,
      x_min, x_max, y_min, y_max, z_min, z_max);
}

void delfem2::GetCenterWidth3DGroup(
    double cw[6],
    //
    const std::vector<double> &aXYZ,
    const std::vector<unsigned int> &aElemInd,
    const std::vector<unsigned int> &aElem,
    int igroup,
    const std::vector<int> &aIndGroup) {
  GetCenterWidthGroup(
      cw[0], cw[1], cw[2], cw[3], cw[4], cw[5],
      aXYZ, aElemInd, aElem, igroup, aIndGroup);
}
