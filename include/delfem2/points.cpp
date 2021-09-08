/*
 * Copyright (c) 2019 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */

#include <cassert>
#include <cmath>
#include <vector>
#include <random>

#include "delfem2/points.h"

#ifndef M_PI
#  define M_PI 3.14159265358979323846
#endif

// ------------------------------------------------

namespace delfem2::points {

//! @details we have "float" and "double" versions Length3 because of sqrtf and sqrt
DFM2_INLINE double Length3(const double p[3]) { return sqrt(p[0] * p[0] + p[1] * p[1] + p[2] * p[2]); }
DFM2_INLINE float Length3(const float p[3]) { return sqrtf(p[0] * p[0] + p[1] * p[1] + p[2] * p[2]); }

//! @details we have "float" and "double" versions Length2 because of sqrtf and sqrt
DFM2_INLINE double Length2(const double p[2]) { return sqrt(p[0] * p[0] + p[1] * p[1]); }
DFM2_INLINE float Length2(const float p[2]) { return sqrtf(p[0] * p[0] + p[1] * p[1]); }

DFM2_INLINE void Cross3D(
    double r[3], const double v1[3], const double v2[3]) {
  r[0] = v1[1] * v2[2] - v2[1] * v1[2];
  r[1] = v1[2] * v2[0] - v2[2] * v1[0];
  r[2] = v1[0] * v2[1] - v2[0] * v1[1];
}

DFM2_INLINE double TriArea2D(
    const double p0[], const double p1[], const double p2[]) {
  return 0.5 * (
      (p1[0] - p0[0]) * (p2[1] - p0[1]) -
          (p2[0] - p0[0]) * (p1[1] - p0[1]));
}

DFM2_INLINE double TriArea3D(
    const double v1[3], const double v2[3], const double v3[3]) {
  double n[3];
  n[0] = (v2[1] - v1[1]) * (v3[2] - v1[2]) - (v3[1] - v1[1]) * (v2[2] - v1[2]);
  n[1] = (v2[2] - v1[2]) * (v3[0] - v1[0]) - (v3[2] - v1[2]) * (v2[0] - v1[0]);
  n[2] = (v2[0] - v1[0]) * (v3[1] - v1[1]) - (v3[0] - v1[0]) * (v2[1] - v1[1]);
  return sqrt(n[0] * n[0] + n[1] * n[1] + n[2] * n[2]) * 0.5;
}

template<typename T>
DFM2_INLINE void UnitNormalAreaTri3(
    T n[3], T &a,
    const T v1[3], const T v2[3], const T v3[3]) {
  n[0] = (v2[1] - v1[1]) * (v3[2] - v1[2]) - (v3[1] - v1[1]) * (v2[2] - v1[2]);
  n[1] = (v2[2] - v1[2]) * (v3[0] - v1[0]) - (v3[2] - v1[2]) * (v2[0] - v1[0]);
  n[2] = (v2[0] - v1[0]) * (v3[1] - v1[1]) - (v3[0] - v1[0]) * (v2[1] - v1[1]);
  a = sqrt(n[0] * n[0] + n[1] * n[1] + n[2] * n[2]) * 0.5;
  const T invlen = 0.5 / a;
  n[0] *= invlen;
  n[1] *= invlen;
  n[2] *= invlen;
}

template<typename T>
DFM2_INLINE void MatVec3(
    T y[3],
    const T m[9], const T x[3]) {
  y[0] = m[0] * x[0] + m[1] * x[1] + m[2] * x[2];
  y[1] = m[3] * x[0] + m[4] * x[1] + m[5] * x[2];
  y[2] = m[6] * x[0] + m[7] * x[1] + m[8] * x[2];
}

//! @details we have "float" and "double" versions
//! Distance3 because of sqrtf and sqrt
DFM2_INLINE double Distance3(const double p0[3], const double p1[3]) {
  return sqrt(
      (p1[0] - p0[0]) * (p1[0] - p0[0]) +
          (p1[1] - p0[1]) * (p1[1] - p0[1]) +
          (p1[2] - p0[2]) * (p1[2] - p0[2]));
}

//! @details we have "float" and "double" versions
//! Distance3 because of sqrtf and sqrt
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

template<typename REAL>
DFM2_INLINE REAL largest(REAL x0, REAL x1, REAL x2) {
  REAL wmax = x0;
  wmax = (x1 > wmax) ? x1 : wmax;
  wmax = (x2 > wmax) ? x2 : wmax;
  return wmax;
}

template<typename T>
DFM2_INLINE T TetVolume3D(
    const T v1[3],
    const T v2[3],
    const T v3[3],
    const T v4[3]) {
  return
      (
          (v2[0] - v1[0]) * ((v3[1] - v1[1]) * (v4[2] - v1[2]) - (v4[1] - v1[1]) * (v3[2] - v1[2])) -
          (v2[1] - v1[1]) * ((v3[0] - v1[0]) * (v4[2] - v1[2]) - (v4[0] - v1[0]) * (v3[2] - v1[2])) +
          (v2[2] - v1[2]) * ((v3[0] - v1[0]) * (v4[1] - v1[1]) - (v4[0] - v1[0]) * (v3[1] - v1[1]))
      ) * 0.16666666666666666666666666666667;
}

DFM2_INLINE void Mat3_Bryant(
    double m[9],
    double rx,
    double ry,
    double rz) {
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

DFM2_INLINE void Mat3_Bryant(
    float m[9],
    float rx,
    float ry,
    float rz) {
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

}

// static function above
// ==============================================
// exposed function below

template<typename T>
void CenterWidth_MinMaxXYZ(
    T &cx, T &cy, T &cz,
    T &wx, T &wy, T &wz,
    //
    T x_min, T x_max,
    T y_min, T y_max,
    T z_min, T z_max) {
  cx = (x_min + x_max) / 2;
  cy = (y_min + y_max) / 2;
  cz = (z_min + z_max) / 2;
  wx = x_max - x_min;
  wy = y_max - y_min;
  wz = z_max - z_min;
}

// -----------------------------------------------------------------------------

template<typename T>
void delfem2::updateMinMaxXYZ(
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
template void delfem2::updateMinMaxXYZ(
    float &x_min, float &x_max,
    float &y_min, float &y_max,
    float &z_min, float &z_max,
    float X, float Y, float Z);
template void delfem2::updateMinMaxXYZ(
    double &x_min, double &x_max,
    double &y_min, double &y_max,
    double &z_min, double &z_max,
    double X, double Y, double Z);
#endif

// -----------------------------

template<typename T>
void delfem2::BoundingBox3_Points3(
    T min3[3],
    T max3[3],
    const T *aXYZ,
    const unsigned int nXYZ) {
  min3[0] = +1;
  max3[0] = -1;
  for (unsigned int ixyz = 0; ixyz < nXYZ; ++ixyz) {
    updateMinMaxXYZ(min3[0], max3[0],
                    min3[1], max3[1],
                    min3[2], max3[2],
                    aXYZ[ixyz * 3 + 0], aXYZ[ixyz * 3 + 1], aXYZ[ixyz * 3 + 2]);
  }
}
#ifdef DFM2_STATIC_LIBRARY
template void delfem2::BoundingBox3_Points3(
    double min3[3], double max3[3],
    const double *aXYZ, const unsigned int nXYZ);
template void delfem2::BoundingBox3_Points3(
    float min3[3], float max3[3],
    const float *aXYZ, const unsigned int nXYZ);
#endif

// --------------------------------------------------------------------------------

template<typename T>
void delfem2::CenterWidth_Point3(
    T &cx, T &cy, T &cz,
    T &wx, T &wy, T &wz,
    const T *paXYZ,
    const unsigned int nXYZ) {
  if (paXYZ == nullptr) {
    cx = cy = cz = 0;
    wx = wy = wz = 1;
    return;
  }
  T x_min = paXYZ[0], x_max = paXYZ[0];
  T y_min = paXYZ[1], y_max = paXYZ[1];
  T z_min = paXYZ[2], z_max = paXYZ[2];
  for (unsigned int ino = 0; ino < nXYZ; ino++) {
    updateMinMaxXYZ(x_min, x_max, y_min, y_max, z_min, z_max,
                    paXYZ[ino * 3 + 0], paXYZ[ino * 3 + 1], paXYZ[ino * 3 + 2]);
  }
  CenterWidth_MinMaxXYZ(cx, cy, cz, wx, wy, wz,
                        x_min, x_max, y_min, y_max, z_min, z_max);
}
#ifdef DFM2_STATIC_LIBRARY
template void delfem2::CenterWidth_Point3(
    float &cx, float &cy, float &cz,
    float &wx, float &wy, float &wz,
    const float *paXYZ, const unsigned int nXYZ);
template void delfem2::CenterWidth_Point3(
    double &cx, double &cy, double &cz,
    double &wx, double &wy, double &wz,
    const double *paXYZ, const unsigned int nXYZ);
#endif


// ---------------------------------

template<typename T>
void delfem2::CenterWidth_Points3(
    T &cx, T &cy, T &cz,
    T &wx, T &wy, T &wz,
    const std::vector<T> &aXYZ) {
  const size_t np = aXYZ.size() / 3;
  if (np == 0) {
    cx = cy = cz = 0;
    wx = wy = wz = 1;
    return;
  }
  T x_min = aXYZ[0], x_max = aXYZ[0];
  T y_min = aXYZ[1], y_max = aXYZ[1];
  T z_min = aXYZ[2], z_max = aXYZ[2];
  for (unsigned int ip = 0; ip < np; ++ip) {
    updateMinMaxXYZ(x_min, x_max, y_min, y_max, z_min, z_max,
                    aXYZ[ip * 3 + 0], aXYZ[ip * 3 + 1], aXYZ[ip * 3 + 2]);
  }
  CenterWidth_MinMaxXYZ(cx, cy, cz, wx, wy, wz,
                        x_min, x_max, y_min, y_max, z_min, z_max);
}
#ifdef DFM2_STATIC_LIBRARY
template void delfem2::CenterWidth_Points3(
    float &cx, float &cy, float &cz,
    float &wx, float &wy, float &wz,
    const std::vector<float> &aXYZ);
template void delfem2::CenterWidth_Points3(
    double &cx, double &cy, double &cz,
    double &wx, double &wy, double &wz,
    const std::vector<double> &aXYZ);
#endif

// -------------------------------------

template<typename T>
void delfem2::CenterWidth_Points3(
    T c[3],
    T w[3],
    const std::vector<T> &aXYZ) {
  delfem2::CenterWidth_Points3(c[0], c[1], c[2],
                               w[0], w[1], w[2],
                               aXYZ);
}
#ifdef DFM2_STATIC_LIBRARY
template void delfem2::CenterWidth_Points3(
    float c[3],
    float w[3],
    const std::vector<float> &aXYZ);
template void delfem2::CenterWidth_Points3(
    double c[3],
    double w[3],
    const std::vector<double> &aXYZ);
#endif

void delfem2::GetCenterWidthLocal(
    double &lcx, double &lcy, double &lcz,
    double &lwx, double &lwy, double &lwz,
    const std::vector<double> &aXYZ,
    const double lex[3],
    const double ley[3],
    const double lez[3]) {
  const size_t nno = aXYZ.size() / 3;
  if (nno == 0) {
    lcx = lcy = lcz = 0;
    lwx = lwy = lwz = 1;
    return;
  }
  const double p0[3] = {aXYZ[0], aXYZ[1], aXYZ[2]};
  double x_min = points::Dot3(p0, lex);
  double x_max = x_min;
  double y_min = points::Dot3(p0, ley);
  double y_max = y_min;
  double z_min = points::Dot3(p0, lez);
  double z_max = z_min;
  for (unsigned int ino = 0; ino < nno; ++ino) {
    const double pi[3] = {
        aXYZ[ino * 3 + 0],
        aXYZ[ino * 3 + 1],
        aXYZ[ino * 3 + 2]};
    updateMinMaxXYZ(x_min, x_max, y_min, y_max, z_min, z_max,
                    points::Dot3(pi, lex),
                    points::Dot3(pi, ley),
                    points::Dot3(pi, lez));
  }
  CenterWidth_MinMaxXYZ(lcx, lcy, lcz, lwx, lwy, lwz,
                        x_min, x_max, y_min, y_max, z_min, z_max);
}


// -------------------------------------

template<typename T>
void delfem2::Scale_PointsX(
    std::vector<T> &aXYZ,
    const T s) {
  const std::size_t n = aXYZ.size();
  for (unsigned int i = 0; i < n; ++i) { aXYZ[i] *= s; }
}
#ifdef DFM2_STATIC_LIBRARY
template void delfem2::Scale_PointsX(std::vector<float> &, float);
template void delfem2::Scale_PointsX(std::vector<double> &, double);
#endif

// --------------

template<typename T>
void delfem2::Scale_Points(
    T *aVec,
    const size_t np,
    const unsigned int ndim,
    const T s) {
  const size_t n = np * ndim;
  for (unsigned int i = 0; i < n; i++) { aVec[i] *= s; }
}
#ifdef DFM2_STATIC_LIBRARY
template void delfem2::Scale_Points(
    float *, size_t, unsigned int, float);
template void delfem2::Scale_Points(
    double *, size_t, unsigned int, double);
#endif

// --------------

template<typename T>
void delfem2::Translate_Points(
    T *pVec,
    const size_t np,
    const unsigned int ndim,
    const T *trns) {
  for (unsigned int ip = 0; ip < np; ip++) {
    for (unsigned int idim = 0; idim < ndim; ++idim) {
      pVec[ip * ndim + idim] += trns[idim];
    }
  }
}
#ifdef DFM2_STATIC_LIBRARY
template void delfem2::Translate_Points(
    float *, size_t, unsigned int, const float *);
template void delfem2::Translate_Points(
    double *, size_t, unsigned int, const double *);
#endif

// --------------

template<typename T>
void delfem2::Translate_Points2(
    std::vector<T> &aXY,
    const T tx,
    const T ty) {
  const T trns[2] = {tx, ty};
  Translate_Points(aXY.data(), aXY.size() / 2, 2, trns);
}
#ifdef DFM2_STATIC_LIBRARY
template void delfem2::Translate_Points2(std::vector<float> &, float, float);
template void delfem2::Translate_Points2(std::vector<double> &, double, double);
#endif

// --------------

template<typename T>
void delfem2::Translate_Points3(
    std::vector<T> &aXYZ,
    const T tx,
    const T ty,
    const T tz) {
  const T trns[3] = {tx, ty, tz};
  Translate_Points(aXYZ.data(), aXYZ.size() / 3, 3, trns);
}
#ifdef DFM2_STATIC_LIBRARY
template void delfem2::Translate_Points3(
    std::vector<float> &, float, float, float);
template void delfem2::Translate_Points3(
    std::vector<double> &, double, double, double);
#endif

// --------------

template<typename T>
void delfem2::Rotate_Points3(
    std::vector<T> &aXYZ,
    const T radx,
    const T rady,
    const T radz) {
  T mat[9];
  points::Mat3_Bryant(mat, radx, rady, radz);
  T *pXYZ = aXYZ.data();
  const size_t nXYZ = aXYZ.size() / 3;
  for (unsigned int ixyz = 0; ixyz < nXYZ; ++ixyz) {
    const T p[3] = {aXYZ[ixyz * 3 + 0], aXYZ[ixyz * 3 + 1], aXYZ[ixyz * 3 + 2]};
    points::MatVec3(pXYZ + ixyz * 3, mat, p);
  }
}
#ifdef DFM2_STATIC_LIBRARY
template void delfem2::Rotate_Points3(
    std::vector<float> &, float, float, float);
template void delfem2::Rotate_Points3(
    std::vector<double> &, double, double, double);
#endif

// -----------------------------------------

double delfem2::Size_Points3D_LongestAABBEdge(
    const std::vector<double> &aXYZ) {
  double c[3], w[3];
  CenterWidth_Points3(c, w,
                      aXYZ);
  return points::largest(w[0], w[1], w[2]);
}


// ---------------------------------------

template<typename T>
DFM2_INLINE void delfem2::Normalize_Points3(
    std::vector<T> &aXYZ,
    T s) {
  T c[3], w[3];
  CenterWidth_Points3(
      c, w,
      aXYZ);
  Translate_Points3(
      aXYZ,
      -c[0], -c[1], -c[2]);
  T wmax = points::largest(w[0], w[1], w[2]);
  Scale_PointsX(
      aXYZ,
      s / wmax);
}
#ifdef DFM2_STATIC_LIBRARY
template void delfem2::Normalize_Points3(std::vector<float> &, float);
template void delfem2::Normalize_Points3(std::vector<double> &, double);
#endif


// ---------------------------------------

template<typename REAL>
DFM2_INLINE void delfem2::NormalizeVector_Points(
    REAL *aVec,
    unsigned int np,
    unsigned int ndim) {
  if (ndim == 2) {
    for (unsigned int ip = 0; ip < np; ++ip) {
      REAL *p = aVec + ip * 2;
      const REAL len = points::Length2(p);
      REAL linv = 1 / len;
      p[0] *= linv;
      p[1] *= linv;
    }
  } else if (ndim == 3) {
    for (unsigned int ip = 0; ip < np; ++ip) {
      REAL *p = aVec + ip * 3;
      const REAL len = points::Length3(p);
      const REAL linv = 1 / len;
      p[0] *= linv;
      p[1] *= linv;
      p[2] *= linv;
    }
  } else {
    assert(0);
  }
}
#ifdef DFM2_STATIC_LIBRARY
template void delfem2::NormalizeVector_Points(
    float *, unsigned int, unsigned int);
template void delfem2::NormalizeVector_Points(
    double *, unsigned int, unsigned int);
#endif

// --------------------------------------


double delfem2::EnergyKinetic(
    const double *aUVW,
    size_t np) {
  double E = 0.0;
  for (unsigned int ip = 0; ip < np; ++ip) {
    double u0 = aUVW[ip * 3 + 0];
    double v0 = aUVW[ip * 3 + 1];
    double w0 = aUVW[ip * 3 + 2];
    E += u0 * u0 + v0 * v0 + w0 * w0;
  }
  return E;
}

// ---------------------------------------

template<typename T>
void delfem2::CG_Point3(
    T *cg,
    const std::vector<T> &aXYZ) {
  cg[0] = cg[1] = cg[2] = 0;
  const size_t nXYZ = aXYZ.size() / 3;
  for (unsigned int ixyz = 0; ixyz < nXYZ; ixyz++) {
    cg[0] += aXYZ[ixyz * 3 + 0];
    cg[1] += aXYZ[ixyz * 3 + 1];
    cg[2] += aXYZ[ixyz * 3 + 2];
  }
  cg[0] /= nXYZ;
  cg[1] /= nXYZ;
  cg[2] /= nXYZ;
}
#ifdef DFM2_STATIC_LIBRARY
template void delfem2::CG_Point3(float *, const std::vector<float> &);
template void delfem2::CG_Point3(double *, const std::vector<double> &);
#endif


// -------------------------------------------------

template<typename T>
void delfem2::Points_RandomUniform(
    T *aCoords,
    size_t np,
    unsigned int ndim,
    const T *minCoords,
    const T *maxCoords) {
  std::random_device rd;
  std::mt19937 rdeng(rd());
  std::uniform_real_distribution<T> dist(0, 1);
  for (unsigned int ip = 0; ip < np; ++ip) {
    for (unsigned int idim = 0; idim < ndim; ++idim) {
      aCoords[ip * ndim + idim]
          = (maxCoords[idim] - minCoords[idim]) * dist(rdeng) + minCoords[idim];
    }
  }
}
#ifdef DFM2_STATIC_LIBRARY
template void delfem2::Points_RandomUniform(
    float *,
    size_t,
    unsigned int,
    const float *,
    const float *);
template void delfem2::Points_RandomUniform(
    double *,
    size_t,
    unsigned int,
    const double *,
    const double *);
#endif

// -------------------------------------------------

void delfem2::TangentVector_Points3(
    std::vector<double> &aOdir,
    const std::vector<double> &aNorm) {
  assert(aOdir.size() == aNorm.size());
  const size_t np = aNorm.size() / 3;
  aOdir.resize(np * 3);
  for (unsigned int ip = 0; ip < np; ++ip) {
    const double *n = aNorm.data() + ip * 3;
    const double *o0 = aOdir.data() + ip * 3;
    const double on = points::Dot3(o0, n);
    const double o1[3] = {
        o0[0] - on * n[0],
        o0[1] - on * n[1],
        o0[2] - on * n[2]};
    const double leninv = points::Length3(o1);
    aOdir[ip * 3 + 0] = o1[0] * leninv;
    aOdir[ip * 3 + 1] = o1[1] * leninv;
    aOdir[ip * 3 + 2] = o1[2] * leninv;
  }
}
