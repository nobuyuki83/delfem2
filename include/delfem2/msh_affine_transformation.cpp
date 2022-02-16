/*
 * Copyright (c) 2019 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */

#include "delfem2/msh_affine_transformation.h"

#include <cassert>
#include <cmath>
#include <vector>
#include <random>

#include "delfem2/msh_boundingbox.h"

#ifndef M_PI
#  define M_PI 3.14159265358979323846
#endif

// ------------------------------------------------

namespace delfem2::points {

//! @details we have "float" and "double" versions Length3 because of sqrtf and sqrt
template <typename T>
DFM2_INLINE T Length3(const T p[3]) {
  return std::sqrt(p[0] * p[0] + p[1] * p[1] + p[2] * p[2]);
}

//! @details we have "float" and "double" versions Length2 because of sqrtf and sqrt
template <typename T>
DFM2_INLINE T Length2(const T p[2]) {
  return std::sqrt(p[0] * p[0] + p[1] * p[1]);
}

template <typename T>
DFM2_INLINE void Cross3D(
    T r[3],
    const T v1[3],
    const T v2[3]) {
  r[0] = v1[1] * v2[2] - v2[1] * v1[2];
  r[1] = v1[2] * v2[0] - v2[2] * v1[0];
  r[2] = v1[0] * v2[1] - v2[0] * v1[1];
}

template <typename T>
DFM2_INLINE double TriArea2D(
    const T p0[2],
    const T p1[2],
    const T p2[2]) {
  return ( (p1[0] - p0[0]) * (p2[1] - p0[1]) - (p2[0] - p0[0]) * (p1[1] - p0[1]) )/2;
}

template <typename T>
DFM2_INLINE double TriArea3D(
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
    T n[3], T &a,
    const T v1[3], const T v2[3], const T v3[3]) {
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
DFM2_INLINE void MatVec3(
    T y[3],
    const T m[9], const T x[3]) {
  y[0] = m[0] * x[0] + m[1] * x[1] + m[2] * x[2];
  y[1] = m[3] * x[0] + m[4] * x[1] + m[5] * x[2];
  y[2] = m[6] * x[0] + m[7] * x[1] + m[8] * x[2];
}

//! @details we have "float" and "double" versions
//! Distance3 because of sqrtf and sqrt
template <typename T>
DFM2_INLINE T Distance3(const T p0[3], const T p1[3]) {
  return std::sqrt(
      (p1[0] - p0[0]) * (p1[0] - p0[0]) +
          (p1[1] - p0[1]) * (p1[1] - p0[1]) +
          (p1[2] - p0[2]) * (p1[2] - p0[2]));
}

template <typename T>
DFM2_INLINE T Distance2(const T p0[3], const T p1[3]) {
  return std::sqrt((p1[0] - p0[0]) * (p1[0] - p0[0]) + (p1[1] - p0[1]) * (p1[1] - p0[1]));
}

template <typename T>
DFM2_INLINE T Dot3(const T p0[3], const T p1[3]) {
  return p0[0]*p1[0] + p0[1]*p1[1] + p0[2]*p1[2];
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

template <typename T>
DFM2_INLINE void Mat3_Bryant(
    T m[9],
    T rx,
    T ry,
    T rz) {
  m[0] = std::cos(rz) * std::cos(ry);
  m[1] = std::cos(rz) * std::sin(ry) * std::sin(rx) - std::sin(rz) * std::cos(rx);
  m[2] = std::cos(rz) * std::sin(ry) * std::cos(rx) + std::sin(rz) * std::sin(rx);
  m[3] = std::sin(rz) * std::cos(ry);
  m[4] = std::sin(rz) * std::sin(ry) * std::sin(rx) + std::cos(rz) * std::cos(rx);
  m[5] = std::sin(rz) * std::sin(ry) * std::cos(rx) - std::cos(rz) * std::sin(rx);
  m[6] = -std::sin(ry);
  m[7] = std::cos(ry) * std::sin(rx);
  m[8] = std::cos(ry) * std::cos(rx);
}

}

// static function above
// ==============================================
// exposed function below


// -------------------------------------

template<typename T>
void delfem2::Scale_PointsX(
    std::vector<T> &vtx_xyz,
    const T scale) {
  const std::size_t n = vtx_xyz.size();
  for (unsigned int i = 0; i < n; ++i) { vtx_xyz[i] *= scale; }
}
#ifdef DFM2_STATIC_LIBRARY
template void delfem2::Scale_PointsX(std::vector<float> &, float);
template void delfem2::Scale_PointsX(std::vector<double> &, double);
#endif

// --------------

template<typename T>
void delfem2::Scale_Points(
    T *vtx_coords,
    const size_t num_vtx,
    const unsigned int ndim,
    const T scale) {
  const size_t n = num_vtx * ndim;
  for (unsigned int i = 0; i < n; i++) { vtx_coords[i] *= scale; }
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
    T *vtx_coords,
    const size_t num_vtx,
    const unsigned int ndim,
    const T *translation) {
  for (unsigned int ip = 0; ip < num_vtx; ip++) {
    for (unsigned int idim = 0; idim < ndim; ++idim) {
      vtx_coords[ip * ndim + idim] += translation[idim];
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

template<typename T>
DFM2_INLINE void delfem2::Normalize_Points3(
    std::vector<T> &vtx_xyz,
    T length_longest_edge_boundingbox) {
  T c[3], w[3];
  CenterWidth_Points3(
      c, w,
      vtx_xyz);
  Translate_Points3(
      vtx_xyz,
      -c[0], -c[1], -c[2]);
  T wmax = points::largest(w[0], w[1], w[2]);
  Scale_PointsX(
      vtx_xyz,
      length_longest_edge_boundingbox / wmax);
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




// ---------------------------------------


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
