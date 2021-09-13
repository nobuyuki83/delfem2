/*
 * Copyright (c) 2019 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */

#include "delfem2/mat3.h"

#include <random>

namespace delfem2 {
namespace mat3 {

// https://en.wikipedia.org/wiki/Quaternions_and_spatial_rotation
// row major matrix
template<typename T>
DFM2_INLINE void SetMatrix3_Quaternion
    (T r[],
     const T q[]) {
  const T x2 = q[0] * q[0] * 2;
  const T y2 = q[1] * q[1] * 2;
  const T z2 = q[2] * q[2] * 2;
  const T xy = q[0] * q[1] * 2;
  const T yz = q[1] * q[2] * 2;
  const T zx = q[2] * q[0] * 2;
  const T xw = q[0] * q[3] * 2;
  const T yw = q[1] * q[3] * 2;
  const T zw = q[2] * q[3] * 2;
  r[0] = 1 - y2 - z2;
  r[1] = xy - zw;
  r[2] = zx + yw;
  r[3] = xy + zw;
  r[4] = 1 - z2 - x2;
  r[5] = yz - xw;
  r[6] = zx - yw;
  r[7] = yz + xw;
  r[8] = 1 - x2 - y2;
}
#ifdef DFM2_STATIC_LIBRARY
template void SetMatrix3_Quaternion(float r[], const float q[]);
template void SetMatrix3_Quaternion(double r[], const double q[]);
#endif

DFM2_INLINE double estimationMaxEigenValue(const double mtm[6]) {
  double maxl = 1;
  {  // estimation of maximum eigen value using Gerschgorin's circle theorem
    maxl = mtm[0] + fabs(mtm[3]) + fabs(mtm[5]);
    const double tmp2 = mtm[1] + fabs(mtm[3]) + fabs(mtm[4]);
    maxl = (tmp2 > maxl) ? tmp2 : maxl;
    const double tmp3 = mtm[2] + fabs(mtm[5]) + fabs(mtm[4]);
    maxl = (tmp3 > maxl) ? tmp3 : maxl;
  }
  return maxl;
}

DFM2_INLINE void Cross3(double r[3], const double v1[3], const double v2[3]) {
  r[0] = v1[1] * v2[2] - v2[1] * v1[2];
  r[1] = v1[2] * v2[0] - v2[2] * v1[0];
  r[2] = v1[0] * v2[1] - v2[0] * v1[1];
}

DFM2_INLINE void Normalize3(double v[3]) {
  double len = sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);
  v[0] /= len;
  v[1] /= len;
  v[2] /= len;
}

DFM2_INLINE double SqLength3(const double v[3]) {
  return v[0] * v[0] + v[1] * v[1] + v[2] * v[2];
}

}
}

// ------------------------

// t is a tmporary buffer size of 9
template<typename T>
void delfem2::Transpose_Mat3(
    T t[9],
    const T a[9]) {
  t[0] = a[0];
  t[1] = a[3];
  t[2] = a[6];
  t[3] = a[1];
  t[4] = a[4];
  t[5] = a[7];
  t[6] = a[2];
  t[7] = a[5];
  t[8] = a[8];
}
#ifdef DFM2_STATIC_LIBRARY
template void delfem2::Transpose_Mat3(float t[], const float a[]);
template void delfem2::Transpose_Mat3(double t[], const double a[]);
#endif

// --------------------------------------

template<typename T0, typename T1>
void delfem2::Inverse_Mat3(
    T0 Ainv[9],
    const T1 A[9]) {
  const T0 det =
      +A[0] * A[4] * A[8] + A[3] * A[7] * A[2] + A[6] * A[1] * A[5]
          - A[0] * A[7] * A[5] - A[6] * A[4] * A[2] - A[3] * A[1] * A[8];
  const T0 inv_det = 1 / det;
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
#ifdef DFM2_STATIC_LIBRARY
template void delfem2::Inverse_Mat3(float Ainv[9], const float A[9]);
template void delfem2::Inverse_Mat3(double Ainv[9], const double A[9]);
#endif


// --------------------------------------

template<typename REAL>
void delfem2::Inverse_Mat3(
    REAL A[9]) {
  const REAL B[9] = {
      A[0], A[1], A[2],
      A[3], A[4], A[5],
      A[6], A[7], A[8]};
  const REAL det =
      +B[0] * B[4] * B[8] + B[3] * B[7] * B[2] + B[6] * B[1] * B[5]
          - B[0] * B[7] * B[5] - B[6] * B[4] * B[2] - B[3] * B[1] * B[8];
  const REAL inv_det = 1 / det;
  A[0] = inv_det * (B[4] * B[8] - B[5] * B[7]);
  A[1] = inv_det * (B[2] * B[7] - B[1] * B[8]);
  A[2] = inv_det * (B[1] * B[5] - B[2] * B[4]);
  A[3] = inv_det * (B[5] * B[6] - B[3] * B[8]);
  A[4] = inv_det * (B[0] * B[8] - B[2] * B[6]);
  A[5] = inv_det * (B[2] * B[3] - B[0] * B[5]);
  A[6] = inv_det * (B[3] * B[7] - B[4] * B[6]);
  A[7] = inv_det * (B[1] * B[6] - B[0] * B[7]);
  A[8] = inv_det * (B[0] * B[4] - B[1] * B[3]);
}
#ifdef DFM2_STATIC_LIBRARY
template void delfem2::Inverse_Mat3(float Ainv[9]);
template void delfem2::Inverse_Mat3(double Ainv[9]);
#endif

// -----------------------------------

// see https://en.wikipedia.org/wiki/Skew-symmetric_matrix
template<typename REAL>
void delfem2::Mat3_Spin(
    REAL *mat,
    const REAL *v) {
  mat[0] = 0;
  mat[1] = -v[2];
  mat[2] = +v[1];
  mat[3] = +v[2];
  mat[4] = 0;
  mat[5] = -v[0];
  mat[6] = -v[1];
  mat[7] = +v[0];
  mat[8] = 0;
}
#ifdef DFM2_STATIC_LIBRARY
template void delfem2::Mat3_Spin(float *mat, const float *v);
template void delfem2::Mat3_Spin(double *mat, const double *v);
#endif

// ---------------------------------------------------------

template<typename REAL>
void delfem2::Mat3_Spin_ScaleAdd(
    REAL *m,
    const REAL *v,
    REAL alpha, REAL beta) {
  m[0] = beta * m[0];
  m[1] = beta * m[1] - v[2] * alpha;
  m[2] = beta * m[2] + v[1] * alpha;
  m[3] = beta * m[3] + v[2] * alpha;
  m[4] = beta * m[4];
  m[5] = beta * m[5] - v[0] * alpha;
  m[6] = beta * m[6] - v[1] * alpha;
  m[7] = beta * m[7] + v[0] * alpha;
  m[8] = beta * m[8];
}
#ifdef DFM2_STATIC_LIBRARY
template void delfem2::Mat3_Spin_ScaleAdd(float *mat, const float *v, float alpha, float beta);
template void delfem2::Mat3_Spin_ScaleAdd(double *mat, const double *v, double alpha, double beta);
#endif

// ---------------------------------------------------------

template<typename REAL>
void delfem2::Mat3_Identity_ScaleAdd(
    REAL *mat,
    REAL alpha, REAL beta) {
  mat[0] = beta * mat[0] + alpha;
  mat[1] = beta * mat[1];
  mat[2] = beta * mat[2];
  mat[3] = beta * mat[3];
  mat[4] = beta * mat[4] + alpha;
  mat[5] = beta * mat[5];
  mat[6] = beta * mat[6];
  mat[7] = beta * mat[7];
  mat[8] = beta * mat[8] + alpha;
}
#ifdef DFM2_STATIC_LIBRARY
template void delfem2::Mat3_Identity_ScaleAdd(float *mat, float alpha, float beta);
template void delfem2::Mat3_Identity_ScaleAdd(double *mat, double alpha, double beta);
#endif

template<typename REAL>
void delfem2::Mat3_Identity(REAL *mat,
                            REAL alpha) {
  mat[0] = alpha;
  mat[1] = 0;
  mat[2] = 0;
  mat[3] = 0;
  mat[4] = alpha;
  mat[5] = 0;
  mat[6] = 0;
  mat[7] = 0;
  mat[8] = alpha;
}
#ifdef DFM2_STATIC_LIBRARY
template void delfem2::Mat3_Identity(float *mat, float alpha);
template void delfem2::Mat3_Identity(double *mat, double alpha);
#endif


// -------------

namespace delfem2 {

template<>
DFM2_INLINE void Mat3_AffineRotation(
    double *mat,
    double theta) {
  mat[0] = +cos(theta);
  mat[1] = -sin(theta);
  mat[2] = 0;
  mat[3] = +sin(theta);
  mat[4] = +cos(theta);
  mat[5] = 0;
  mat[6] = 0;
  mat[7] = 0;
  mat[8] = 1;
}

template<>
DFM2_INLINE void Mat3_AffineRotation(
    float *mat,
    float theta) {
  mat[0] = +cosf(theta);
  mat[1] = -sinf(theta);
  mat[2] = 0;
  mat[3] = +sinf(theta);
  mat[4] = +cosf(theta);
  mat[5] = 0;
  mat[6] = 0;
  mat[7] = 0;
  mat[8] = 1;
}

}

// --------------------

template<typename T>
void delfem2::Mat3_AffineTranslation(
    T *mat,
    const T transl[2]) {
  mat[0] = 1;
  mat[1] = 0;
  mat[2] = transl[0];
  mat[3] = 0;
  mat[4] = 1;
  mat[5] = transl[1];
  mat[6] = 0;
  mat[7] = 0;
  mat[8] = 1;
}
#ifdef DFM2_STATIC_LIBRARY
template void delfem2::Mat3_AffineTranslation(float *, const float [2]);
template void delfem2::Mat3_AffineTranslation(double *, const double [2]);
#endif

// --------------------------------------------------------


template<typename T0, typename T1, typename T2>
void delfem2::MatVec3(
    T0 y[3],
    const T1 m[9],
    const T2 x[3]) {
  y[0] = m[0] * x[0] + m[1] * x[1] + m[2] * x[2];
  y[1] = m[3] * x[0] + m[4] * x[1] + m[5] * x[2];
  y[2] = m[6] * x[0] + m[7] * x[1] + m[8] * x[2];
}
#ifdef DFM2_STATIC_LIBRARY
template void delfem2::MatVec3(float y[3], const float m[9], const float x[3]);
template void delfem2::MatVec3(double y[3], const double m[9], const double x[3]);
#endif

// ----------------------------

template<typename T>
void delfem2::MatVec3_ScaleAdd(
    T y[3],
    const T m[9], const T x[3],
    T alpha, T beta) {
  y[0] = beta * y[0] + alpha * (m[0] * x[0] + m[1] * x[1] + m[2] * x[2]);
  y[1] = beta * y[1] + alpha * (m[3] * x[0] + m[4] * x[1] + m[5] * x[2]);
  y[2] = beta * y[2] + alpha * (m[6] * x[0] + m[7] * x[1] + m[8] * x[2]);
}
#ifdef DFM2_STATIC_LIBRARY
template void delfem2::MatVec3_ScaleAdd(
    float y[3],
    const float m[9], const float x[3], float alpha, float beta);
template void delfem2::MatVec3_ScaleAdd(
    double y[3],
    const double m[9], const double x[3], double alpha, double beta);
#endif

DFM2_INLINE void delfem2::VecMat3
    (double y[3],
     const double x[3], const double m[9]) {
  y[0] = m[0] * x[0] + m[3] * x[1] + m[6] * x[2];
  y[1] = m[1] * x[0] + m[4] * x[1] + m[7] * x[2];
  y[2] = m[2] * x[0] + m[5] * x[1] + m[8] * x[2];
}

template<typename T0, typename T1, typename T2>
DFM2_INLINE void delfem2::MatTVec3(
    T0 y[3],
    const T1 m[9],
    const T2 x[3]) {
  y[0] = m[0] * x[0] + m[3] * x[1] + m[6] * x[2];
  y[1] = m[1] * x[0] + m[4] * x[1] + m[7] * x[2];
  y[2] = m[2] * x[0] + m[5] * x[1] + m[8] * x[2];
}
#ifdef DFM2_STATIC_LIBRARY
template void delfem2::MatTVec3(float y[3], const float m[9], const float x[3]);
template void delfem2::MatTVec3(double y[3], const double m[9], const double x[3]);
#endif

template<typename T0, typename T1, typename T2, typename T3, typename T4>
void delfem2::MatTVec3_ScaleAdd(
    T0 y[3],
    const T1 m[9], const T2 x[3],
    T3 alpha, T4 beta) {
  y[0] = beta * y[0] + alpha * (m[0] * x[0] + m[3] * x[1] + m[6] * x[2]);
  y[1] = beta * y[1] + alpha * (m[1] * x[0] + m[4] * x[1] + m[7] * x[2]);
  y[2] = beta * y[2] + alpha * (m[2] * x[0] + m[5] * x[1] + m[8] * x[2]);
}
#ifdef DFM2_STATIC_LIBRARY
template void delfem2::MatTVec3_ScaleAdd
    (float y[3],
     const float m[9], const float x[3], float alpha, float beta);
template void delfem2::MatTVec3_ScaleAdd
    (double y[3],
     const double m[9], const double x[3], double alpha, double beta);
#endif

template<typename T>
DFM2_INLINE void delfem2::Vec2_Mat3Vec2_AffineProjection(
    T y[2],
    const T A[9],
    const T x[2]) {
  y[0] = A[0] * x[0] + A[1] * x[1] + A[2];
  y[1] = A[3] * x[0] + A[4] * x[1] + A[5];
  const T w = A[6] * x[0] + A[7] * x[1] + A[8];
  y[0] /= w;
  y[1] /= w;
}
#ifdef DFM2_STATIC_LIBRARY
template void delfem2::Vec2_Mat3Vec2_AffineProjection(float [2], const float [9], const float [2]);
template void delfem2::Vec2_Mat3Vec2_AffineProjection(double [2], const double [9], const double [2]);
#endif

// ---------


template<typename T>
DFM2_INLINE void delfem2::Vec2_Mat3Vec2_AffineDirection(
    T y[2],
    const T A[9],
    const T x[2]) {
  y[0] = A[0] * x[0] + A[1] * x[1];
  y[1] = A[3] * x[0] + A[4] * x[1];
}
#ifdef DFM2_STATIC_LIBRARY
template void delfem2::Vec2_Mat3Vec2_AffineDirection(float [2], const float [9], const float [2]);
template void delfem2::Vec2_Mat3Vec2_AffineDirection(double [2], const double [9], const double [2]);
#endif

// above: mat3 and vec3
// ------------------------------------
// below: mat3 and quaternion

template<typename REAL>
DFM2_INLINE void delfem2::Mat3_Quat(
    REAL r[],
    const REAL q[]) {
  const REAL x2 = q[0] * q[0] * 2;
  const REAL y2 = q[1] * q[1] * 2;
  const REAL z2 = q[2] * q[2] * 2;
  const REAL xy = q[0] * q[1] * 2;
  const REAL yz = q[1] * q[2] * 2;
  const REAL zx = q[2] * q[0] * 2;
  const REAL xw = q[0] * q[3] * 2;
  const REAL yw = q[1] * q[3] * 2;
  const REAL zw = q[2] * q[3] * 2;
  r[0] = 1 - y2 - z2;
  r[1] = xy - zw;
  r[2] = zx + yw;
  r[3] = xy + zw;
  r[4] = 1 - z2 - x2;
  r[5] = yz - xw;
  r[6] = zx - yw;
  r[7] = yz + xw;
  r[8] = 1 - x2 - y2;
}
#ifdef DFM2_STATIC_LIBRARY
template void delfem2::Mat3_Quat(float r[], const float q[]);
template void delfem2::Mat3_Quat(double r[], const double q[]);
#endif

// ----------------------------------------------

template<typename T0, typename T1, typename T2>
void delfem2::MatMat3(
    T0 *AB,
    const T1 *A,
    const T2 *B) {
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      AB[i * 3 + j] =
          A[i * 3 + 0] * B[0 * 3 + j] +
          A[i * 3 + 1] * B[1 * 3 + j] +
          A[i * 3 + 2] * B[2 * 3 + j];
    }
  }
}
#ifdef DFM2_STATIC_LIBRARY
template void delfem2::MatMat3(float *C, const float *A, const float *B);
template void delfem2::MatMat3(double *C, const double *A, const double *B);
#endif

// ---------------------------------------

template<typename T0, typename T1, typename T2>
void delfem2::MatMatT3(
    T0 *ABt,
    const T1 *A,
    const T2 *B) {
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      ABt[i * 3 + j] =
          A[i * 3 + 0] * B[j * 3 + 0] +
          A[i * 3 + 1] * B[j * 3 + 1] +
          A[i * 3 + 2] * B[j * 3 + 2];
    }
  }
}
#ifdef DFM2_STATIC_LIBRARY
template void delfem2::MatMatT3(float *C, const float *A, const float *B);
template void delfem2::MatMatT3(double *C, const double *A, const double *B);
#endif

// ---------------------------------------

template<typename T0, typename T1, typename T2>
void delfem2::MatTMat3(
    T0 *AtB,
    const T1 *A, const T2 *B) {
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      AtB[i * 3 + j] =
          A[0 * 3 + i] * B[0 * 3 + j] +
          A[1 * 3 + i] * B[1 * 3 + j] +
          A[2 * 3 + i] * B[2 * 3 + j];
    }
  }
}
#ifdef DFM2_STATIC_LIBRARY
template void delfem2::MatTMat3(float *C, const float *A, const float *B);
template void delfem2::MatTMat3(double *C, const double *A, const double *B);
#endif

// --------------------------------------

template<typename T>
void delfem2::MatTMat3_ScaleAdd(
    T *C,
    const T *A, const T *B,
    T alpha, T beta) {
  for (int i = 0; i < 9; ++i) { C[i] *= beta; }
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      C[i * 3 + j] += alpha * (A[0 * 3 + i] * B[0 * 3 + j] + A[1 * 3 + i] * B[1 * 3 + j] + A[2 * 3 + i] * B[2 * 3 + j]);
    }
  }
}
#ifdef DFM2_STATIC_LIBRARY
template void delfem2::MatTMat3_ScaleAdd
    (float *C,
     const float *A, const float *B, float alpha, float beta);
template void delfem2::MatTMat3_ScaleAdd
    (double *C,
     const double *A, const double *B, double alpha, double beta);
#endif

// ---------------------------------------

template<typename T>
T delfem2::Det_Mat3(const T U[9]) {
  return +U[0] * U[4] * U[8] + U[3] * U[7] * U[2] + U[6] * U[1] * U[5]
      - U[0] * U[7] * U[5] - U[6] * U[4] * U[2] - U[3] * U[1] * U[8];
}
#ifdef DFM2_STATIC_LIBRARY
template float delfem2::Det_Mat3(const float U[9]);
template double delfem2::Det_Mat3(const double U[9]);
#endif

// ---------------------------------------

template<typename T>
T delfem2::SquareNormFrobenius_SymMat3
    (const T sm[6]) {
  return sm[0] * sm[0] + sm[1] * sm[1] + sm[2] * sm[2] + 2 * (sm[3] * sm[3] + sm[4] * sm[4] + sm[5] * sm[5]);
}
#ifdef DFM2_STATIC_LIBRARY
template float delfem2::SquareNormFrobenius_SymMat3(const float sm[6]);
template double delfem2::SquareNormFrobenius_SymMat3(const double sm[6]);
#endif

// ---------------------------------------

template<typename REAL>
void delfem2::Mat3_Rotation_Cartesian(
    REAL mat[9],
    const REAL vec[3]) {
  REAL sqt = vec[0] * vec[0] + vec[1] * vec[1] + vec[2] * vec[2];
  if (sqt < 1.0e-20) { // infinitesmal rotation approximation
    mat[0] = 1;
    mat[1] = -vec[2];
    mat[2] = +vec[1];
    mat[3] = +vec[2];
    mat[4] = 1;
    mat[5] = -vec[0];
    mat[6] = -vec[1];
    mat[7] = +vec[0];
    mat[8] = 1;
    return;
  }
  REAL t = std::sqrt(sqt);
  REAL invt = 1 / t;
  REAL n[3] = {vec[0] * invt, vec[1] * invt, vec[2] * invt};
  const REAL c0 = std::cos(t);
  const REAL s0 = std::sin(t);
  mat[0 * 3 + 0] = c0 + (1 - c0) * n[0] * n[0];
  mat[0 * 3 + 1] = -n[2] * s0 + (1 - c0) * n[0] * n[1];
  mat[0 * 3 + 2] = +n[1] * s0 + (1 - c0) * n[0] * n[2];
  mat[1 * 3 + 0] = +n[2] * s0 + (1 - c0) * n[1] * n[0];
  mat[1 * 3 + 1] = c0 + (1 - c0) * n[1] * n[1];
  mat[1 * 3 + 2] = -n[0] * s0 + (1 - c0) * n[1] * n[2];
  mat[2 * 3 + 0] = -n[1] * s0 + (1 - c0) * n[2] * n[0];
  mat[2 * 3 + 1] = +n[0] * s0 + (1 - c0) * n[2] * n[1];
  mat[2 * 3 + 2] = c0 + (1 - c0) * n[2] * n[2];
}
#ifdef DFM2_STATIC_LIBRARY
template void delfem2::Mat3_Rotation_Cartesian(float mat[9], const float vec[3]);
template void delfem2::Mat3_Rotation_Cartesian(double mat[9], const double vec[3]);
#endif


// ---------------------------------------

// compute eigen value & vector for symmmetric matrix
// sm[6] = (M_00,M_11,M_22,M_12,M_20,M_01)
// M = ULU^T
// u[9] = (U_00,U_01,U_02, U_10,U_11,U_12, U_20,U_21,U_22)
DFM2_INLINE bool delfem2::eigenSym3(
    double u[9],
    double l[3],
    const double sm[6],
    int nitr) {
  // initialize u as identity matrix
  u[0] = u[4] = u[8] = 1.0;
  u[1] = u[2] = u[3] = u[5] = u[6] = u[7] = 0.0;
  l[0] = l[1] = l[2] = 0.0;
  double dnrm = SquareNormFrobenius_SymMat3(sm);
  if (dnrm < 1.0e-30) { return false; } // this matrix is too small
  const double scale = sqrt(dnrm);
  const double invscl = 1.0 / scale;
  double sms[6] = {sm[0] * invscl, sm[1] * invscl, sm[2] * invscl, sm[3] * invscl, sm[4] * invscl, sm[5] * invscl};
  for (int itr = 0; itr < nitr; itr++) {
    const double m[6] = {sms[0], sms[1], sms[2], sms[3], sms[4], sms[5]};
    const double v[9] = {u[0], u[1], u[2], u[3], u[4], u[5], u[6], u[7], u[8]};
    const double a12 = fabs(sms[3]);
    const double a20 = fabs(sms[4]);
    const double a01 = fabs(sms[5]);
    if (a12 >= a20 && a12 >= a01) {
      // when a12 sms[3] is the biggest
      const double t = 0.5 * atan2(2 * m[3], m[2] - m[1]);
      const double ct = cos(t);
      const double st = sin(t);
      sms[1] = ct * ct * m[1] + st * st * m[2] - 2 * st * ct * m[3];
      sms[2] = ct * ct * m[2] + st * st * m[1] + 2 * st * ct * m[3];
      sms[3] = 0; // (ct*ct-st*st)*m[3]+st*ct*(m[1]-m[2]);
      sms[4] = st * m[5] + ct * m[4];
      sms[5] = ct * m[5] - st * m[4];
      //
      u[1] = +ct * v[1] - st * v[2];
      u[2] = +st * v[1] + ct * v[2];
      u[4] = +ct * v[4] - st * v[5];
      u[5] = +st * v[4] + ct * v[5];
      u[7] = +ct * v[7] - st * v[8];
      u[8] = +st * v[7] + ct * v[8];
    } else if (a20 >= a01 && a20 >= a12) {
      // when a20 sms[4] is the biggest
      // the above condition statement shoud pass exactly once for each iteration.
      const double t = 0.5 * atan2(2 * m[4], m[2] - m[0]);
      const double ct = cos(t);
      const double st = sin(t);
      sms[0] = ct * ct * m[0] + st * st * m[2] - 2 * st * ct * m[4];
      sms[2] = ct * ct * m[2] + st * st * m[0] + 2 * st * ct * m[4];
      sms[3] = st * m[5] + ct * m[3];
      sms[4] = 0; // (ct*ct-st*st)*m[4]+st*ct*(m[0]-m[2]);
      sms[5] = ct * m[5] - st * m[3];
      //
      u[0] = +ct * v[0] - st * v[2];
      u[2] = +st * v[0] + ct * v[2];
      u[3] = +ct * v[3] - st * v[5];
      u[5] = +st * v[3] + ct * v[5];
      u[6] = +ct * v[6] - st * v[8];
      u[8] = +st * v[6] + ct * v[8];
    } else {
      // when a01 sms[5] is the biggest
      // the condition statement shoud pass exactly once for each iteration.
      const double t = 0.5 * atan2(2 * m[5], m[1] - m[0]);
      const double ct = cos(t);
      const double st = sin(t);
      sms[0] = ct * ct * m[0] + st * st * m[1] - 2 * st * ct * m[5];
      sms[1] = ct * ct * m[1] + st * st * m[0] + 2 * st * ct * m[5];
      sms[3] = st * m[4] + ct * m[3];
      sms[4] = ct * m[4] - st * m[3];
      sms[5] = 0; // (ct*ct-st*st)*m[5]+st*ct*(m[0]-m[1]);
      //
      u[0] = +ct * v[0] - st * v[1];
      u[1] = +st * v[0] + ct * v[1];
      u[3] = +ct * v[3] - st * v[4];
      u[4] = +st * v[3] + ct * v[4];
      u[6] = +ct * v[6] - st * v[7];
      u[7] = +st * v[6] + ct * v[7];
    }
  }
  l[0] = scale * sms[0];
  l[1] = scale * sms[1];
  l[2] = scale * sms[2];
  return true;
}

DFM2_INLINE void SortEigen3
    (double G[3], double V[9]) {
  double t;
  if (G[1] > G[0]) {  // swap 01
    t = G[0];
    G[0] = G[1];
    G[1] = t;
    t = V[0];
    V[0] = V[1];
    V[1] = t;
    t = V[3];
    V[3] = V[4];
    V[4] = t;
    t = V[6];
    V[6] = V[7];
    V[7] = t;
  }
  if (G[2] > G[1]) {
    t = G[1];
    G[1] = G[2];
    G[2] = t;
    t = V[1];
    V[1] = V[2];
    V[2] = t;
    t = V[4];
    V[4] = V[5];
    V[5] = t;
    t = V[7];
    V[7] = V[8];
    V[8] = t;
  }
  if (G[1] > G[0]) { // swap 01
    t = G[0];
    G[0] = G[1];
    G[1] = t;
    t = V[0];
    V[0] = V[1];
    V[1] = t;
    t = V[3];
    V[3] = V[4];
    V[4] = t;
    t = V[6];
    V[6] = V[7];
    V[7] = t;
  }
}

// m = UGV^T
DFM2_INLINE void delfem2::svd3(
    double U[9],
    double G[3],
    double V[9],
    const double m[9],
    int nitr) {
  // M^TM = VGGV^T
  const double mtm[6] = {
      m[0] * m[0] + m[3] * m[3] + m[6] * m[6],
      m[1] * m[1] + m[4] * m[4] + m[7] * m[7],
      m[2] * m[2] + m[5] * m[5] + m[8] * m[8],
      m[1] * m[2] + m[4] * m[5] + m[7] * m[8],
      m[2] * m[0] + m[5] * m[3] + m[8] * m[6],
      m[0] * m[1] + m[3] * m[4] + m[6] * m[7]};
  double lv[3];
  eigenSym3(V, lv,
            mtm, nitr);
  SortEigen3(lv, V);
  if (lv[0] < 0) { lv[0] = 0.0; }
  if (lv[1] < 0) { lv[1] = 0.0; }
  if (lv[2] < 0) { lv[2] = 0.0; }
  G[0] = sqrt(lv[0]);
  G[1] = sqrt(lv[1]);
  G[2] = sqrt(lv[2]);

  double u0[3] = {
      m[0] * V[0] + m[1] * V[3] + m[2] * V[6],
      m[3] * V[0] + m[4] * V[3] + m[5] * V[6],
      m[6] * V[0] + m[7] * V[3] + m[8] * V[6]};
  double u1[3] = {
      m[0] * V[1] + m[1] * V[4] + m[2] * V[7],
      m[3] * V[1] + m[4] * V[4] + m[5] * V[7],
      m[6] * V[1] + m[7] * V[4] + m[8] * V[7]};
  double u2[3] = {
      m[0] * V[2] + m[1] * V[5] + m[2] * V[8],
      m[3] * V[2] + m[4] * V[5] + m[5] * V[8],
      m[6] * V[2] + m[7] * V[5] + m[8] * V[8]};

  if (Det_Mat3(V) < 0) {  // making right hand coordinate
    V[0 * 3 + 2] *= -1;
    V[1 * 3 + 2] *= -1;
    V[2 * 3 + 2] *= -1;
    G[2] *= -1;
    U[0 * 3 + 2] *= -1;
    U[1 * 3 + 2] *= -1;
    U[2 * 3 + 2] *= -1;
  }

  const double sql0 = mat3::SqLength3(u0);
  if (sql0 > 1.0e-20) {
    mat3::Normalize3(u0);
    const double sql1 = mat3::SqLength3(u1);
    if (sql1 < 1.0e-20) {
      u1[0] = 1.0 - fabs(u0[0]);
      u1[1] = 1.0 - fabs(u0[1]);
      u1[2] = 1.0 - fabs(u0[2]);
    } else {
      mat3::Normalize3(u1);
    }
    const double d01 = u0[0] * u1[0] + u0[1] * u1[1] + u0[2] * u1[2];
    u1[0] -= d01 * u0[0];
    u1[1] -= d01 * u0[1];
    u1[2] -= d01 * u0[2];
    mat3::Normalize3(u1);
    double s2[3];
    mat3::Cross3(s2, u0, u1);
    const double d22 = u2[0] * s2[0] + u2[1] * s2[1] + u2[2] * s2[2];
    u2[0] = s2[0];
    u2[1] = s2[1];
    u2[2] = s2[2];
    if (d22 < 0) {
      G[2] *= -1;
    }
  } else {
    u0[0] = 1;
    u0[1] = 0;
    u0[2] = 0;
    u1[0] = 0;
    u1[1] = 1;
    u1[2] = 0;
    u2[0] = 0;
    u2[1] = 0;
    u2[2] = 1;
  }
  U[0] = u0[0];
  U[1] = u1[0];
  U[2] = u2[0];
  U[3] = u0[1];
  U[4] = u1[1];
  U[5] = u2[1];
  U[6] = u0[2];
  U[7] = u1[2];
  U[8] = u2[2];
}

DFM2_INLINE void delfem2::GetRotPolarDecomp(
    double R[9],
    //
    const double am[9],
    int nitr) {
  double U[9], G[3], V[9];
  // am = UGV^T
  svd3(U, G, V,
       am, nitr);
  // R = UV^T
  MatMatT3(R,
           U, V);
}

// -----------------------------------

// https://en.wikipedia.org/wiki/Axis%E2%80%93angle_representation
template<typename T>
DFM2_INLINE void delfem2::AxisAngleVectorCartesian_Mat3(
    T v[3],
    const T m[9]) {
  const T cos_t0 = (m[0] + m[4] + m[8] - 1) / 2;
  if (std::fabs(cos_t0 - 1) < 1.0e-5) {  // very small rotation
    v[0] = (m[3 * 2 + 1] - m[3 * 1 + 2]) / 2;
    v[1] = (m[3 * 0 + 2] - m[3 * 2 + 0]) / 2;
    v[2] = (m[3 * 1 + 0] - m[3 * 0 + 1]) / 2;
    return;
  }
  const T t0 = std::acos(cos_t0);
  const T c0 = t0 / (2 * std::sin(t0));
  v[0] = c0 * (m[3 * 2 + 1] - m[3 * 1 + 2]);
  v[1] = c0 * (m[3 * 0 + 2] - m[3 * 2 + 0]);
  v[2] = c0 * (m[3 * 1 + 0] - m[3 * 0 + 1]);
}

template<typename T>
DFM2_INLINE void delfem2::AxisAngleVectorCRV_Mat3(
    T crv[3],
    const T mat[9]) {
  T eparam2[4];
  const CMat3<T> m(mat);
  m.GetQuat_RotMatrix(eparam2);
  crv[0] = 4 * eparam2[0] / (1 + eparam2[3]);
  crv[1] = 4 * eparam2[1] / (1 + eparam2[3]);
  crv[2] = 4 * eparam2[2] / (1 + eparam2[3]);
}

// -----------------------------------

namespace delfem2 {

template<typename T>
CMat3<T> operator*(double d, const CMat3<T> &rhs) {
  CMat3<T> temp = rhs;
  temp *= d;
  return temp;
}
#ifdef DFM2_STATIC_LIBRARY
template CMat3<double> operator*(double d, const CMat3<double> &rhs);
#endif

template<typename T>
CMat3<T> operator*(const CMat3<T> &m, T d) {
  CMat3<T> t = m;
  t *= d;
  return t;
}
#ifdef DFM2_STATIC_LIBRARY
template CMat3<float> operator*(const CMat3<float> &m, float d);
template CMat3<double> operator*(const CMat3<double> &m, double d);
#endif

// --------------------

template<typename T>
CMat3<T> operator/(const CMat3<T> &m, T d) {
  CMat3<T> temp = m;
  temp /= d;
  return temp;
}
#ifdef DFM2_STATIC_LIBRARY
template CMat3<float> operator/(const CMat3<float> &m, float d);
template CMat3<double> operator/(const CMat3<double> &m, double d);
#endif

// ----------------------

template<typename T>
CMat3<T> operator+(const CMat3<T> &lhs, const CMat3<T> &rhs) {
  CMat3<T> temp = lhs;
  temp += rhs;
  return temp;
}
#ifdef DFM2_STATIC_LIBRARY
template CMat3<float> operator+(const CMat3<float> &lhs, const CMat3<float> &rhs);
template CMat3<double> operator+(const CMat3<double> &lhs, const CMat3<double> &rhs);
#endif

// ------------------

template<typename T>
CMat3<T> operator*(const CMat3<T> &lhs, const CMat3<T> &rhs) {
  return lhs.MatMat(rhs);
}
#ifdef DFM2_STATIC_LIBRARY
template CMat3<float> operator*(const CMat3<float> &lhs, const CMat3<float> &rhs);
template CMat3<double> operator*(const CMat3<double> &lhs, const CMat3<double> &rhs);
#endif

// ------------------------------

template<typename T>
CMat3<T> operator-(const CMat3<T> &lhs, const CMat3<T> &rhs) {
  CMat3<T> temp = lhs;
  temp -= rhs;
  return temp;
}
#ifdef DFM2_STATIC_LIBRARY
template CMat3<double> operator-(const CMat3<double> &, const CMat3<double> &);
template CMat3<float> operator-(const CMat3<float> &, const CMat3<float> &);
#endif

// ------------------------------

template<typename T>
std::ostream &operator<<(std::ostream &output, const CMat3<T> &m) {
  output.setf(std::ios::scientific);
  output << m.p_[0 * 3 + 0] << " " << m.p_[0 * 3 + 1] << " " << m.p_[0 * 3 + 2] << " ";
  output << m.p_[1 * 3 + 0] << " " << m.p_[1 * 3 + 1] << " " << m.p_[1 * 3 + 2] << " ";
  output << m.p_[2 * 3 + 0] << " " << m.p_[2 * 3 + 1] << " " << m.p_[2 * 3 + 2] << " ";
  return output;
}

template<typename T>
std::istream &operator>>(std::istream &input, CMat3<T> &m) {
  input >> m.p_[0 * 3 + 0] >> m.p_[0 * 3 + 1] >> m.p_[0 * 3 + 2];
  input >> m.p_[1 * 3 + 0] >> m.p_[1 * 3 + 1] >> m.p_[1 * 3 + 2];
  input >> m.p_[2 * 3 + 0] >> m.p_[2 * 3 + 1] >> m.p_[2 * 3 + 2];
  return input;
}

}

// -------------------------------------------------------------------

template<typename T>
delfem2::CMat3<T>::CMat3(): p_{0, 0, 0, 0, 0, 0, 0, 0, 0} {}
#ifdef DFM2_STATIC_LIBRARY
template delfem2::CMat3<float>::CMat3();
template delfem2::CMat3<double>::CMat3();
#endif

// ---------------------

template<typename T>
delfem2::CMat3<T>::CMat3(const T s): p_{s, 0, 0, 0, s, 0, 0, 0, s} {}
#ifdef DFM2_STATIC_LIBRARY
template delfem2::CMat3<float>::CMat3(float);
template delfem2::CMat3<double>::CMat3(double);
#endif

// ----------------------

template<typename T>
delfem2::CMat3<T>::CMat3
    (T v00, T v01, T v02,
     T v10, T v11, T v12,
     T v20, T v21, T v22):
    p_{v00, v01, v02, v10, v11, v12, v20, v21, v22} {}
#ifdef DFM2_STATIC_LIBRARY
template delfem2::CMat3<float>::CMat3(float v00, float v01, float v02,
                                      float v10, float v11, float v12,
                                      float v20, float v21, float v22);
template delfem2::CMat3<double>::CMat3(double v00, double v01, double v02,
                                       double v10, double v11, double v12,
                                       double v20, double v21, double v22);
#endif


// ----------------------

template<typename T>
delfem2::CMat3<T>::CMat3(T x, T y, T z):
    p_{x, 0, 0, 0, y, 0, 0, 0, z} {}
#ifdef DFM2_STATIC_LIBRARY
template delfem2::CMat3<float>::CMat3(float, float, float);
template delfem2::CMat3<double>::CMat3(double, double, double);
#endif

// ----------------------

template<typename T>
delfem2::CMat3<T>::CMat3(const T m[9]):
    p_{m[0], m[1], m[2], m[3], m[4], m[5], m[6], m[7], m[8]} {}
#ifdef DFM2_STATIC_LIBRARY
template delfem2::CMat3<double>::CMat3(const double m[9]);
template delfem2::CMat3<float>::CMat3(const float m[9]);
#endif

template<typename T>
void delfem2::CMat3<T>::MatVec(const T vec0[], T vec1[]) const {
  ::delfem2::MatVec3(vec1, p_, vec0);
}
#ifdef DFM2_STATIC_LIBRARY
template void delfem2::CMat3<float>::MatVec(const float vec0[], float vec1[]) const;
template void delfem2::CMat3<double>::MatVec(const double vec0[], double vec1[]) const;
#endif

// -------------------------------

template<typename T>
void delfem2::CMat3<T>::MatVecTrans(const T vec0[], T vec1[]) const {
  vec1[0] = p_[0] * vec0[0] + p_[3] * vec0[1] + p_[6] * vec0[2];
  vec1[1] = p_[1] * vec0[0] + p_[4] * vec0[1] + p_[7] * vec0[2];
  vec1[2] = p_[2] * vec0[0] + p_[5] * vec0[1] + p_[8] * vec0[2];
}
#ifdef DFM2_STATIC_LIBRARY
template void delfem2::CMat3f::MatVecTrans(const float vec0[], float vec1[]) const;
template void delfem2::CMat3d::MatVecTrans(const double vec0[], double vec1[]) const;
#endif

// --------------------------------

template<typename T>
delfem2::CMat3<T> delfem2::CMat3<T>::MatMat(const CMat3<T> &mat0) const {
  CMat3 m;
  ::delfem2::MatMat3(
      m.p_,
      this->p_, mat0.p_);
  return m;
}
#ifdef DFM2_STATIC_LIBRARY
template delfem2::CMat3<double> delfem2::CMat3<double>::MatMat(const CMat3<double> &mat0) const;
#endif

// --------------------------------

template<typename T>
delfem2::CMat3<T> delfem2::CMat3<T>::MatMatTrans(const CMat3<T> &mat0) const {
  CMat3 m;
  for (unsigned int i = 0; i < 3; i++) {
    for (unsigned int j = 0; j < 3; j++) {
      m.mat[i * 3 + j] =
          p_[0 * 3 + i] * mat0.p_[0 * 3 + j] +
          p_[1 * 3 + i] * mat0.p_[1 * 3 + j] +
          p_[2 * 3 + i] * mat0.p_[2 * 3 + j];
    }
  }
  return m;
}

template<typename T>
delfem2::CMat3<T> delfem2::CMat3<T>::Inverse() const {
  CMat3 mi = *this;
  mi.SetInverse();
  return mi;
}
#ifdef DFM2_STATIC_LIBRARY
template delfem2::CMat3<double> delfem2::CMat3<double>::Inverse() const;
#endif

// ------------------------------------------------------------------

template<typename T>
void delfem2::CMat3<T>::SetInverse() {
  ::delfem2::Inverse_Mat3(p_);
}
#ifdef DFM2_STATIC_LIBRARY
template void delfem2::CMat3<double>::SetInverse();
template void delfem2::CMat3<float>::SetInverse();
#endif

template<typename T>
void delfem2::CMat3<T>::SetSymetric(const T sm[6]) {
  p_[0] = sm[0];
  p_[1] = sm[5];
  p_[2] = sm[4];
  p_[3] = sm[5];
  p_[4] = sm[1];
  p_[5] = sm[3];
  p_[6] = sm[4];
  p_[7] = sm[3];
  p_[8] = sm[2];
}
#ifdef DFM2_STATIC_LIBRARY
template void delfem2::CMat3<float>::SetSymetric(const float sm[6]);
template void delfem2::CMat3<double>::SetSymetric(const double sm[6]);
#endif

// --------------------------------

template<typename T>
void delfem2::CMat3<T>::setZero() {
  for (auto &v : p_) { v = 0.0; }
}
#ifdef DFM2_STATIC_LIBRARY
template void delfem2::CMat3<float>::setZero();
template void delfem2::CMat3<double>::setZero();
#endif

// --------------------

namespace delfem2 {

template<>
DFM2_INLINE void CMat3<double>::SetRandom() {
  std::random_device rd;
  std::mt19937 mt(rd());
  std::uniform_real_distribution<double> dist(-50.0, 50.0);
  for (double &v : p_) { v = dist(mt); }
}

}

// -----------------------------------------

template<typename REAL>
void delfem2::CMat3<REAL>::SetRotMatrix_Cartesian(const REAL vec[]) {
  REAL sqt = vec[0] * vec[0] + vec[1] * vec[1] + vec[2] * vec[2];
  if (sqt < 1.0e-20) { // infinitesmal rotation approximation
    p_[0] = 1;
    p_[1] = -vec[2];
    p_[2] = +vec[1];
    p_[3] = +vec[2];
    p_[4] = 1;
    p_[5] = -vec[0];
    p_[6] = -vec[1];
    p_[7] = +vec[0];
    p_[8] = 1;
    return;
  }
  REAL t = std::sqrt(sqt);
  REAL invt = 1 / t;
  REAL n[3] = {vec[0] * invt, vec[1] * invt, vec[2] * invt};
  const REAL c0 = cos(t);
  const REAL s0 = sin(t);
  p_[0 * 3 + 0] = c0 + (1 - c0) * n[0] * n[0];
  p_[0 * 3 + 1] = -n[2] * s0 + (1 - c0) * n[0] * n[1];
  p_[0 * 3 + 2] = +n[1] * s0 + (1 - c0) * n[0] * n[2];
  p_[1 * 3 + 0] = +n[2] * s0 + (1 - c0) * n[1] * n[0];
  p_[1 * 3 + 1] = c0 + (1 - c0) * n[1] * n[1];
  p_[1 * 3 + 2] = -n[0] * s0 + (1 - c0) * n[1] * n[2];
  p_[2 * 3 + 0] = -n[1] * s0 + (1 - c0) * n[2] * n[0];
  p_[2 * 3 + 1] = +n[0] * s0 + (1 - c0) * n[2] * n[1];
  p_[2 * 3 + 2] = c0 + (1 - c0) * n[2] * n[2];
}
#ifdef DFM2_STATIC_LIBRARY
template void delfem2::CMat3<float>::SetRotMatrix_Cartesian(const float vec[]);
template void delfem2::CMat3<double>::SetRotMatrix_Cartesian(const double vec[]);
#endif

// -------------------------------

template<typename T>
void delfem2::CMat3<T>::SetRotMatrix_Cartesian(T x, T y, T z) {
  const T vec[3] = {x, y, z};
  this->SetRotMatrix_Cartesian(vec);
}
#ifdef DFM2_STATIC_LIBRARY
template void delfem2::CMat3<float>::SetRotMatrix_Cartesian(float x, float y, float z);
template void delfem2::CMat3<double>::SetRotMatrix_Cartesian(double x, double y, double z);
#endif

// ----------------------------------

template<typename T>
void delfem2::CMat3<T>::SetRotMatrix_Rodrigues(const T vec[]) {
  constexpr T half = static_cast<T>(0.5);
  constexpr T quarter = static_cast<T>(0.25);
  const T sqlen = vec[0] * vec[0] + vec[1] * vec[1] + vec[2] * vec[2];
  const T tmp1 = 1 / (1 + quarter * sqlen);
  p_[0] = 1 + tmp1 * (+half * vec[0] * vec[0] - half * sqlen);
  p_[1] = +tmp1 * (-vec[2] + half * vec[0] * vec[1]);
  p_[2] = +tmp1 * (+vec[1] + half * vec[0] * vec[2]);
  p_[3] = +tmp1 * (+vec[2] + half * vec[1] * vec[0]);
  p_[4] = 1 + tmp1 * (+half * vec[1] * vec[1] - half * sqlen);
  p_[5] = +tmp1 * (-vec[0] + half * vec[1] * vec[2]);
  p_[6] = +tmp1 * (-vec[1] + half * vec[2] * vec[0]);
  p_[7] = +tmp1 * (+vec[0] + half * vec[2] * vec[1]);
  p_[8] = 1 + tmp1 * (+half * vec[2] * vec[2] - half * sqlen);
}
#ifdef DFM2_STATIC_LIBRARY
template void delfem2::CMat3<float>::SetRotMatrix_Rodrigues(const float vec[]);
template void delfem2::CMat3<double>::SetRotMatrix_Rodrigues(const double vec[]);
#endif


// ----------------------------------

template<typename T>
void delfem2::CMat3<T>::SetRotMatrix_CRV(const T crv[]) {
  constexpr T one8th = static_cast<T>(1.0 / 8.0);
  const T c0 = one8th * (16 - crv[0] * crv[0] - crv[1] * crv[1] - crv[2] * crv[2]);
  const T tmp = 1 / ((4 - c0) * (4 - c0));
  p_[0 * 3 + 0] = tmp * ((c0 * c0 + 8 * c0 - 16) + 2 * crv[0] * crv[0]);
  p_[0 * 3 + 1] = tmp * (2 * crv[0] * crv[1] - 2 * c0 * crv[2]);
  p_[0 * 3 + 2] = tmp * (2 * crv[0] * crv[2] + 2 * c0 * crv[1]);
  p_[1 * 3 + 0] = tmp * (2 * crv[1] * crv[0] + 2 * c0 * crv[2]);
  p_[1 * 3 + 1] = tmp * ((c0 * c0 + 8 * c0 - 16) + 2 * crv[1] * crv[1]);
  p_[1 * 3 + 2] = tmp * (2 * crv[1] * crv[2] - 2 * c0 * crv[0]);
  p_[2 * 3 + 0] = tmp * (2 * crv[2] * crv[0] - 2 * c0 * crv[1]);
  p_[2 * 3 + 1] = tmp * (2 * crv[2] * crv[1] + 2 * c0 * crv[0]);
  p_[2 * 3 + 2] = tmp * ((c0 * c0 + 8 * c0 - 16) + 2 * crv[2] * crv[2]);
}
#ifdef DFM2_STATIC_LIBRARY
template void delfem2::CMat3<float>::SetRotMatrix_CRV(const float crv[]);
template void delfem2::CMat3<double>::SetRotMatrix_CRV(const double crv[]);
#endif

// ----------------------------------

template<typename T>
void delfem2::CMat3<T>::SetRotMatrix_Quaternion(const T quat[]) {
  mat3::SetMatrix3_Quaternion(p_, quat);
}
#ifdef DFM2_STATIC_LIBRARY
template void delfem2::CMat3<float>::SetRotMatrix_Quaternion(const float quat[]);
template void delfem2::CMat3<double>::SetRotMatrix_Quaternion(const double quat[]);
#endif

// ----------------------------------

template<typename T>
void delfem2::CMat3<T>::SetRotMatrix_BryantAngle(T rx, T ry, T rz) {
  CMat3 mx;
  T rvx[3] = {rx, 0, 0};
  mx.SetRotMatrix_Cartesian(rvx);
  CMat3 my;
  T rvy[3] = {0, ry, 0};
  my.SetRotMatrix_Cartesian(rvy);
  CMat3 mz;
  T rvz[3] = {0, 0, rz};
  mz.SetRotMatrix_Cartesian(rvz);
  CMat3 m = mz;
  m = m.MatMat(my);
  m = m.MatMat(mx);
  *this = m;
}
#ifdef DFM2_STATIC_LIBRARY
template void delfem2::CMat3f::SetRotMatrix_BryantAngle(float rx, float ry, float rz);
template void delfem2::CMat3d::SetRotMatrix_BryantAngle(double rx, double ry, double rz);
#endif

// ---------------------------------

template<typename T>
void delfem2::CMat3<T>::GetQuat_RotMatrix(
    T quat[]) const {
  constexpr T one4th = static_cast<T>(0.25);
  const T smat[16] = {
      1 + p_[0 * 3 + 0] - p_[1 * 3 + 1] - p_[2 * 3 + 2],   // 00
      p_[0 * 3 + 1] + p_[1 * 3 + 0],  // 01
      p_[0 * 3 + 2] + p_[2 * 3 + 0],  // 02
      p_[2 * 3 + 1] - p_[1 * 3 + 2],  // 03
      p_[1 * 3 + 0] + p_[0 * 3 + 1],  // 10
      1 - p_[0 * 3 + 0] + p_[1 * 3 + 1] - p_[2 * 3 + 2],  // 11
      p_[1 * 3 + 2] + p_[2 * 3 + 1],  // 12
      p_[0 * 3 + 2] - p_[2 * 3 + 0],  // 13
      p_[0 * 3 + 2] + p_[2 * 3 + 0],  // 20
      p_[1 * 3 + 2] + p_[2 * 3 + 1],  // 21
      1 - p_[0 * 3 + 0] - p_[1 * 3 + 1] + p_[2 * 3 + 2],  // 22
      p_[1 * 3 + 0] - p_[0 * 3 + 1],  // 23
      p_[2 * 3 + 1] - p_[1 * 3 + 2],  // 30
      p_[0 * 3 + 2] - p_[2 * 3 + 0],  // 31
      p_[1 * 3 + 0] - p_[0 * 3 + 1],  // 32
      1 + p_[0 * 3 + 0] + p_[1 * 3 + 1] + p_[2 * 3 + 2],  // 33
  };

  unsigned int imax;
  imax = (smat[0 * 4 + 0] > smat[1 * 4 + 1]) ? 0 : 1;
  imax = (smat[imax * 4 + imax] > smat[2 * 4 + 2]) ? imax : 2;
  imax = (smat[imax * 4 + imax] > smat[3 * 4 + 3]) ? imax : 3;

  quat[imax] = std::sqrt(smat[imax * 4 + imax]) / 2;
  for (unsigned int k = 0; k < 4; k++) {
    if (k == imax) continue;
    quat[k] = smat[imax * 4 + k] * one4th / quat[imax];
  }
}
#ifdef DFM2_STATIC_LIBRARY
template void delfem2::CMat3<float>::GetQuat_RotMatrix(float quat[]) const;
template void delfem2::CMat3<double>::GetQuat_RotMatrix(double quat[]) const;
#endif

// -------------------------------

template<typename T>
void delfem2::CMat3<T>::SetIdentity(T scale) {
  p_[0] = scale;
  p_[1] = 0;
  p_[2] = 0;
  p_[3] = 0;
  p_[4] = scale;
  p_[5] = 0;
  p_[6] = 0;
  p_[7] = 0;
  p_[8] = scale;
}
#ifdef DFM2_STATIC_LIBRARY
template void delfem2::CMat3<double>::SetIdentity(double scale);
template void delfem2::CMat3<float>::SetIdentity(float scale);
#endif

