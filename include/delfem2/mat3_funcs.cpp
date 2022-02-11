/*
 * Copyright (c) 2019 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */

#include "delfem2/mat3_funcs.h"

#include <random>

#include "delfem2/vec3_funcs.h"

namespace delfem2::mat3 {

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

// --------------

template<typename T>
DFM2_INLINE void delfem2::Vec2_Mat3Vec2_Homography(
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
template void delfem2::Vec2_Mat3Vec2_Homography(float [2], const float [9], const float [2]);
template void delfem2::Vec2_Mat3Vec2_Homography(double [2], const double [9], const double [2]);
#endif
// --------------

template<typename T>
DFM2_INLINE std::array<T,2> delfem2::Vec2_Mat3Vec2_Homography(
    const T A[9],
    const T x[2]) {
  T y0 = A[0] * x[0] + A[1] * x[1] + A[2];
  T y1 = A[3] * x[0] + A[4] * x[1] + A[5];
  const T w = A[6] * x[0] + A[7] * x[1] + A[8];
  return {y0 / w, y1 / w};
}
#ifdef DFM2_STATIC_LIBRARY
template std::array<float,2> delfem2::Vec2_Mat3Vec2_Homography(const float [9], const float [2]);
template std::array<double,2> delfem2::Vec2_Mat3Vec2_Homography(const double [9], const double [2]);
#endif

// -------------

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
// --------------------------
// below: mat3 and quaternion

template <typename T>
DFM2_INLINE void delfem2::Quat_Mat3(
  T quat[4],
  const T p_[9]){
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
template void delfem2::Quat_Mat3(float quat[4], const float p_[9]);
template void delfem2::Quat_Mat3(double quat[4], const double p_[9]);
#endif

// ----------------

// https://en.wikipedia.org/wiki/Quaternions_and_spatial_rotation
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
T delfem2::SquareNormFrobenius_SymMat3(
    const T sm[6]) {
  return sm[0] * sm[0] + sm[1] * sm[1] + sm[2] * sm[2] + 2 * (sm[3] * sm[3] + sm[4] * sm[4] + sm[5] * sm[5]);
}
#ifdef DFM2_STATIC_LIBRARY
template float delfem2::SquareNormFrobenius_SymMat3(const float sm[6]);
template double delfem2::SquareNormFrobenius_SymMat3(const double sm[6]);
#endif

// ---------------------------------------

template<typename REAL>
void delfem2::Mat3_RotMatFromAxisAngleVec(
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
template void delfem2::Mat3_RotMatFromAxisAngleVec(float mat[9], const float vec[3]);
template void delfem2::Mat3_RotMatFromAxisAngleVec(double mat[9], const double vec[3]);
#endif

// ---------------------------------------

template <typename VEC, typename REAL>
std::array<REAL,9> delfem2::Mat3_RotMatFromAxisAngleVec(const VEC &vec) {
  REAL sqt = vec[0] * vec[0] + vec[1] * vec[1] + vec[2] * vec[2];
  if (sqt < 1.0e-20) { // infinitesmal rotation approximation
    return {
        1,
        -vec[2],
        +vec[1],
        +vec[2],
        1,
        -vec[0],
        -vec[1],
        +vec[0],
        1 };
  }
  REAL t = std::sqrt(sqt);
  REAL invt = 1 / t;
  REAL n[3] = {vec[0] * invt, vec[1] * invt, vec[2] * invt};
  const REAL c0 = cos(t);
  const REAL s0 = sin(t);
  return {
      c0 + (1 - c0) * n[0] * n[0],
      -n[2] * s0 + (1 - c0) * n[0] * n[1],
      +n[1] * s0 + (1 - c0) * n[0] * n[2],
      +n[2] * s0 + (1 - c0) * n[1] * n[0],
      c0 + (1 - c0) * n[1] * n[1],
      -n[0] * s0 + (1 - c0) * n[1] * n[2],
      -n[1] * s0 + (1 - c0) * n[2] * n[0],
      +n[0] * s0 + (1 - c0) * n[2] * n[1],
      c0 + (1 - c0) * n[2] * n[2] };
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
  Quat_Mat3(eparam2, mat);
//  const CMat3<T> m(mat);
//  m.GetQuaternion(eparam2);
  crv[0] = 4 * eparam2[0] / (1 + eparam2[3]);
  crv[1] = 4 * eparam2[1] / (1 + eparam2[3]);
  crv[2] = 4 * eparam2[2] / (1 + eparam2[3]);
}

template<typename T>
void delfem2::EulerAngle_Mat3(
  T ea[3],
  const T m[9],
  const std::array<int, 3> &axis_idxs) {
  if (axis_idxs[0] == 2 && axis_idxs[1] == 1 && axis_idxs[2] == 0) {
    const T y1 = -std::asin(m[2*3+0]);
    const T inv_cos_y1 = 1 / std::cos(y1);
    assert(std::fabs(cos(y1))>1.0e-10);
    ea[1] = y1;
    ea[0] = std::atan2(m[1*3+0] * inv_cos_y1, m[0*3+0] * inv_cos_y1);
    ea[2] = std::atan2(m[2*3+1] * inv_cos_y1, m[2*3+2] * inv_cos_y1);
  }
  else if (axis_idxs[0] == 2 && axis_idxs[1] == 0 && axis_idxs[2] == 1) {
    const T y1 = std::asin(m[2*3+1]);
    assert(std::fabs(cos(y1))>1.0e-10);
    const T inv_cos_y1 = 1 / std::cos(y1);
    ea[1] = y1;
    ea[0] = std::atan2(-m[0*3+1] * inv_cos_y1, m[1*3+1] * inv_cos_y1);
    ea[2] = std::atan2(-m[2*3+0] * inv_cos_y1, m[2*3+2] * inv_cos_y1);
  }
  else{
    assert(0);
  }
}
#ifdef DFM2_STATIC_LIBRARY
template void delfem2::EulerAngle_Mat3(
  double ea[3],
  const double m[9],
  const std::array<int, 3> &axis_idxs);
template void delfem2::EulerAngle_Mat3(
  float ea[3],
  const float m[9],
  const std::array<int, 3> &axis_idxs);
#endif

// ---------------

template<typename VEC, typename REAL>
std::array<REAL, 9> delfem2::Mat3_MinimumRotation(
    const VEC &V,
    const VEC &v) {
  VEC ep = V.normalized();
  VEC eq = v.normalized();
  VEC n = ep.cross(eq);
  const REAL st2 = n.dot(n);
  const REAL ct = ep.dot(eq);
  constexpr REAL half = static_cast<REAL>(0.5);

  if (st2 < 1.0e-8f) { // very small angle or n is zero
    // inifinitesimal rotation
    if (ct > 0.99) {
      return {
          1 + half * (n.x * n.x - st2),
          -n.z + half * (n.x * n.y),
          +n.y + half * (n.x * n.z),
          +n.z + half * (n.y * n.x),
          1 + half * (n.y * n.y - st2),
          -n.x + half * (n.y * n.z),
          -n.y + half * (n.z * n.x),
          +n.x + half * (n.z * n.y),
          1 + half * (n.z * n.z - st2)};
    } else {
      VEC epx, epy;
      FrameFromVectorZ(epx, epy, ep);
      const VEC eqx = epx - eq.dot(epx) * eq; // vector orthogonal to eq
      const VEC eqy = eq.cross(eqx);
      return {
          eqx.dot(epx), eqy.dot(epx), eq.dot(epx),
          eqx.dot(epy), eqy.dot(epy), eq.dot(epy),
          eqx.dot(ep), eqy.dot(ep), eq.dot(ep)};
      /*
      const CMat3<REAL> Rp = Mat3_3Bases(epx, epy, ep);
      const CMat3<REAL> Rq = Mat3_3Bases(eqx, eqy, eq);
      return Rq * Rp.transpose();
       */
    }
  }
  const REAL st = std::sqrt(st2);
  n.normalize();
  // Rodoriguez's rotation formula
  return {
      ct + (1 - ct) * n.x * n.x,
      -n.z * st + (1 - ct) * n.x * n.y,
      +n.y * st + (1 - ct) * n.x * n.z,
      +n.z * st + (1 - ct) * n.y * n.x,
      ct + (1 - ct) * n.y * n.y,
      -n.x * st + (1 - ct) * n.y * n.z,
      -n.y * st + (1 - ct) * n.z * n.x,
      +n.x * st + (1 - ct) * n.z * n.y,
      ct + (1 - ct) * n.z * n.z};
}

// ============================================

#ifdef DFM2_STATIC_LIBRARY

#include "delfem2/vec3.h"
//
namespace delfem2 {
using f0 = float [3];
using d0 = double [3];
using f1 = float *;
using d1 = double *;
using f2 = std::array<float, 3>;
using d2 = std::array<double, 3>;
using f3 = CVec3f;
using d3 = CVec3d;
//
template std::array<float, 9> Mat3_MinimumRotation(const f3 &, const f3 &);
template std::array<double, 9> Mat3_MinimumRotation(const d3 &, const d3 &);
//
template std::array<double,9> Mat3_RotMatFromAxisAngleVec(const d0 &vec);
template std::array<double,9> Mat3_RotMatFromAxisAngleVec(const d1 &vec);
template std::array<double,9> Mat3_RotMatFromAxisAngleVec(const d2 &vec);
template std::array<double,9> Mat3_RotMatFromAxisAngleVec(const d3 &vec);
template std::array<float,9> Mat3_RotMatFromAxisAngleVec(const f0 &vec);
template std::array<float,9> Mat3_RotMatFromAxisAngleVec(const f1 &vec);
template std::array<float,9> Mat3_RotMatFromAxisAngleVec(const f2 &vec);
template std::array<float,9> Mat3_RotMatFromAxisAngleVec(const f3 &vec);
}
#endif
