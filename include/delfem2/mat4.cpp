/*
 * Copyright (c) 2019 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */

#include "delfem2/mat4.h"

#include <cstring>
#include <climits>

// ------------------------

namespace delfem2::mat4 {

template<typename REAL>
DFM2_INLINE void CalcInvMat(
    REAL *a,
    const unsigned int n,
    int &info) {
  REAL tmp1;

  info = 0;
  for (unsigned int i = 0; i < n; i++) {
    if (fabs(a[i * n + i]) < 1.0e-30) {
      info = 1;
      return;
    }
    if (a[i * n + i] < 0.0) {
      info--;
    }
    tmp1 = 1 / a[i * n + i];
    a[i * n + i] = 1.0;
    for (unsigned int k = 0; k < n; k++) {
      a[i * n + k] *= tmp1;
    }
    for (unsigned int j = 0; j < n; j++) {
      if (j != i) {
        tmp1 = a[j * n + i];
        a[j * n + i] = 0.0;
        for (unsigned int k = 0; k < n; k++) {
          a[j * n + k] -= tmp1 * a[i * n + k];
        }
      }
    }
  }
}

template<typename REAL>
bool CalcInvMatPivot(REAL *a, unsigned int n, unsigned int *tmp) {
  unsigned int *row = tmp;
  for (unsigned int ipv = 0; ipv < n; ipv++) {
    // find maximum
    REAL big = 0.0;
    unsigned int pivot_row = UINT_MAX;
    for (unsigned int i = ipv; i < n; i++) {
      if (fabs(a[i * n + ipv]) < big) { continue; }
      big = fabs(a[i * n + ipv]);
      pivot_row = i;
    }
    if (pivot_row == UINT_MAX) { return false; }
    row[ipv] = pivot_row;

    // swapping column
    if (ipv != pivot_row) {
      for (unsigned int i = 0; i < n; i++) {
        REAL temp = a[ipv * n + i];
        a[ipv * n + i] = a[pivot_row * n + i];
        a[pivot_row * n + i] = temp;
      }
    }

    // set diagonal 1 for for pivotting column
    REAL inv_pivot = 1 / a[ipv * n + ipv];
    a[ipv * n + ipv] = 1.0;
    for (unsigned int j = 0; j < n; j++) {
      a[ipv * n + j] *= inv_pivot;
    }

    // set pivot column 0 except for pivotting column
    for (unsigned int i = 0; i < n; i++) {
      if (i == ipv) { continue; }
      REAL temp = a[i * n + ipv];
      a[i * n + ipv] = 0.0;
      for (unsigned int j = 0; j < n; j++) {
        a[i * n + j] -= temp * a[ipv * n + j];
      }
    }

  }

  // swaping column
  for (int j = int(n - 1); j >= 0; j--) {
    if ((unsigned int) j == row[j]) { continue; }
    for (unsigned int i = 0; i < n; i++) {
      REAL temp = a[i * n + j];
      a[i * n + j] = a[i * n + row[j]];
      a[i * n + row[j]] = temp;
    }
  }
  return true;
}

template<typename REAL>
DFM2_INLINE void Normalize3D(
    REAL vec[3]) {
  const REAL len = std::sqrt(vec[0] * vec[0] + vec[1] * vec[1] + vec[2] * vec[2]);
  const REAL leninv = 1 / len;
  vec[0] *= leninv;
  vec[1] *= leninv;
  vec[2] *= leninv;
}

template<typename REAL>
DFM2_INLINE void Cross3D(
    REAL r[3], const REAL v1[3], const REAL v2[3]) {
  r[0] = v1[1] * v2[2] - v2[1] * v1[2];
  r[1] = v1[2] * v2[0] - v2[2] * v1[0];
  r[2] = v1[0] * v2[1] - v2[0] * v1[1];
}

}


// ------------------------------
// below: mat4



template<typename REAL>
DFM2_INLINE void delfem2::Mat4_AffineTransProjectionOrtho(
    REAL mP[16],
    double xmin, double xmax, // -x, +x
    double ymin, double ymax, // -y, +y
    double zmin, double zmax) // -z, +z
{
  // column 0
  mP[0 * 4 + 0] = static_cast<REAL>(2.0 / (xmax - xmin));
  mP[0 * 4 + 1] = 0;
  mP[0 * 4 + 2] = 0;
  mP[0 * 4 + 3] = 0;
  // column 1
  mP[1 * 4 + 0] = 0;
  mP[1 * 4 + 1] = static_cast<REAL>(2.0 / (ymax - ymin));
  mP[1 * 4 + 2] = 0;
  mP[1 * 4 + 3] = 0;
  // column 2
  mP[2 * 4 + 0] = 0;
  mP[2 * 4 + 1] = 0;
  mP[2 * 4 + 2] = static_cast<REAL>(2.0 / (zmax - zmin));
  mP[2 * 4 + 3] = 0;
  // collumn 3
  mP[3 * 4 + 0] = static_cast<REAL>(-(xmin + xmax) / (xmax - xmin));
  mP[3 * 4 + 1] = static_cast<REAL>(-(ymax + ymin) / (ymax - ymin));
  mP[3 * 4 + 2] = static_cast<REAL>(-(zmax + zmin) / (zmax - zmin));
  mP[3 * 4 + 3] = 1;
}
#ifdef DFM2_STATIC_LIBRARY
template void delfem2::Mat4_AffineTransProjectionOrtho(
    double mP[16],
    double xmin, double xmax, // -x, +x
    double ymin, double ymax, // -y, +y
    double zmin, double zmax); // -z, +z
template void delfem2::Mat4_AffineTransProjectionOrtho(
    float mP[16],
    double xmin, double xmax,  // -x, +x
    double ymin, double ymax,  // -y, +y
    double zmin, double zmax);  // -z, +z
#endif

// -----------------------

template<typename REAL>
DFM2_INLINE void delfem2::Mat4_AffineTransProjectionFrustum(
    REAL mP[16],
    REAL fovyInRad,
    REAL aspectRatio,
    REAL zmin,
    REAL zmax) {
  const REAL yratio = 1 / std::tan(fovyInRad / 2); // how z change w.r.t. the y change
  const REAL xratio = yratio / aspectRatio;
  // column 0
  mP[0 * 4 + 0] = xratio;
  mP[0 * 4 + 1] = 0.0;
  mP[0 * 4 + 2] = 0.0;
  mP[0 * 4 + 3] = 0.0;
  // column 1
  mP[1 * 4 + 0] = 0.0;
  mP[1 * 4 + 1] = yratio;
  mP[1 * 4 + 2] = 0.0;
  mP[1 * 4 + 3] = 0.0;
  // column 2
  mP[2 * 4 + 0] = 0.0;
  mP[2 * 4 + 1] = 0.0;
  mP[2 * 4 + 2] = -(zmin + zmax) / (zmax - zmin);
  mP[2 * 4 + 3] = -1.0;
  // column 3
  mP[3 * 4 + 0] = 0.0;
  mP[3 * 4 + 1] = 0.0;
  mP[3 * 4 + 2] = +(zmin * zmax * 2) / (zmax - zmin);
  mP[3 * 4 + 3] = 0.0;
}
#ifdef DFM2_STATIC_LIBRARY
template DFM2_INLINE void delfem2::Mat4_AffineTransProjectionFrustum(
    float mP[16],
    float fovyInRad,
    float aspectRatio,
    float zmin,
    float zmax);
template DFM2_INLINE void delfem2::Mat4_AffineTransProjectionFrustum(
    double mP[16],
    double fovyInRad,
    double aspectRatio,
    double zmin,
    double zmax);
#endif


// -----------------------------------------




template<typename T0, typename T1>
DFM2_INLINE void delfem2::MultMat4AffineTransTranslateFromRight(
    T0 *matrix,
    T1 x,
    T1 y,
    T1 z) {
  matrix[12] = matrix[0] * x + matrix[4] * y + matrix[8] * z + matrix[12];
  matrix[13] = matrix[1] * x + matrix[5] * y + matrix[9] * z + matrix[13];
  matrix[14] = matrix[2] * x + matrix[6] * y + matrix[10] * z + matrix[14];
  matrix[15] = matrix[3] * x + matrix[7] * y + matrix[11] * z + matrix[15];
}

template<typename REAL>
DFM2_INLINE void delfem2::Mat4_AffineTransLookAt(
    REAL *Mr,
    REAL eyex, REAL eyey, REAL eyez,
    REAL cntx, REAL cnty, REAL cntz,
    REAL upx, REAL upy, REAL upz) {
  const REAL eyePosition3D[3] = {eyex, eyey, eyez};
  const REAL center3D[3] = {cntx, cnty, cntz};
  const REAL upVector3D[3] = {upx, upy, upz};
  // ------------------
  REAL forward[3] = {
      center3D[0] - eyePosition3D[0],
      center3D[1] - eyePosition3D[1],
      center3D[2] - eyePosition3D[2]};
  mat4::Normalize3D(forward);
  // ------------------
  // Side = forward x up
  REAL side[3] = {1, 0, 0};
  mat4::Cross3D(side, forward, upVector3D);
  mat4::Normalize3D(side);
  // ------------------
  //Recompute up as: up = side x forward
  REAL up[3] = {0, 1, 0};
  mat4::Cross3D(up, side, forward);
  // ------------------
  Mr[0] = side[0];
  Mr[4] = side[1];
  Mr[8] = side[2];
  Mr[12] = 0.0;
  // ------------------
  Mr[1] = up[0];
  Mr[5] = up[1];
  Mr[9] = up[2];
  Mr[13] = 0.0;
  // ------------------
  Mr[2] = -forward[0];
  Mr[6] = -forward[1];
  Mr[10] = -forward[2];
  Mr[14] = 0.0;
  // ------------------
  Mr[3] = 0.0;
  Mr[7] = 0.0;
  Mr[11] = 0.0;
  Mr[15] = 1.0;
  // ------------------
  delfem2::MultMat4AffineTransTranslateFromRight(
      Mr,
      -eyePosition3D[0], -eyePosition3D[1], -eyePosition3D[2]);
}
#ifdef DFM2_STATIC_LIBRARY
template DFM2_INLINE void delfem2::Mat4_AffineTransLookAt(
    double *Mr,
    double eyex, double eyey, double eyez,
    double cntx, double cnty, double cntz,
    double upx, double upy, double upz);
template DFM2_INLINE void delfem2::Mat4_AffineTransLookAt(
    float *Mr,
    float eyex, float eyey, float eyez,
    float cntx, float cnty, float cntz,
    float upx, float upy, float upz);
#endif


// ------------------------
// below: mat vec

template<typename T>
DFM2_INLINE void delfem2::Mat4Vec3(
    T vo[3],
    const T M[16],
    const T vi[3]) {
  vo[0] = M[0 * 4 + 0] * vi[0] + M[0 * 4 + 1] * vi[1] + M[0 * 4 + 2] * vi[2];
  vo[1] = M[1 * 4 + 0] * vi[0] + M[1 * 4 + 1] * vi[1] + M[1 * 4 + 2] * vi[2];
  vo[2] = M[2 * 4 + 0] * vi[0] + M[2 * 4 + 1] * vi[1] + M[2 * 4 + 2] * vi[2];
}
#ifdef DFM2_STATIC_LIBRARY
template void delfem2::Mat4Vec3(float vo[3], const float M[16], const float vi[3]);
template void delfem2::Mat4Vec3(double vo[3], const double M[16], const double vi[3]);
#endif

DFM2_INLINE void delfem2::Vec3Mat4(
    double vo[3],
    const double vi[3],
    const double M[16]) {
  vo[0] = vi[0] * M[0 * 4 + 0] + vi[1] * M[1 * 4 + 0] + vi[2] * M[2 * 4 + 0];
  vo[1] = vi[0] * M[0 * 4 + 1] + vi[1] * M[1 * 4 + 1] + vi[2] * M[2 * 4 + 1];
  vo[2] = vi[0] * M[0 * 4 + 2] + vi[1] * M[1 * 4 + 2] + vi[2] * M[2 * 4 + 2];
}

template<typename T>
DFM2_INLINE void delfem2::MatVec4(
    T v[4],
    const T A[16],
    const T x[4]) {
  v[0] = A[0 * 4 + 0] * x[0] + A[0 * 4 + 1] * x[1] + A[0 * 4 + 2] * x[2] + A[0 * 4 + 3] * x[3];
  v[1] = A[1 * 4 + 0] * x[0] + A[1 * 4 + 1] * x[1] + A[1 * 4 + 2] * x[2] + A[1 * 4 + 3] * x[3];
  v[2] = A[2 * 4 + 0] * x[0] + A[2 * 4 + 1] * x[1] + A[2 * 4 + 2] * x[2] + A[2 * 4 + 3] * x[3];
  v[3] = A[3 * 4 + 0] * x[0] + A[3 * 4 + 1] * x[1] + A[3 * 4 + 2] * x[2] + A[3 * 4 + 3] * x[3];
}
#ifdef DFM2_STATIC_LIBRARY
template void delfem2::MatVec4(
    float v[4], const float A[16], const float x[4]);
template void delfem2::MatVec4(
    double v[4], const double A[16], const double x[4]);
#endif

template<typename T>
DFM2_INLINE void delfem2::VecMat4(
    T v[4],
    const T x[4],
    const T A[16]) {
  v[0] = A[0 * 4 + 0] * x[0] + A[1 * 4 + 0] * x[1] + A[2 * 4 + 0] * x[2] + A[3 * 4 + 0] * x[3];
  v[1] = A[0 * 4 + 1] * x[0] + A[1 * 4 + 1] * x[1] + A[2 * 4 + 1] * x[2] + A[3 * 4 + 1] * x[3];
  v[2] = A[0 * 4 + 2] * x[0] + A[1 * 4 + 2] * x[1] + A[2 * 4 + 2] * x[2] + A[3 * 4 + 2] * x[3];
  v[3] = A[0 * 4 + 3] * x[0] + A[1 * 4 + 3] * x[1] + A[2 * 4 + 3] * x[2] + A[3 * 4 + 3] * x[3];
}
#ifdef DFM2_STATIC_LIBRARY
template void delfem2::VecMat4(
    float v[4], const float x[4], const float A[16]);
template void delfem2::VecMat4(
    double v[4], const double x[4], const double A[16]);
#endif


// ---------------------------

template<typename T0, typename T1, typename T2>
DFM2_INLINE void delfem2::Vec3_Mat4Vec3_AffineProjection(
    T0 y0[3],
    const T1 a[16],
    const T2 x0[3]) {
  const T1 x1[4] = {(T1) x0[0], (T1) x0[1], (T1) x0[2], (T1) 1};
  T1 y1[4];
  MatVec4(y1, a, x1);
  y0[0] = y1[0] / y1[3];
  y0[1] = y1[1] / y1[3];
  y0[2] = y1[2] / y1[3];
}
#ifdef DFM2_STATIC_LIBRARY
template void delfem2::Vec3_Mat4Vec3_AffineProjection(
    float y0[3], const float a[16], const float x0[3]);
template void delfem2::Vec3_Mat4Vec3_AffineProjection(
    double y0[3], const double a[16], const double x0[3]);
#endif

// ---------------------------

template<typename T>
DFM2_INLINE
std::array<T,2> delfem2::Vec2_Mat4Vec3_AffineProjection(
    const T a[16],
    const T x0[3]) {
  const T x1[4] = { x0[0], x0[1], x0[2], 1};
  T y1[4];
  MatVec4(y1, a, x1);
  return { y1[0] / y1[3], y1[1] / y1[3] };
}
#ifdef DFM2_STATIC_LIBRARY
template std::array<float,2> delfem2::Vec2_Mat4Vec3_AffineProjection(
    const float a[16],
    const float x0[3]);
template std::array<double,2> delfem2::Vec2_Mat4Vec3_AffineProjection(
    const double a[16],
    const double x0[3]);
#endif

// ----------------------------

template<typename T0, typename T1, typename T2>
DFM2_INLINE void delfem2::Vec3_Vec3Mat4_AffineProjection(
    T0 y0[3],
    const T1 x0[3],
    const T2 a[16]) {
  const T2 x1[4] = {(T2) x0[0], (T2) x0[1], (T2) x0[2], 1};
  T2 y1[4];
  VecMat4(y1, x1, a);
  y0[0] = y1[0] / y1[3];
  y0[1] = y1[1] / y1[3];
  y0[2] = y1[2] / y1[3];
}
#ifdef DFM2_STATIC_LIBRARY
template void delfem2::Vec3_Vec3Mat4_AffineProjection(
    float y0[3], const float x0[3], const float a[16]);
template void delfem2::Vec3_Vec3Mat4_AffineProjection(
    double y0[3], const double x0[3], const double a[16]);
#endif

// ----------------------

template<typename T>
DFM2_INLINE void delfem2::Vec3_Mat4Vec3_Affine(
    T y0[3],
    const T a[16],
    const T x0[3]) {
  const T x1[4] = {x0[0], x0[1], x0[2], 1.0};
  T y1[4];
  MatVec4(y1, a, x1);
  y0[0] = y1[0];
  y0[1] = y1[1];
  y0[2] = y1[2];
}
#ifdef DFM2_STATIC_LIBRARY
template void delfem2::Vec3_Mat4Vec3_Affine(
    float y0[3], const float a[16], const float x0[3]);
template void delfem2::Vec3_Mat4Vec3_Affine(
    double y0[3], const double a[16], const double x0[3]);
#endif


// ----------------------

template<typename T>
DFM2_INLINE void delfem2::Mat4_AffineScale
    (T A[16],
     T s) {
  for (int i = 0; i < 16; ++i) { A[i] = 0.0; }
  A[0 * 4 + 0] = s;
  A[1 * 4 + 1] = s;
  A[2 * 4 + 2] = s;
  A[3 * 4 + 3] = 1.0;
}
#ifdef DFM2_STATIC_LIBRARY
template void delfem2::Mat4_AffineScale(float A[16], float s);
template void delfem2::Mat4_AffineScale(double A[16], double s);
#endif

// ------------------------

template<typename T>
DFM2_INLINE void delfem2::Mat4_AffineTranslation
    (T A[16],
     T dx, T dy, T dz) {
  for (auto i = 0; i < 16; ++i) { A[i] = 0.0; }
  for (int i = 0; i < 4; ++i) { A[i * 4 + i] = 1.0; }
  A[0 * 4 + 3] = dx;
  A[1 * 4 + 3] = dy;
  A[2 * 4 + 3] = dz;
}
#ifdef DFM2_STATIC_LIBRARY
template void delfem2::Mat4_AffineTranslation(
    float A[16],
    float dx, float dy, float dz);
template void delfem2::Mat4_AffineTranslation(
    double A[16],
    double dx, double dy, double dz);
#endif

template<typename T>
DFM2_INLINE void
delfem2::Mat4_AffineTranslation(
    T A[16],
    const T v[3]) {
  A[0 * 4 + 0] = 1;
  A[0 * 4 + 1] = 0;
  A[0 * 4 + 2] = 0;
  A[0 * 4 + 3] = v[0];
  A[1 * 4 + 0] = 0;
  A[1 * 4 + 1] = 1;
  A[1 * 4 + 2] = 0;
  A[1 * 4 + 3] = v[1];
  A[2 * 4 + 0] = 0;
  A[2 * 4 + 1] = 0;
  A[2 * 4 + 2] = 1;
  A[2 * 4 + 3] = v[2];
  A[3 * 4 + 0] = 0;
  A[3 * 4 + 1] = 0;
  A[3 * 4 + 2] = 0;
  A[3 * 4 + 3] = 1;
}
#ifdef DFM2_STATIC_LIBRARY
template void delfem2::Mat4_AffineTranslation(
    float A[16],
    const float v[3]);
template void delfem2::Mat4_AffineTranslation(
    double A[16],
    const double v[3]);
#endif

template<typename REAL>
DFM2_INLINE void delfem2::Mat4_AffineTransTranslate(
    REAL r[],
    const REAL t[]) {
  // column 0
  r[0] = 1;
  r[1] = 0;
  r[2] = 0;
  r[3] = 0;
  // column 1
  r[4] = 0;
  r[5] = 1;
  r[6] = 0;
  r[7] = 0;
  // column 2
  r[8] = 0;
  r[9] = 0;
  r[10] = 1;
  r[11] = 0;
  // column 3
  r[12] = t[0];
  r[13] = t[1];
  r[14] = t[2];
  r[15] = 1;
}
#ifdef DFM2_STATIC_LIBRARY
template void delfem2::Mat4_AffineTransTranslate(float r[], const float t[]);
template void delfem2::Mat4_AffineTransTranslate(double r[], const double t[]);
#endif


// --------------------------

template<typename T>
DFM2_INLINE void delfem2::Mat4_AffineRotationRodriguez(
    T A[16],
    T dx, T dy, T dz) {
  constexpr T half = static_cast<T>(0.5);
  constexpr T one4th = static_cast<T>(0.25);
  for (int i = 0; i < 16; ++i) { A[i] = 0; }
  //
  const T sqlen = dx * dx + dy * dy + dz * dz;
  const T tmp1 = 1 / (1 + one4th * sqlen);
  A[0 * 4 + 0] = 1 + tmp1 * (+half * dx * dx - half * sqlen);
  A[0 * 4 + 1] = +tmp1 * (-dz + half * dx * dy);
  A[0 * 4 + 2] = +tmp1 * (+dy + half * dx * dz);
  A[0 * 4 + 3] = 0;
  //
  A[1 * 4 + 0] = +tmp1 * (+dz + half * dy * dx);
  A[1 * 4 + 1] = 1 + tmp1 * (+half * dy * dy - half * sqlen);
  A[1 * 4 + 2] = +tmp1 * (-dx + half * dy * dz);
  A[1 * 4 + 3] = 0;
  //
  A[2 * 4 + 0] = +tmp1 * (-dy + half * dz * dx);
  A[2 * 4 + 1] = +tmp1 * (+dx + half * dz * dy);
  A[2 * 4 + 2] = 1 + tmp1 * (+half * dz * dz - half * sqlen);
  A[2 * 4 + 3] = 0;
  //
  A[3 * 4 + 0] = 0;
  A[3 * 4 + 1] = 0;
  A[3 * 4 + 2] = 0;
  A[3 * 4 + 3] = 1;
}
#ifdef DFM2_STATIC_LIBRARY
template void delfem2::Mat4_AffineRotationRodriguez(
    float A[16],
    float dx, float dy, float dz);
template void delfem2::Mat4_AffineRotationRodriguez(
    double A[16],
    double dx, double dy, double dz);
#endif

// ------------------------------------------------

template<typename REAL>
void delfem2::Mat4_Identity(
    REAL A[16]) {
  for (int i = 0; i < 16; ++i) { A[i] = 0; }
  A[0 * 4 + 0] = 1;
  A[1 * 4 + 1] = 1;
  A[2 * 4 + 2] = 1;
  A[3 * 4 + 3] = 1;
}
#ifdef DFM2_STATIC_LIBRARY
template void delfem2::Mat4_Identity(float A[16]);
template void delfem2::Mat4_Identity(double A[16]);
#endif


// ------------------------------------------------

template<typename REAL>
void delfem2::Rotate_Mat4AffineRodriguez(
    REAL A[16],
    const REAL V[3]) {
  REAL B[16];
  Mat4_AffineRotationRodriguez(B,
                               V[0], V[1], V[2]);
  REAL C[16];
  MatMat4(C,
          B, A);

  for (int i = 0; i < 16; ++i) { A[i] = C[i]; }
}
#ifdef DFM2_STATIC_LIBRARY
template void delfem2::Rotate_Mat4AffineRodriguez(
    float A[16], const float V[3]);
template void delfem2::Rotate_Mat4AffineRodriguez(
    double A[16], const double V[3]);
#endif

// -----------------------------------

template<typename REAL>
void delfem2::Mat4_Rotation_Cartesian(
    REAL mat[16],
    const REAL vec[3]) {
  const REAL sqt = vec[0] * vec[0] + vec[1] * vec[1] + vec[2] * vec[2];
  if (sqt < 1.0e-20) {  // infinitesmal rotation approximation
    // row0
    mat[0 * 4 + 0] = 1;
    mat[0 * 4 + 1] = -vec[2];
    mat[0 * 4 + 2] = +vec[1];
    mat[0 * 4 + 3] = 0;
    // row1
    mat[1 * 4 + 0] = +vec[2];
    mat[1 * 4 + 1] = 1;
    mat[1 * 4 + 2] = -vec[0];
    mat[1 * 4 + 3] = 0;
    // row2
    mat[2 * 4 + 0] = -vec[1];
    mat[2 * 4 + 1] = +vec[0];
    mat[2 * 4 + 2] = 1;
    mat[2 * 4 + 3] = 0;
    // row3
    mat[3 * 4 + 0] = 0;
    mat[3 * 4 + 1] = 0;
    mat[3 * 4 + 2] = 0;
    mat[3 * 4 + 3] = 1;
    return;
  }
  REAL t = std::sqrt(sqt);
  REAL invt = 1 / t;
  REAL n[3] = {vec[0] * invt, vec[1] * invt, vec[2] * invt};
  const REAL c0 = std::cos(t);
  const REAL s0 = std::sin(t);
  // row0
  mat[0 * 4 + 0] = c0 + (1 - c0) * n[0] * n[0];
  mat[0 * 4 + 1] = -n[2] * s0 + (1 - c0) * n[0] * n[1];
  mat[0 * 4 + 2] = +n[1] * s0 + (1 - c0) * n[0] * n[2];
  mat[0 * 4 + 3] = 0;
  // row1
  mat[1 * 4 + 0] = +n[2] * s0 + (1 - c0) * n[1] * n[0];
  mat[1 * 4 + 1] = c0 + (1 - c0) * n[1] * n[1];
  mat[1 * 4 + 2] = -n[0] * s0 + (1 - c0) * n[1] * n[2];
  mat[1 * 4 + 3] = 0;
  // row2
  mat[2 * 4 + 0] = -n[1] * s0 + (1 - c0) * n[2] * n[0];
  mat[2 * 4 + 1] = +n[0] * s0 + (1 - c0) * n[2] * n[1];
  mat[2 * 4 + 2] = c0 + (1 - c0) * n[2] * n[2];
  mat[2 * 4 + 3] = 0;
  // row3
  mat[3 * 4 + 0] = 0;
  mat[3 * 4 + 1] = 0;
  mat[3 * 4 + 2] = 0;
  mat[3 * 4 + 3] = 1;
}
#ifdef DFM2_STATIC_LIBRARY
template void delfem2::Mat4_Rotation_Cartesian(float mat[16], const float vec[3]);
template void delfem2::Mat4_Rotation_Cartesian(double mat[16], const double vec[3]);
#endif

// ----------------------------

template<typename REAL>
void delfem2::Translate_Mat4Affine(
    REAL A[16],
    const REAL V[3]) {
  A[0 * 4 + 3] += V[0];
  A[1 * 4 + 3] += V[1];
  A[2 * 4 + 3] += V[2];
}
#ifdef DFM2_STATIC_LIBRARY
template void delfem2::Translate_Mat4Affine(float A[16], const float V[3]);
template void delfem2::Translate_Mat4Affine(double A[16], const double V[3]);
#endif

// -------------------

DFM2_INLINE void delfem2::Mat4_ScaleRotTrans(
    double m[16],
    double scale,
    const double quat[4],
    const double trans[3]) {
  delfem2::Mat4_Quat(m, quat);
  for (int i = 0; i < 3; ++i) {
    for (int j = 0; j < 3; ++j) {
      m[i * 4 + j] *= scale;
    }
  }
  m[0 * 4 + 3] = trans[0];
  m[1 * 4 + 3] = trans[1];
  m[2 * 4 + 3] = trans[2];
}

// ----------------------------

template<typename REAL>
DFM2_INLINE void delfem2::Mat4_Quat(
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
  r[3] = 0;
  r[4] = xy + zw;
  r[5] = 1 - z2 - x2;
  r[6] = yz - xw;
  r[7] = 0;
  r[8] = zx - yw;
  r[9] = yz + xw;
  r[10] = 1 - x2 - y2;
  r[11] = 0;
  r[12] = 0;
  r[13] = 0;
  r[14] = 0;
  r[15] = 1;
}
#ifdef DFM2_STATIC_LIBRARY
template void delfem2::Mat4_Quat(float r[], const float q[]);
template void delfem2::Mat4_Quat(double r[], const double q[]);
#endif

// --------------------------------

// return transpose matrix of Mat4_Quat
template<typename REAL>
DFM2_INLINE void delfem2::Mat4_QuatConj(
    REAL *r,
    const REAL *q) {
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
  r[1] = xy + zw;
  r[2] = zx - yw;
  r[3] = 0;
  r[4] = xy - zw;
  r[5] = 1 - z2 - x2;
  r[6] = yz + xw;
  r[7] = 0;
  r[8] = zx + yw;
  r[9] = yz - xw;
  r[10] = 1 - x2 - y2;
  r[11] = 0;
  r[12] = 0;
  r[13] = 0;
  r[14] = 0;
  r[15] = 1;
}
#ifdef DFM2_STATIC_LIBRARY
template DFM2_INLINE void delfem2::Mat4_QuatConj(float *r, const float *q);
template DFM2_INLINE void delfem2::Mat4_QuatConj(double *r, const double *q);
#endif

/*
// maybe same as above ?
template <typename REAL>
DFM2_INLINE void delfem2::Mat4_AffineTransQuat(
    REAL r[],
    const REAL q[])
{
  REAL x2 = q[0] * q[0] * 2;
  REAL y2 = q[1] * q[1] * 2;
  REAL z2 = q[2] * q[2] * 2;
  REAL xy = q[0] * q[1] * 2;
  REAL yz = q[1] * q[2] * 2;
  REAL zx = q[2] * q[0] * 2;
  REAL xw = q[0] * q[3] * 2;
  REAL yw = q[1] * q[3] * 2;
  REAL zw = q[2] * q[3] * 2;
  // column 0
  r[ 0] = 1 - y2 - z2;
  r[ 1] = xy + zw;
  r[ 2] = zx - yw;
  r[ 3] = 0;
  // column 1
  r[ 4] = xy - zw;
  r[ 5] = 1 - z2 - x2;
  r[ 6] = yz + xw;
  r[ 7] = 0;
  // column 2
  r[ 8] = zx + yw;
  r[ 9] = yz - xw;
  r[10] = 1 - x2 - y2;
  r[11] = 0;
  // column 3
  r[12] = 0;
  r[13] = 0;
  r[14] = 0;
  r[15] = 1;
}
#ifdef DFM2_STATIC_LIBRARY
template void delfem2::Mat4_AffineTransQuat(float r[], const float q[]);
template void delfem2::Mat4_AffineTransQuat(double r[], const double q[]);
#endif
*/

/*
void delfem2::MatMat4(
    double m01[16],
    const double m0[16],
    const double m1[16])
{
  for(int i=0;i<4;++i){
    for(int j=0;j<4;++j){
      m01[i*4+j] = m0[i*4+0]*m1[0*4+j] + m0[i*4+1]*m1[1*4+j] + m0[i*4+2]*m1[2*4+j] + m0[i*4+3]*m1[3*4+j];
    }
  }
}

void delfem2::Copy_Mat4(double m1[16], const double m0[16])
{
  for(int i=0;i<16;++i){ m1[i] = m0[i]; }
}
 */




template<typename REAL>
DFM2_INLINE void delfem2::Inverse_Mat4(
    REAL minv[16],
    const REAL m[16]) {
  for (int i = 0; i < 16; ++i) { minv[i] = m[i]; }
  int info;
  mat4::CalcInvMat(minv, 4, info);
  if (info != 0) {
    for (int i = 0; i < 16; ++i) { minv[i] = m[i]; }
    unsigned int tmp[4];
    mat4::CalcInvMatPivot(minv, 4, tmp);
  }
}
#ifdef DFM2_STATIC_LIBRARY
template void delfem2::Inverse_Mat4(float minv[], const float m[]);
template void delfem2::Inverse_Mat4(double minv[], const double m[]);
#endif

// ------------------------------------------------------------------

template<typename T>
delfem2::CMat4<T> delfem2::CMat4<T>::MatMat(const CMat4<T> &mat0) const {
  CMat4 m;
  ::delfem2::MatMat4(m.mat,
                     this->mat, mat0.mat);
  return m;
}
#ifdef DFM2_STATIC_LIBRARY
template delfem2::CMat4<float> delfem2::CMat4<float>::MatMat(
    const CMat4<float> &mat0) const;
template delfem2::CMat4<double> delfem2::CMat4<double>::MatMat(
    const CMat4<double> &mat0) const;
#endif

// -----------------------------------

template<typename REAL>
delfem2::CMat4<REAL> delfem2::CMat4<REAL>::Quat(const REAL *q) {
  CMat4<REAL> m;
  const REAL x2 = q[0] * q[0] * 2;
  const REAL y2 = q[1] * q[1] * 2;
  const REAL z2 = q[2] * q[2] * 2;
  const REAL xy = q[0] * q[1] * 2;
  const REAL yz = q[1] * q[2] * 2;
  const REAL zx = q[2] * q[0] * 2;
  const REAL xw = q[0] * q[3] * 2;
  const REAL yw = q[1] * q[3] * 2;
  const REAL zw = q[2] * q[3] * 2;
  m.setZero();
  m.mat[0 * 4 + 0] = 1 - y2 - z2;
  m.mat[0 * 4 + 1] = xy - zw;
  m.mat[0 * 4 + 2] = zx + yw;
  m.mat[1 * 4 + 0] = xy + zw;
  m.mat[1 * 4 + 1] = 1 - z2 - x2;
  m.mat[1 * 4 + 2] = yz - xw;
  m.mat[2 * 4 + 0] = zx - yw;
  m.mat[2 * 4 + 1] = yz + xw;
  m.mat[2 * 4 + 2] = 1 - x2 - y2;
  m.mat[3 * 4 + 3] = 1;
  return m;
}
#ifdef DFM2_STATIC_LIBRARY
template delfem2::CMat4<float> delfem2::CMat4<float>::Quat(const float *q);
template delfem2::CMat4<double> delfem2::CMat4<double>::Quat(const double *q);
#endif

// ---------------------------

namespace delfem2 {

template<typename T>
CMat4<T> operator*(const CMat4<T> &lhs, const CMat4<T> &rhs) {
  CMat4<T> q;
  MatMat4(q.mat, lhs.mat, rhs.mat);
  return q;
}
#ifdef DFM2_STATIC_LIBRARY
template CMat4d operator*(const CMat4d &lhs, const CMat4d &rhs);
template CMat4f operator*(const CMat4f &lhs, const CMat4f &rhs);
#endif

template<typename T>
CMat4<T> operator-(const CMat4<T> &lhs, const CMat4<T> &rhs) {
  CMat4<T> q;
  for (int i = 0; i < 16; ++i) { q.mat[i] = lhs.mat[i] - rhs.mat[i]; }
  return q;
}
#ifdef DFM2_STATIC_LIBRARY
template CMat4d operator-(const CMat4d &lhs, const CMat4d &rhs);
template CMat4f operator-(const CMat4f &lhs, const CMat4f &rhs);
#endif

template<typename T>
CMat4<T> operator+(const CMat4<T> &lhs, const CMat4<T> &rhs) {
  CMat4<T> q;
  for (int i = 0; i < 16; ++i) { q.mat[i] = lhs.mat[i] + rhs.mat[i]; }
  return q;
}
#ifdef DFM2_STATIC_LIBRARY
template CMat4d operator+(const CMat4d &lhs, const CMat4d &rhs);
template CMat4f operator+(const CMat4f &lhs, const CMat4f &rhs);
#endif

}

// --------------------------------------------

template<typename REAL>
delfem2::CMat4<REAL> delfem2::CMat4<REAL>::Inverse() const {
  CMat4<REAL> m;
  std::memcpy(m.mat, mat, sizeof(REAL) * 16);
  Inverse_Mat4(m.mat, this->mat);
  return m;
}
#ifdef DFM2_STATIC_LIBRARY
template delfem2::CMat4d delfem2::CMat4d::Inverse() const;
template delfem2::CMat4f delfem2::CMat4f::Inverse() const;
#endif
