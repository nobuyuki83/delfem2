/*
 * Copyright (c) 2019 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */

#include "delfem2/vec3.h"

#include <cmath>
#include <stack>

#ifndef M_PI
#  define M_PI 3.14159265358979323846
#endif

// =====================================
// below: unexposed 

namespace delfem2::vec3 {

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

// there is another impelemntation in quat.h so this is "static function"
// transform vector with quaternion
// quaternion order (x,y,z,w)
template<typename REAL>
DFM2_INLINE void MyQuatVec(
    REAL vo[],
    const REAL q[],
    const REAL vi[]) {
  REAL x2 = q[0] * q[0] * 2;
  REAL y2 = q[1] * q[1] * 2;
  REAL z2 = q[2] * q[2] * 2;
  REAL xy = q[0] * q[1] * 2;
  REAL yz = q[1] * q[2] * 2;
  REAL zx = q[2] * q[0] * 2;
  REAL xw = q[0] * q[3] * 2;
  REAL yw = q[1] * q[3] * 2;
  REAL zw = q[2] * q[3] * 2;
  vo[0] = (1 - y2 - z2) * vi[0] + (xy - zw) * vi[1] + (zx + yw) * vi[2];
  vo[1] = (xy + zw) * vi[0] + (1 - z2 - x2) * vi[1] + (yz - xw) * vi[2];
  vo[2] = (zx - yw) * vi[0] + (yz + xw) * vi[1] + (1 - x2 - y2) * vi[2];
}
#ifdef DFM2_STATIC_LIBRARY
template void MyQuatVec(float vo[], const float q[], const float vi[]);
template void MyQuatVec(double vo[], const double q[], const double vi[]);
#endif

// ----------------------

template<typename REAL>
DFM2_INLINE void MyMat4Vec3(
    REAL vo[3],
    const REAL M[16], const REAL vi[3]) {
  vo[0] = M[0 * 4 + 0] * vi[0] + M[0 * 4 + 1] * vi[1] + M[0 * 4 + 2] * vi[2];
  vo[1] = M[1 * 4 + 0] * vi[0] + M[1 * 4 + 1] * vi[1] + M[1 * 4 + 2] * vi[2];
  vo[2] = M[2 * 4 + 0] * vi[0] + M[2 * 4 + 1] * vi[1] + M[2 * 4 + 2] * vi[2];
}
#ifdef DFM2_STATIC_LIBRARY
template void MyMat4Vec3(float vo[3],
    const float M[16], const float vi[3]);
template void MyMat4Vec3(double vo[3],
    const double M[16], const double vi[3]);
#endif

// ----------------------

// there is formal implementation in quat.cpp so this is static to avoid dumplicated
// quaternion order (x,y,z,w)
template<typename REAL>
DFM2_INLINE void MyQuatConjVec(
    REAL vo[3],
    const REAL q[4],
    const REAL vi[3]) {
  const REAL x2 = q[0] * q[0] * 2;
  const REAL y2 = q[1] * q[1] * 2;
  const REAL z2 = q[2] * q[2] * 2;
  const REAL xy = q[0] * q[1] * 2;
  const REAL yz = q[1] * q[2] * 2;
  const REAL zx = q[2] * q[0] * 2;
  const REAL xw = q[0] * q[3] * 2;
  const REAL yw = q[1] * q[3] * 2;
  const REAL zw = q[2] * q[3] * 2;
  vo[0] = (1 - y2 - z2) * vi[0] + (xy + zw) * vi[1] + (zx - yw) * vi[2];
  vo[1] = (xy - zw) * vi[0] + (1 - z2 - x2) * vi[1] + (yz + xw) * vi[2];
  vo[2] = (zx + yw) * vi[0] + (yz - xw) * vi[1] + (1 - x2 - y2) * vi[2];
}
#ifdef DFM2_STATIC_LIBRARY
template void MyQuatConjVec(float vo[3], const float q[4], const float vi[3]);
template void MyQuatConjVec(double vo[3], const double q[4], const double vi[3]);
#endif

// --------------------------

template<typename REAL>
DFM2_INLINE void MyInverse_Mat3(
    REAL Ainv[9],
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

}

// ===========================================================

namespace delfem2 { // template specialization need to be done in the namespace

template<>
DFM2_INLINE double Distance3
    (const double p0[3],
     const double p1[3]) {
  return sqrt(
      (p1[0] - p0[0]) * (p1[0] - p0[0]) +
          (p1[1] - p0[1]) * (p1[1] - p0[1]) +
          (p1[2] - p0[2]) * (p1[2] - p0[2]));
}

template<>
DFM2_INLINE float Distance3
    (const float p0[3],
     const float p1[3]) {
  return sqrtf(
      (p1[0] - p0[0]) * (p1[0] - p0[0]) +
          (p1[1] - p0[1]) * (p1[1] - p0[1]) +
          (p1[2] - p0[2]) * (p1[2] - p0[2]));
}

}

// -------------------------

namespace delfem2 { // template specialization need to be done in the namespace

template<>
DFM2_INLINE double Length3(const double v[3]) {
  return sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);
}

template<>
DFM2_INLINE float Length3(const float v[3]) {
  return sqrtf(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);
}

}

// ---------------------------------

template<typename T>
DFM2_INLINE T delfem2::Dot3(const T a[3], const T b[3]) {
  return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}
#ifdef DFM2_STATIC_LIBRARY
template float delfem2::Dot3(const float a[3], const float b[3]);
template double delfem2::Dot3(const double a[3], const double b[3]);
#endif


// ---------------------------------

template<typename REAL>
DFM2_INLINE void delfem2::Add3(
    REAL vo[3],
    const REAL vi[3]) {
  vo[0] += vi[0];
  vo[1] += vi[1];
  vo[2] += vi[2];
}
#ifdef DFM2_STATIC_LIBRARY
template void delfem2::Add3(float vo[3], const float vi[3]);
template void delfem2::Add3(double vo[3], const double vi[3]);
#endif

// ---------------------------------

// ---------------------------------

template<typename T0, typename T1, typename T2>
void delfem2::Cross3(
    T0 r[3],
    const T1 v1[3],
    const T2 v2[3]) {
  r[0] = v1[1] * v2[2] - v2[1] * v1[2];
  r[1] = v1[2] * v2[0] - v2[2] * v1[0];
  r[2] = v1[0] * v2[1] - v2[0] * v1[1];
}
#ifdef DFM2_STATIC_LIBRARY
template void delfem2::Cross3(float r[3], const float v1[3], const float v2[3]);
template void delfem2::Cross3(double r[3], const double v1[3], const double v2[3]);
#endif


// ---------------------------------

template<typename T>
DFM2_INLINE void delfem2::Normalize3(T v[3]) {
  T len = std::sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);
  v[0] /= len;
  v[1] /= len;
  v[2] /= len;
}
#ifdef DFM2_STATIC_LIBRARY
template void delfem2::Normalize3(float v[3]);
template void delfem2::Normalize3(double v[3]);
#endif

// ---------------------------------

template<typename T>
T delfem2::Area_Tri3(const T v1[3], const T v2[3], const T v3[3]) {
  const T n[3] = {
      (v2[1] - v1[1]) * (v3[2] - v1[2]) - (v3[1] - v1[1]) * (v2[2] - v1[2]),
      (v2[2] - v1[2]) * (v3[0] - v1[0]) - (v3[2] - v1[2]) * (v2[0] - v1[0]),
      (v2[0] - v1[0]) * (v3[1] - v1[1]) - (v3[0] - v1[0]) * (v2[1] - v1[1])};
  return Length3(n) / 2;
}
#ifdef DFM2_STATIC_LIBRARY
template float delfem2::Area_Tri3(const float v1[3], const float v2[3], const float v3[3]);
template double delfem2::Area_Tri3(const double v1[3], const double v2[3], const double v3[3]);
#endif

// ----------------------------------

template<typename T>
DFM2_INLINE T delfem2::SquareDistance3(const T p0[3], const T p1[3]) {
  return
      (p1[0] - p0[0]) * (p1[0] - p0[0]) +
          (p1[1] - p0[1]) * (p1[1] - p0[1]) +
          (p1[2] - p0[2]) * (p1[2] - p0[2]);
}
#ifdef DFM2_STATIC_LIBRARY
template float delfem2::SquareDistance3(const float p0[3], const float p1[3]);
template double delfem2::SquareDistance3(const double p0[3], const double p1[3]);
#endif

// ------------------------------------

template<typename T>
DFM2_INLINE T delfem2::SquareLength3(const T v[3]) {
  return v[0] * v[0] + v[1] * v[1] + v[2] * v[2];
}
#ifdef DFM2_STATIC_LIBRARY
template float delfem2::SquareLength3(const float v[3]);
template double delfem2::SquareLength3(const double v[3]);
#endif

// ------------------------------------

template<typename T>
T delfem2::ScalarTripleProduct3(const T a[], const T b[], const T c[]) {
  return a[0] * (b[1] * c[2] - b[2] * c[1])
      + a[1] * (b[2] * c[0] - b[0] * c[2])
      + a[2] * (b[0] * c[1] - b[1] * c[0]);
}
#ifdef DFM2_STATIC_LIBRARY
template float delfem2::ScalarTripleProduct3(const float a[], const float b[], const float c[]);
template double delfem2::ScalarTripleProduct3(const double a[], const double b[], const double c[]);
#endif

// -----------------------------------------

template<typename REAL>
void delfem2::NormalTri3(
    REAL n[3],
    const REAL v1[3],
    const REAL v2[3],
    const REAL v3[3]) {
  n[0] = (v2[1] - v1[1]) * (v3[2] - v1[2]) - (v3[1] - v1[1]) * (v2[2] - v1[2]);
  n[1] = (v2[2] - v1[2]) * (v3[0] - v1[0]) - (v3[2] - v1[2]) * (v2[0] - v1[0]);
  n[2] = (v2[0] - v1[0]) * (v3[1] - v1[1]) - (v3[0] - v1[0]) * (v2[1] - v1[1]);
}
#ifdef DFM2_STATIC_LIBRARY
template void delfem2::NormalTri3(float n[3], const float v1[3], const float v2[3], const float v3[3]);
template void delfem2::NormalTri3(double n[3], const double v1[3], const double v2[3], const double v3[3]);
#endif

// ------------------------------------------

template<typename REAL>
void delfem2::UnitNormalAreaTri3(
    REAL n[3],
    REAL &a,
    const REAL v1[3], const REAL v2[3], const REAL v3[3]) {
  NormalTri3(n,
             v1, v2, v3);
  a = Length3(n) / 2;
  const REAL invlen = 1 / (a * 2);
  n[0] *= invlen;
  n[1] *= invlen;
  n[2] *= invlen;
}
#ifdef DFM2_STATIC_LIBRARY
template void delfem2::UnitNormalAreaTri3(float n[3], float& a,
    const float v1[3], const float v2[3], const float v3[3]);
template void delfem2::UnitNormalAreaTri3(double n[3], double& a,
    const double v1[3], const double v2[3], const double v3[3]);
#endif

// ------------------------------------------

DFM2_INLINE void delfem2::GetVertical2Vector3D(
    const double vec_n[3],
    double vec_x[3],
    double vec_y[3]) {
  const double vec_s[3] = {0, 1, 0};
  Cross3(vec_x, vec_s, vec_n);
  const double len = Length3(vec_x);
  if (len < 1.0e-10) {
    const double vec_t[3] = {1, 0, 0};
    Cross3(vec_x, vec_t, vec_n);  // z????
    Cross3(vec_y, vec_n, vec_x);  // x????
  } else {
    const double invlen = 1.0 / len;
    vec_x[0] *= invlen;
    vec_x[1] *= invlen;
    vec_x[2] *= invlen;
    Cross3(vec_y, vec_n, vec_x);
  }
}

// ---------------------------------------------------------------------------

template<typename T>
T delfem2::Dot(const CVec3<T> &arg1, const CVec3<T> &arg2) {
  return Dot3(arg1.p, arg2.p);
}
#ifdef DFM2_STATIC_LIBRARY
template float delfem2::Dot(const CVec3f& arg1, const CVec3f& arg2);
template double delfem2::Dot(const CVec3d& arg1, const CVec3d& arg2);
#endif

// ----------------------------

// cross product
template<typename T>
delfem2::CVec3<T> delfem2::Cross(const CVec3<T> &arg1, const CVec3<T> &arg2) {
  CVec3<T> temp;
  Cross3(temp.p, arg1.p, arg2.p);
  return temp;
}
#ifdef DFM2_STATIC_LIBRARY
template delfem2::CVec3f delfem2::Cross(const CVec3f& arg1, const CVec3f& arg2);
template delfem2::CVec3d delfem2::Cross(const CVec3d& arg1, const CVec3d& arg2);
#endif

// ----------------------------

template<typename REAL>
void delfem2::AverageTwo3(
    REAL po[3],
    const REAL p0[3], const REAL p1[3]) {
  constexpr REAL half = static_cast<REAL>(0.5);
  po[0] = (p0[0] + p1[0]) * half;
  po[1] = (p0[1] + p1[1]) * half;
  po[2] = (p0[2] + p1[2]) * half;
}
#ifdef DFM2_STATIC_LIBRARY
template void delfem2::AverageTwo3(float po[3], const float p0[3], const float p1[3]);
template void delfem2::AverageTwo3(double po[3], const double p0[3], const double p1[3]);
#endif

// ----------------------------

template<typename REAL>
void delfem2::AverageFour3(
    REAL po[3],
    const REAL p0[3],
    const REAL p1[3],
    const REAL p2[3],
    const REAL p3[3]) {
  constexpr REAL quarter(0.25);
  po[0] = (p0[0] + p1[0] + p2[0] + p3[0]) * quarter;
  po[1] = (p0[1] + p1[1] + p2[1] + p3[1]) * quarter;
  po[2] = (p0[2] + p1[2] + p2[2] + p3[2]) * quarter;
}
#ifdef DFM2_STATIC_LIBRARY
template void delfem2::AverageFour3(float po[3],
    const float p0[3], const float p1[3], const float p2[3], const float p3[3]);
template void delfem2::AverageFour3(double po[3],
    const double p0[3], const double p1[3], const double p2[3], const double p3[3]);
#endif



// above: without CVec3
// --------------------------------------------------------------------------------------------
// below: with CVec3

namespace delfem2 {

//! scale
template<typename T0, typename T1>
CVec3<T0> operator*(T1 d, const CVec3<T0> &rhs) {
  return CVec3<T0>(
      static_cast<T0>(rhs.x * d),
      static_cast<T0>(rhs.y * d),
      static_cast<T0>(rhs.z * d));
}
#ifdef DFM2_STATIC_LIBRARY
template CVec3f operator* (float d, const CVec3f& rhs);
template CVec3d operator* (double d, const CVec3d& rhs);
template CVec3d operator* (int d, const CVec3d& rhs);
#endif

// -----------------------

//! divide by real number
template<typename T>
CVec3<T> operator/(const CVec3<T> &vec, T d) {
  CVec3<T> temp = vec;
  temp /= d;
  return temp;
}
#ifdef DFM2_STATIC_LIBRARY
template CVec3f operator/ (const CVec3f& vec, float d);
template CVec3d operator/ (const CVec3d& vec, double d);
#endif

// ----------------


// ------------------

template<typename T>
std::ostream &operator<<(std::ostream &output, const CVec3<T> &v) {
  output.setf(std::ios::scientific);
  output << v.p[0] << " " << v.p[1] << " " << v.p[2];
  return output;
}
#ifdef DFM2_STATIC_LIBRARY
template std::ostream &operator<<(std::ostream &output, const CVec3d& v);
template std::ostream &operator<<(std::ostream &output, const CVec3f& v);
#endif

// ---------------------

template<typename T>
std::istream &operator>>(std::istream &input, CVec3<T> &v) {
  input >> v.p[0] >> v.p[1] >> v.p[2];
  return input;
}
#ifdef DFM2_STATIC_LIBRARY
template std::istream &operator>>(std::istream &input, CVec3d& v);
template std::istream &operator>>(std::istream &input, CVec3f& v);
#endif

} // namespace delfem2

// ----------------------

// ------------------------------------------------------------------------------------------------------------------------------

template<typename T>
delfem2::CVec3<T> delfem2::Mat3Vec(const T mat[9], const CVec3<T> &v) {
  CVec3<T> u;
  vec3::MyMatVec3(u.p, mat, v.p);
  return u;
}
#ifdef DFM2_STATIC_LIBRARY
template delfem2::CVec3<float> delfem2::Mat3Vec(
    const float mat[9], const CVec3<float>& v);
template delfem2::CVec3<double> delfem2::Mat3Vec(
    const double mat[9], const CVec3<double>& v);
#endif

// -------------------------

template<typename T>
delfem2::CVec3<T> delfem2::Mat4Vec(const T mat[16], const CVec3<T> &v) {
  CVec3<T> u;
  vec3::MyMat4Vec3(u.p, mat, v.p);
  return u;
}
#ifdef DFM2_STATIC_LIBRARY
template delfem2::CVec3f delfem2::Mat4Vec(const float mat[16], const CVec3f& v);
template delfem2::CVec3d delfem2::Mat4Vec(const double mat[16], const CVec3d& v);
#endif

// ------------------------

template<typename T>
DFM2_INLINE delfem2::CVec3<T> delfem2::QuatVec(
    const T quat[4],
    const CVec3<T> &v0) {
  T v1a[3];
  vec3::MyQuatVec(v1a, quat, v0.p);
  return CVec3<T>(v1a[0], v1a[1], v1a[2]);
}
#ifdef DFM2_STATIC_LIBRARY
template delfem2::CVec3f delfem2::QuatVec (const float quat[4], const CVec3f& v0);
template delfem2::CVec3d delfem2::QuatVec (const double quat[4], const CVec3d& v0);
#endif

// ----------------------------

template<typename REAL>
delfem2::CVec3<REAL> delfem2::QuatConjVec
    (const REAL quat[4],
     const CVec3<REAL> &v0) {
  REAL v1a[3];
  vec3::MyQuatConjVec(v1a,
                      quat, v0.p);
  return CVec3<REAL>(v1a[0], v1a[1], v1a[2]);
}
#ifdef DFM2_STATIC_LIBRARY
template delfem2::CVec3f delfem2::QuatConjVec(const float quat[4], const CVec3f& v0);
template delfem2::CVec3d delfem2::QuatConjVec(const double quat[4], const CVec3d& v0);
#endif


// ------------------------------------------------------------

template<typename T>
T delfem2::ScalarTripleProduct(
    const CVec3<T> &a,
    const CVec3<T> &b,
    const CVec3<T> &c) {
//  return a.x*(b.y*c.z - b.z*c.y) + a.y*(b.z*c.x - b.x*c.z) + a.z*(b.x*c.y - b.y*c.x);
  return
  a.x * (b.y * c.z - b.z * c.y) +
  a.y * (b.z * c.x - b.x * c.z) +
  a.z * (b.x * c.y - b.y * c.x);
}
#ifdef DFM2_STATIC_LIBRARY
template double delfem2::ScalarTripleProduct(
    const CVec3<double> &a,
    const CVec3<double> &b,
    const CVec3<double> &c);
template float delfem2::ScalarTripleProduct(
    const CVec3<float> &a,
    const CVec3<float> &b,
    const CVec3<float> &c);
#endif

// --------------------------

namespace delfem2 {

template<typename T>
bool operator==(const CVec3<T> &lhs, const CVec3<T> &rhs) {
  if (fabs(lhs.p[0] - rhs.p[0]) < NEARLY_ZERO
      && fabs(lhs.p[1] - rhs.p[1]) < NEARLY_ZERO
      && fabs(lhs.p[2] - rhs.p[2]) < NEARLY_ZERO) { return true; }
  else { return false; }
}

template<typename T>
bool operator!=(const CVec3<T> &lhs, const CVec3<T> &rhs) {
  return !(lhs == rhs);
}

} // namespace delfem2


// ----------------------------

template<typename T>
void delfem2::CVec3<T>::normalize() {
  T invmag = 1 / norm();
  p[0] *= invmag;
  p[1] *= invmag;
  p[2] *= invmag;
}
#ifdef DFM2_STATIC_LIBRARY
template void delfem2::CVec3<float>::normalize();
template void delfem2::CVec3<double>::normalize();
#endif

// ----------------------------

template<typename T>
void delfem2::CVec3<T>::setZero() {
  p[0] = 0;
  p[1] = 0;
  p[2] = 0;
}
#ifdef DFM2_STATIC_LIBRARY
template void delfem2::CVec3<float>::setZero();
template void delfem2::CVec3<double>::setZero();
template void delfem2::CVec3<int>::setZero();
#endif


// -----------------------
/*
template <typename T, typename DIST>
void delfem2::CVec3<T>::SetRandom(DIST dist)
{
  p[0] = dist();
  p[1] = dist();
  p[2] = dist();
}
#ifdef DFM2_STATIC_LIBRARY
template void delfem2::CVec3<float>::SetRandom();
template void delfem2::CVec3<double>::SetRandom();
#endif
 */

// -------------------------------------------------------------------------

template<typename T>
void delfem2::GetVertical2Vector
    (const CVec3<T> &vec_n,
     CVec3<T> &vec_x,
     CVec3<T> &vec_y) {
  vec_x = Cross(CVec3<T>(0, 1, 0), vec_n);
  const T len = vec_x.norm();
  if (len < 1.0e-10) {
    vec_x = Cross(CVec3<T>(1, 0, 0), vec_n);  // z????
    vec_x.normalize();
    vec_y = Cross(vec_n, vec_x);  // x????
  } else {
    const T invlen = 1 / len;
    vec_x *= invlen;
    vec_y = Cross(vec_n, vec_x);
  }
}
#ifdef DFM2_STATIC_LIBRARY
template void delfem2::GetVertical2Vector(const CVec3f& vec_n,
                                          CVec3f& vec_x,
                                          CVec3f& vec_y);
template void delfem2::GetVertical2Vector(const CVec3d& vec_n,
                                          CVec3d& vec_x,
                                          CVec3d& vec_y);
#endif

// ------------------------------------------------------------------------------

// matrix are column major
template<typename T>
delfem2::CVec3<T> delfem2::mult_GlAffineMatrix
    (const float *m,
     const CVec3<T> &p) {
  CVec3<T> v;
  v.p[0] = m[0 * 4 + 0] * p.p[0] + m[1 * 4 + 0] * p.p[1] + m[2 * 4 + 0] * p.p[2] + m[3 * 4 + 0];
  v.p[1] = m[0 * 4 + 1] * p.p[0] + m[1 * 4 + 1] * p.p[1] + m[2 * 4 + 1] * p.p[2] + m[3 * 4 + 1];
  v.p[2] = m[0 * 4 + 2] * p.p[0] + m[1 * 4 + 2] * p.p[1] + m[2 * 4 + 2] * p.p[2] + m[3 * 4 + 2];
  return v;
}

template<typename T>
delfem2::CVec3<T> delfem2::solve_GlAffineMatrix
    (const float *m,
     const CVec3<T> &p) {
  CVec3<T> v = p - CVec3<T>(m[3 * 4 + 0], m[3 * 4 + 1], m[3 * 4 + 2]);
  double M[9] = {
      m[0 * 4 + 0], m[1 * 4 + 0], m[2 * 4 + 0],
      m[0 * 4 + 1], m[1 * 4 + 1], m[2 * 4 + 1],
      m[0 * 4 + 2], m[1 * 4 + 2], m[2 * 4 + 2]};
  double Minv[9];
  vec3::MyInverse_Mat3(Minv, M);
  return Mat3Vec(Minv, v);
//  CMatrix3 Minv = M.Inverse();  
//  return Minv*v;
}

template<typename T>
delfem2::CVec3<T> delfem2::solve_GlAffineMatrixDirection
    (const float *m,
     const CVec3<T> &v) {
  double M[9] = {
      m[0 * 4 + 0], m[1 * 4 + 0], m[2 * 4 + 0],
      m[0 * 4 + 1], m[1 * 4 + 1], m[2 * 4 + 1],
      m[0 * 4 + 2], m[1 * 4 + 2], m[2 * 4 + 2]};
  double Minv[9];
  vec3::MyInverse_Mat3(Minv, M);
  return Mat3Vec(Minv, v);
  /*
  CMatrix3 M(m[0*4+0],m[1*4+0],m[2*4+0],
             m[0*4+1],m[1*4+1],m[2*4+1],
             m[0*4+2],m[1*4+2],m[2*4+2]);
   */
  /*
   CMatrix3 M(m[0*4+0], m[0*4+1], m[0*4+2],
   m[1*4+0], m[1*4+1], m[1*4+2],
   m[2*4+0], m[2*4+1], m[2*4+2]);
   */
//  CMatrix3 Minv = M.Inverse();
//  return Minv*v;
}

// ----------------------

template<typename T>
delfem2::CVec3<T> delfem2::screenProjection
    (const CVec3<T> &v,
     const float *mMV,
     const float *mPj) {
  CVec3<T> v0 = mult_GlAffineMatrix(mMV, v);
  CVec3<T> v1 = mult_GlAffineMatrix(mPj, v0);
  float w1 = mPj[11] * (float) v0.p[2] + mPj[15];
  return CVec3<T>(v1.p[0] / w1, v1.p[1] / w1, 0.0);
}
#ifdef DFM2_STATIC_LIBRARY
template delfem2::CVec3d delfem2::screenProjection (
    const CVec3d& v,
    const float* mMV,
    const float* mPj);
#endif

template<typename T>
delfem2::CVec3<T> delfem2::screenUnProjection(
    const CVec3<T> &v,
    const float *mMV,
    const float *mPj) {
  float D = mPj[11] + mPj[15]; // z is 1 after model view
  CVec3<T> v0(D * v.p[0], D * v.p[1], 0.0);
  CVec3<T> v1 = solve_GlAffineMatrix(mPj, v0);
  v1.p[2] = 1;
  CVec3<T> v2 = solve_GlAffineMatrix(mMV, v1);
  return v2;
}
#ifdef DFM2_STATIC_LIBRARY
template delfem2::CVec3d delfem2::screenUnProjection(
    const CVec3d& v,
    const float* mMV, const float* mPj);
#endif

template<typename T>
delfem2::CVec3<T> delfem2::screenUnProjectionDirection(
    const CVec3<T> &v,
    const float *mMV,
    const float *mPj) {
  CVec3<T> v0 = solve_GlAffineMatrixDirection(mPj, v);
  CVec3<T> v1 = solve_GlAffineMatrixDirection(mMV, v0);
  v1.normalize();
  return v1;
}
#ifdef DFM2_STATIC_LIBRARY
template delfem2::CVec3d delfem2::screenUnProjectionDirection(
    const CVec3d& v,
    const float* mMV,
    const float* mPj);
#endif

// -----------------------

template<typename T>
delfem2::CVec3<T> delfem2::screenDepthDirection(
    const CVec3<T> &v,
    const float *mMV,
    const float *mPj) {
  float Dv = mPj[11] + mPj[15]; // z is 1 after model view
  CVec3<T> v0(Dv * v.p[0], Dv * v.p[1], 0.0);
  CVec3<T> v1 = solve_GlAffineMatrix(mPj, v0);
  v1.p[2] = 1;
  //
  float Du = mPj[11] * 2.f + mPj[15]; // z is 1 after model view
  CVec3<T> u0(Du * v.p[0], Du * v.p[1], 0.0);
  CVec3<T> u1 = solve_GlAffineMatrix(mPj, u0);
  u1.p[2] = 2;
  //
  CVec3<T> v2 = solve_GlAffineMatrixDirection(mMV, (v1 - u1));
  v2.normalize();
  return v2;
}

// ----------------------------------------------------------------------------

template<typename T>
void delfem2::Cross(
    CVec3<T> &lhs,
    const CVec3<T> &v1,
    const CVec3<T> &v2) {
  lhs.p[0] = v1.p[1] * v2.p[2] - v2.p[1] * v1.p[2];
  lhs.p[1] = v1.p[2] * v2.p[0] - v2.p[2] * v1.p[0];
  lhs.p[2] = v1.p[0] * v2.p[1] - v2.p[0] * v1.p[1];
}
#ifdef DFM2_STATIC_LIBRARY
template void delfem2::Cross(CVec3d& lhs, const CVec3d& v1, const CVec3d& v2 );
#endif

// ----------------------

template<typename T>
T delfem2::Area_Tri(
    const CVec3<T> &v1,
    const CVec3<T> &v2,
    const CVec3<T> &v3) {
  const T x = (v2.p[1] - v1.p[1]) * (v3.p[2] - v1.p[2]) - (v3.p[1] - v1.p[1]) * (v2.p[2] - v1.p[2]);
  const T y = (v2.p[2] - v1.p[2]) * (v3.p[0] - v1.p[0]) - (v3.p[2] - v1.p[2]) * (v2.p[0] - v1.p[0]);
  const T z = (v2.p[0] - v1.p[0]) * (v3.p[1] - v1.p[1]) - (v3.p[0] - v1.p[0]) * (v2.p[1] - v1.p[1]);
  return 0.5 * std::sqrt(x * x + y * y + z * z);
}
#ifdef DFM2_STATIC_LIBRARY
template double delfem2::Area_Tri(const CVec3<double>& v1, const CVec3<double>& v2, const CVec3<double>& v3);
#endif

// ---------------------------

template<typename T>
T delfem2::SquareTriArea(
    const CVec3<T> &v1,
    const CVec3<T> &v2,
    const CVec3<T> &v3) {
  const T dtmp_x = (v2.p[1] - v1.p[1]) * (v3.p[2] - v1.p[2]) - (v2.p[2] - v1.p[2]) * (v3.p[1] - v1.p[1]);
  const T dtmp_y = (v2.p[2] - v1.p[2]) * (v3.p[0] - v1.p[0]) - (v2.p[0] - v1.p[0]) * (v3.p[2] - v1.p[2]);
  const T dtmp_z = (v2.p[0] - v1.p[0]) * (v3.p[1] - v1.p[1]) - (v2.p[1] - v1.p[1]) * (v3.p[0] - v1.p[0]);
  return (dtmp_x * dtmp_x + dtmp_y * dtmp_y + dtmp_z * dtmp_z) * 0.25;
}

// ---------------------

template<typename T>
double delfem2::SquareDistance(
    const CVec3<T> &ipo0,
    const CVec3<T> &ipo1) {
  return (ipo1.p[0] - ipo0.p[0]) * (ipo1.p[0] - ipo0.p[0])
      + (ipo1.p[1] - ipo0.p[1]) * (ipo1.p[1] - ipo0.p[1])
      + (ipo1.p[2] - ipo0.p[2]) * (ipo1.p[2] - ipo0.p[2]);
}
#ifdef DFM2_STATIC_LIBRARY
template double delfem2::SquareDistance(
    const CVec3d& ipo0,
    const CVec3d& ipo1);
template double delfem2::SquareDistance(
    const CVec3f& ipo0,
    const CVec3f& ipo1);
#endif

// ---------------------

template<typename T>
double delfem2::SquareLength(
    const CVec3<T> &point) {
  return point.p[0] * point.p[0] + point.p[1] * point.p[1] + point.p[2] * point.p[2];
}

// ---------------------

//! length of vector
template<typename T>
double delfem2::Length(const CVec3<T> &point) {
  return delfem2::Length3(point.p);
}

//! distance between two points
template<typename T>
double delfem2::Distance(
    const CVec3<T> &p0,
    const CVec3<T> &p1) {
  return delfem2::Distance3(p0.p, p1.p);
}
#ifdef DFM2_STATIC_LIBRARY
template double delfem2::Distance(const CVec3<double>& p0, const CVec3<double>& p1);
#endif

// -------------------------------------------

template<typename T>
void delfem2::Normal(
    CVec3<T> &vnorm,
    const CVec3<T> &v1,
    const CVec3<T> &v2,
    const CVec3<T> &v3) {
  vnorm.p[0] = (v2.p[1] - v1.p[1]) * (v3.p[2] - v1.p[2]) - (v2.p[2] - v1.p[2]) * (v3.p[1] - v1.p[1]);
  vnorm.p[1] = (v2.p[2] - v1.p[2]) * (v3.p[0] - v1.p[0]) - (v2.p[0] - v1.p[0]) * (v3.p[2] - v1.p[2]);
  vnorm.p[2] = (v2.p[0] - v1.p[0]) * (v3.p[1] - v1.p[1]) - (v2.p[1] - v1.p[1]) * (v3.p[0] - v1.p[0]);
}
#ifdef DFM2_STATIC_LIBRARY
template void delfem2::Normal(CVec3d& vnorm, const CVec3d& v1, const CVec3d& v2, const CVec3d& v3);
#endif

// -------------------------------------------

template<typename T>
delfem2::CVec3<T> delfem2::Normal(
    const CVec3<T> &v1,
    const CVec3<T> &v2,
    const CVec3<T> &v3) {
  CVec3<T> vnorm;
  vnorm.p[0] = (v2.p[1] - v1.p[1]) * (v3.p[2] - v1.p[2]) - (v2.p[2] - v1.p[2]) * (v3.p[1] - v1.p[1]);
  vnorm.p[1] = (v2.p[2] - v1.p[2]) * (v3.p[0] - v1.p[0]) - (v2.p[0] - v1.p[0]) * (v3.p[2] - v1.p[2]);
  vnorm.p[2] = (v2.p[0] - v1.p[0]) * (v3.p[1] - v1.p[1]) - (v2.p[1] - v1.p[1]) * (v3.p[0] - v1.p[0]);
  return vnorm;
}
#ifdef DFM2_STATIC_LIBRARY
template delfem2::CVec3d delfem2::Normal(const CVec3d& v1, const CVec3d& v2, const CVec3d& v3);
#endif


// --------------------------------------------

template<typename T>
void delfem2::UnitNormal(
    CVec3<T> &vnorm,
    const CVec3<T> &v1,
    const CVec3<T> &v2,
    const CVec3<T> &v3) {
  vnorm.p[0] = (v2.p[1] - v1.p[1]) * (v3.p[2] - v1.p[2]) - (v2.p[2] - v1.p[2]) * (v3.p[1] - v1.p[1]);
  vnorm.p[1] = (v2.p[2] - v1.p[2]) * (v3.p[0] - v1.p[0]) - (v2.p[0] - v1.p[0]) * (v3.p[2] - v1.p[2]);
  vnorm.p[2] = (v2.p[0] - v1.p[0]) * (v3.p[1] - v1.p[1]) - (v2.p[1] - v1.p[1]) * (v3.p[0] - v1.p[0]);
  const T dtmp1 = 1 / vnorm.norm();
  vnorm.p[0] *= dtmp1;
  vnorm.p[1] *= dtmp1;
  vnorm.p[2] *= dtmp1;
}
#ifdef DFM2_STATIC_LIBRARY
template void delfem2::UnitNormal(CVec3f& vnorm, const CVec3f& v1, const CVec3f& v2, const CVec3f& v3);
template void delfem2::UnitNormal(CVec3d& vnorm, const CVec3d& v1, const CVec3d& v2, const CVec3d& v3);
#endif

// ---------------------------------------

template<typename T>
delfem2::CVec3<T> delfem2::UnitNormal(
    const CVec3<T> &v1,
    const CVec3<T> &v2,
    const CVec3<T> &v3) {
  CVec3<T> vnorm;
  vnorm.p[0] = (v2.p[1] - v1.p[1]) * (v3.p[2] - v1.p[2]) - (v2.p[2] - v1.p[2]) * (v3.p[1] - v1.p[1]);
  vnorm.p[1] = (v2.p[2] - v1.p[2]) * (v3.p[0] - v1.p[0]) - (v2.p[0] - v1.p[0]) * (v3.p[2] - v1.p[2]);
  vnorm.p[2] = (v2.p[0] - v1.p[0]) * (v3.p[1] - v1.p[1]) - (v2.p[1] - v1.p[1]) * (v3.p[0] - v1.p[0]);
  const T dtmp1 = 1 / vnorm.norm();
  vnorm.p[0] *= dtmp1;
  vnorm.p[1] *= dtmp1;
  vnorm.p[2] *= dtmp1;
  return vnorm;
}
#ifdef DFM2_STATIC_LIBRARY
template delfem2::CVec3f delfem2::UnitNormal
 (const CVec3f& v1,
  const CVec3f& v2,
  const CVec3f& v3);
template delfem2::CVec3d delfem2::UnitNormal
(const CVec3d& v1,
 const CVec3d& v2,
 const CVec3d& v3);
#endif

// ---------------------------------------------------

template<typename T>
delfem2::CVec3<T> delfem2::RotateVector(
    const CVec3<T> &vec0,
    const CVec3<T> &rot) {
  const double theta = rot.norm();
  if (theta < 1.0e-30) {
    return vec0;
  }
  CVec3<T> e0 = rot;
  e0.normalize();
  CVec3<T> e2 = delfem2::Cross(e0, vec0);
  if (e2.norm() < 1.0e-30) {
    return vec0;
  }
  e2.normalize();
  CVec3<T> e1 = delfem2::Cross(e2, e0);
  assert(fabs(e1.norm() - 1) < 1.0e-10);
  //	assert( e2.p[0]*vec_0.p[0] + e2.p[1]*vec_0.p[1] + e2.p[2]*vec_0.p[2] < 1.0e-10 );
  const double dot00 = Dot(vec0, e0);
  const double dot01 = Dot(vec0, e1);
  const double cost = cos(theta);
  const double sint = sin(theta);
  CVec3<T> vec1;
  vec1.p[0] = dot00 * e0.p[0] + dot01 * cost * e1.p[0] + dot01 * sint * e2.p[0];
  vec1.p[1] = dot00 * e0.p[1] + dot01 * cost * e1.p[1] + dot01 * sint * e2.p[1];
  vec1.p[2] = dot00 * e0.p[2] + dot01 * cost * e1.p[2] + dot01 * sint * e2.p[2];
  return vec1;
}

template<typename T>
delfem2::CVec3<T> delfem2::RandVector() {
  CVec3<T> r;
  r.p[0] = (2 * (double) rand() / (RAND_MAX + 1.0) - 1);
  r.p[1] = (2 * (double) rand() / (RAND_MAX + 1.0) - 1);
  r.p[2] = (2 * (double) rand() / (RAND_MAX + 1.0) - 1);
  return r;
}

template<typename T>
delfem2::CVec3<T> delfem2::RandUnitVector() {
  for (int itr = 0; itr < 100; itr++) {
    CVec3<T> r = RandVector<T>();
    double l = r.norm();
    if ((l <= 1 || itr == 9) && l > 1.0e-5) {
      r.normalize();
      return r;
    }
  }
  return CVec3<T>(1, 0, 0);
}

template<typename T>
delfem2::CVec3<T> delfem2::RandGaussVector() {
  double a0 = rand() / (RAND_MAX + 1.0);
  double a1 = rand() / (RAND_MAX + 1.0);
  double a2 = rand() / (RAND_MAX + 1.0);
  double a3 = rand() / (RAND_MAX + 1.0);

  double x = sqrt(-2.0 * log(a0)) * cos(3.1415 * 2 * a1);
  double y = sqrt(-2.0 * log(a0)) * sin(3.1415 * 2 * a1);
  double z = sqrt(-2.0 * log(a2)) * cos(3.1415 * 2 * a3);
  return CVec3<T>(x, y, z);
}



// ----------------------------------------------------------------------------------------
// using <vector> from here

// --------------------------------------------------
// TODO: following should be move to mesh class?

namespace delfem2 {

template<typename T>
std::ostream &operator<<(
    std::ostream &output,
    const std::vector<CVec3 < T>>
& aV) {
output<<aV.
size()
<<
std::endl;
for (
int iv = 0;
iv<(int)
aV.
size();
++iv) {
output<<"  "<<aV[iv]<<
std::endl;
}
return
output;
}

template<typename T>
std::istream &operator>>(std::istream &input, std::vector<CVec3 < T>>
& aV){
int nV;
input>>
nV;
aV.
resize(nV);
for (
int iv = 0;
iv<nV;
iv++){
input>>aV[iv];
}
return
input;
}

} // end namespace delfem2

template<typename T>
double delfem2::Area_Tri(
    const int iv1,
    const int iv2,
    const int iv3,
    const std::vector<CVec3 < T>>
&aPoint ) {
return
Area_Tri(aPoint[iv1], aPoint[iv2], aPoint[iv3]
);
}
