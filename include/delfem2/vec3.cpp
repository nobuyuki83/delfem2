/*
 * Copyright (c) 2019 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */

#include "delfem2/vec3.h"

#include <cmath>
#include <stack>

#include "delfem2/geo_vec3.h"

#ifndef M_PI
#  define M_PI 3.14159265358979323846
#endif

// =====================================
// below: unexposed 

namespace delfem2::vec3 {

DFM2_INLINE bool MyIsnan(double x) { return x != x; }

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

// ------------------------------------------

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
template void delfem2::AverageFour3(
    float po[3],
    const float p0[3], const float p1[3], const float p2[3], const float p3[3]);
template void delfem2::AverageFour3(
    double po[3],
    const double p0[3], const double p1[3], const double p2[3], const double p3[3]);
#endif


// above: without CVec3 (move to "geo_vec3_raw.cpp"?)
// ======================================================
// below: with CVec



// ---------------------

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
template<typename T0, typename T1>
CVec3<T0> operator/(const CVec3<T0> &vec, T1 d) {
  CVec3<T0> temp = vec;
  temp /= d;
  return temp;
}
#ifdef DFM2_STATIC_LIBRARY
template CVec3f operator/ (const CVec3f& vec, float d);
template CVec3f operator/ (const CVec3f& vec, int d);
template CVec3d operator/ (const CVec3d& vec, double d);
template CVec3d operator/ (const CVec3d& vec, int d);
#endif

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


// -----------------------------------------------------------

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

// ----------------------------------------------------------------------------

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

template<typename REAL>
std::array<REAL, 9> delfem2::Mat3_MinimumRotation(
    const CVec3<REAL> &V,
    const CVec3<REAL> &v) {
  CVec3<REAL> ep = V.normalized();
  CVec3<REAL> eq = v.normalized();
  CVec3<REAL> n = ep.cross(eq);
  const REAL st2 = n.dot(n);
  const REAL ct = ep.dot(eq);

  if (st2 < 1.0e-8f) { // very small angle or n is zero
    // inifinitesimal rotation
    if (ct > 0.99) {
      return {
          1.f + 0.5f * (n.x * n.x - st2),
          -n.z + 0.5f * (n.x * n.y),
          +n.y + 0.5f * (n.x * n.z),
          +n.z + 0.5f * (n.y * n.x),
          1.f + 0.5f * (n.y * n.y - st2),
          -n.x + 0.5f * (n.y * n.z),
          -n.y + 0.5f * (n.z * n.x),
          +n.x + 0.5f * (n.z * n.y),
          1.f + 0.5f * (n.z * n.z - st2)};
    } else {
      CVec3<REAL> epx, epy;
      delfem2::FrameFromVectorZ(epx, epy, ep);
      const CVec3<REAL> eqx = epx - eq.dot(epx) * eq; // vector orthogonal to eq
      const CVec3<REAL> eqy = eq.cross(eqx);
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
#ifdef DFM2_STATIC_LIBRARY
template std::array<double,9> delfem2::Mat3_MinimumRotation(const CVec3d& V, const CVec3d& v);
template std::array<float,9> delfem2::Mat3_MinimumRotation(const CVec3f& V, const CVec3f& v);
#endif


// --------------------------

template<typename REAL>
DFM2_INLINE std::array<REAL, 9> delfem2::Mat3_ParallelTransport(
    const CVec3<REAL> &p0,
    const CVec3<REAL> &p1,
    const CVec3<REAL> &q0,
    const CVec3<REAL> &q1) {
  return Mat3_MinimumRotation(p1 - p0, q1 - q0);
}
