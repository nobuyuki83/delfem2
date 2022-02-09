/*
 * Copyright (c) 2019 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */

#include "delfem2/mat3.h"

#include <random>

#include "delfem2/mat3_funcs.h"

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
std::array<T,3> delfem2::CMat3<T>::MatVec(const T vec0[3]) const {
  std::array<T,3> vec1;
  ::delfem2::MatVec3(vec1.data(), p_, vec0);
  return vec1;
}
#ifdef DFM2_STATIC_LIBRARY
template std::array<float,3> delfem2::CMat3f::MatVec(const float vec0[]) const;
template std::array<double,3> delfem2::CMat3d::MatVec(const double vec0[]) const;
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
  Mat3_Quat(p_, quat);
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
std::array<T,4> delfem2::CMat3<T>::GetQuaternion() const {
  std::array<T,4> q;
  Quat_Mat3(q.data(), p_);
  return q;
}
#ifdef DFM2_STATIC_LIBRARY
template std::array<float,4> delfem2::CMat3<float>::GetQuaternion() const;
template std::array<double,4> delfem2::CMat3<double>::GetQuaternion() const;
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

