/*
 * Copyright (c) 2019 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */

/**
 * @file 3x3 matrix class (CMat3) and functions
 */


#ifndef DFM2_MAT3_H
#define DFM2_MAT3_H

#include <vector>
#include <cassert>
#include <cmath>
#include <iostream>
#include <limits> // using NaN Check

#include "delfem2/dfm2_inline.h"

#define NEARLY_ZERO 1.e-16

// -----------------------------

namespace delfem2 {

/**
 * @func Set spin 3x3 matrix (skew asymetric matrix)
 * @tparam REAL float and double
 * @param mat pointer of 3x3 matrix
 * @param v ponter of 3 vector
 */
template<typename REAL>
void Mat3_Spin(
    REAL *mat,
    const REAL *v);

template<typename REAL>
void Mat3_Spin_ScaleAdd(
    REAL *mat,
    const REAL *v,
    REAL alpha, REAL beta);

template<typename REAL>
void Mat3_Identity(
    REAL *mat,
    REAL alpha);

template<typename REAL>
DFM2_INLINE void Mat3_Identity_ScaleAdd(
    REAL *mat,
    REAL alpha = 1, REAL beta = 0);

template<typename T>
DFM2_INLINE void Mat3_AffineRotation(
    T *mat,
    T theta);

template<typename T>
void Mat3_AffineTranslation(
    T *mat,
    const T transl[2]);

template<typename T0, typename T1>
void Copy_Mat3(
    T0 m0[9],
    const T1 m1[9]) {
  for (int i = 0; i < 9; ++i) { m0[i] = m1[i]; }
}

// ------------

template<typename T>
void Transpose_Mat3(
    T At[],
    const T A[]);

template<typename T0, typename T1>
void Inverse_Mat3(
    T0 Ainv[],
    const T1 A[]);

template<typename REAL>
void Inverse_Mat3(
    REAL Ainv[9]);

/**
 * compute C = A*B
 * @tparam T0
 * @tparam T1
 * @tparam T2
 * @param C
 * @param A
 * @param B
 */
template<typename T0, typename T1, typename T2>
void MatMat3(
    T0 *AB,
    const T1 *A,
    const T2 *B);

template<typename T0, typename T1, typename T2>
void MatMatT3(
    T0 *ABt,
    const T1 *A,
    const T2 *B);

/**
 * @func product of a transposed 3x3 matrix and another 3x3 matrix.
 * [C] = [A]^T[B}
 * @details row major data structure
 */
template<typename T0, typename T1, typename T2>
void MatTMat3(
    T0 *AtB,
    const T1 *A,
    const T2 *B);

/**
 * @func adding scaled product of a transposed 3x3 matrix and another 3x3 matrix.
 * [C] = alpha * [A]^T[B} + beta* [C]
 * @details row major data structure
 */
template<typename T>
void MatTMat3_ScaleAdd(
    T *C,
    const T *A,
    const T *B,
    T alpha,
    T beta);

template<typename T>
T Det_Mat3(
    const T U[9]);

template<typename T>
T SquareNormFrobenius_SymMat3(
    const T sm[6]);

template<typename REAL>
void Mat3_Rotation_Cartesian(
    REAL mat[9],
    const REAL vec[3]);

/**
 * @func compute eigen value & vector for symmmetric matrix
 * @details
 * sm[6] = (M_00,M_11,M_22,M_12,M_20,M_01)
 * M = ULU^T
 * u[9] = (U_00,U_01,U_02, U_10,U_11,U_12, U_20,U_21,U_22)
 */
DFM2_INLINE bool eigenSym3(
    double u[9],
    double l[3],
    const double sm[6],
    int nitr);

DFM2_INLINE void svd3(
    double U[9],
    double G[3],
    double V[9],
    const double m[9],
    int nitr);

DFM2_INLINE void GetRotPolarDecomp(
    double R[9],
    const double am[9],
    int nitr);

template<typename T>
DFM2_INLINE void AxisAngleVectorCartesian_Mat3(
    T v[3],
    const T m[9]);

template<typename T>
DFM2_INLINE void AxisAngleVectorCRV_Mat3(
    T crv[3],
    const T mat[9]);

// ------------------------------------------------
// below: mat3 and vec3

template<typename T0, typename T1, typename T2>
DFM2_INLINE void MatTVec3(
    T0 y[3],
    const T1 m[9],
    const T2 x[3]);

/**
 * @func {y} = beta*{y} + alpha*[M]^T{x}
 */
template<typename T0, typename T1, typename T2, typename T3, typename T4>
void MatTVec3_ScaleAdd(
    T0 y[3],
    const T1 m[9],
    const T2 x[3],
    T3 alpha,
    T4 beta);

/**
 * @func matrix vector product for 3x3 matrix {y} := [m]{x}
 */
template<typename T0, typename T1, typename T2>
void MatVec3(
    T0 y[3],
    const T1 m[9],
    const T2 x[3]);

template<typename T>
void MatVec3_ScaleAdd(
    T y[3],
    const T m[9],
    const T x[3],
    T alpha,
    T beta);

DFM2_INLINE void VecMat3(
    double y[3],
    const double x[3],
    const double m[9]);

template<typename T>
DFM2_INLINE void Vec2_Mat3Vec2_AffineProjection(
    T y[2],
    const T Z[9],
    const T x[2]);

template<typename T>
DFM2_INLINE void Vec2_Mat3Vec2_AffineDirection(
    T y[2],
    const T A[9],
    const T x[2]);

// --------------------------------
// below: mat3 and quat

template<typename REAL>
DFM2_INLINE void Mat3_Quat(
    REAL r[],
    const REAL q[]);

template<typename T>
class CMat3; // this pre-definition is needed for following functions

template<typename T>
CMat3<T> operator+(const CMat3<T> &lhs, const CMat3<T> &rhs);

template<typename T>
CMat3<T> operator-(const CMat3<T> &lhs, const CMat3<T> &rhs);

template<typename T>
CMat3<T> operator*(double d, const CMat3<T> &rhs);

template<typename T>
CMat3<T> operator*(const CMat3<T> &m, T d);

template<typename T>
CMat3<T> operator*(const CMat3<T> &lhs, const CMat3<T> &rhs);

template<typename T>
CMat3<T> operator/(const CMat3<T> &m, T d);

template<typename T>
std::ostream &operator<<(std::ostream &output, const CMat3<T> &m);

template<typename T>
std::istream &operator>>(std::istream &output, CMat3<T> &m);

static inline bool myIsNAN_Matrix3(double d) { return !(d > d - 1); }

/**
 * @class class of 3x3 matrix
 */
template<typename REAL>
class CMat3 {
 public:
  CMat3();
  explicit CMat3(REAL s);
  explicit CMat3(const REAL m[9]);
  CMat3(REAL v00, REAL v01, REAL v02,
        REAL v10, REAL v11, REAL v12,
        REAL v20, REAL v21, REAL v22);
  CMat3(REAL x, REAL y, REAL z);
  // ---------------
  REAL *data() { return p_; }
  [[nodiscard]] const REAL *data() const { return p_; }
  // ---------------
  void GetElements(REAL m[9]) const { for (unsigned int i = 0; i < 9; i++) { m[i] = p_[i]; }}
  [[nodiscard]] double Get(int i, int j) const { return p_[i * 3 + j]; }
  // ---------
  void AffineMatrixTrans(REAL m[16]) const {
    m[0 * 4 + 0] = p_[0];
    m[1 * 4 + 0] = p_[1];
    m[2 * 4 + 0] = p_[2];
    m[3 * 4 + 0] = 0;
    m[0 * 4 + 1] = p_[3];
    m[1 * 4 + 1] = p_[4];
    m[2 * 4 + 1] = p_[5];
    m[3 * 4 + 1] = 0;
    m[0 * 4 + 2] = p_[6];
    m[1 * 4 + 2] = p_[7];
    m[2 * 4 + 2] = p_[8];
    m[3 * 4 + 2] = 0;
    m[0 * 4 + 3] = 0;
    m[1 * 4 + 3] = 0;
    m[2 * 4 + 3] = 0;
    m[3 * 4 + 3] = 1;
  }
  void CopyToMat4(REAL m[16]) const {
    m[0 * 4 + 0] = p_[0];
    m[0 * 4 + 1] = p_[1];
    m[0 * 4 + 2] = p_[2];
    m[1 * 4 + 0] = p_[3];
    m[1 * 4 + 1] = p_[4];
    m[1 * 4 + 2] = p_[5];
    m[2 * 4 + 0] = p_[6];
    m[2 * 4 + 1] = p_[7];
    m[2 * 4 + 2] = p_[8];
  }
  void CopyTo(REAL *ptr) const {
    for (int i = 0; i < 9; ++i) { ptr[i] = p_[i]; }
  }
  void CopyToScale(REAL *ptr, REAL s) const {
    for (int i = 0; i < 9; ++i) { ptr[i] = p_[i] * s; }
  }
  void AddToScale(REAL *ptr, REAL s) const {
    for (int i = 0; i < 9; ++i) { ptr[i] += p_[i] * s; }
  }
  // ---------------
//  CVector3 MatVec(const CVector3& vec0) const;
  void MatVec(const REAL vec0[], REAL vec1[]) const;
  void MatVecTrans(const REAL vec0[], REAL vec1[]) const;

//  CVector3 MatVecTrans(const CVector3& vec0) const;
  [[nodiscard]] CMat3 MatMat(const CMat3 &mat0) const;
  [[nodiscard]] CMat3 MatMatTrans(const CMat3 &mat0) const;
  // ----------------
  [[nodiscard]] CMat3 Sym() const {
    CMat3 m;
    for (unsigned int i = 0; i < 3; i++) {
      for (unsigned int j = 0; j < 3; j++) {
        m.p_[i * 3 + j] = (p_[i * 3 + j] + p_[j * 3 + i]) * 0.5;
      }
    }
    return m;
  }
  inline CMat3 operator-() const { return (*this) * static_cast<REAL>(-1); }
  inline CMat3 operator+() const { return (*this); }
  inline CMat3 &operator+=(const CMat3 &rhs) {
    for (unsigned int i = 0; i < 9; i++) { p_[i] += rhs.p_[i]; }
    return *this;
  }
  inline CMat3 &operator-=(const CMat3 &rhs) {
    for (unsigned int i = 0; i < 9; i++) { p_[i] -= rhs.p_[i]; }
    return *this;
  }
  inline CMat3 &operator*=(REAL d) {
    for (auto &m : p_) { m *= d; }
    return *this;
  }
  inline CMat3 &operator/=(REAL d) {
    REAL invd = (REAL) 1.0 / d;
    for (auto &m : p_) { m *= invd; }
    return *this;
  }
  inline REAL operator[](int i) const {
    return this->p_[i];
  }
  inline REAL &operator()(int i, int j) {
    return this->p_[i * 3 + j];
  }
  inline const REAL &operator()(int i, int j) const {
    return this->p_[i * 3 + j];
  }
  // -------------------------
  [[nodiscard]] CMat3 Inverse() const;
  // -------------------------
  // function whose name starts with "Set" changes itself
  void SetInverse();
  void SetSymetric(const REAL sm[6]);
  void setZero();
  void SetRandom();
  void SetRotMatrix_Cartesian(const REAL vec[]);
  void SetRotMatrix_Cartesian(REAL x, REAL y, REAL z);
  void SetRotMatrix_Rodrigues(const REAL vec[]);
  void SetRotMatrix_CRV(const REAL crv[]);
  void SetRotMatrix_Quaternion(const REAL quat[]);
  void SetRotMatrix_BryantAngle(REAL rx, REAL ry, REAL rz);
  void SetIdentity(REAL scale = 1);
  void SetMat4(const REAL m[16]) {
    for (int i = 0; i < 3; ++i) {
      for (int j = 0; j < 3; ++j) {
        this->p_[i * 3 + j] = m[i * 4 + j];
      }
    }
  }
  // ------------------------
  void GetQuat_RotMatrix(REAL quat[]) const;
  // ------------------------
  [[nodiscard]] CMat3 transpose() const {
    CMat3 m;
    m.p_[0] = p_[0];
    m.p_[1] = p_[3];
    m.p_[2] = p_[6];
    m.p_[3] = p_[1];
    m.p_[4] = p_[4];
    m.p_[5] = p_[7];
    m.p_[6] = p_[2];
    m.p_[7] = p_[5];
    m.p_[8] = p_[8];
    return m;
  }
  [[nodiscard]] bool isNaN() const {
    double s = p_[0] + p_[1] + p_[2] + p_[3] + p_[4] + p_[5] + p_[6] + p_[7] + p_[8];
    return myIsNAN_Matrix3(s) != 0;
  }
  [[nodiscard]] double determinant() const {
    return
        p_[0] * p_[4] * p_[8] +
        p_[3] * p_[7] * p_[2] +
        p_[6] * p_[1] * p_[5] -
        p_[0] * p_[7] * p_[5] -
        p_[6] * p_[4] * p_[2] -
        p_[3] * p_[1] * p_[8];
  }
  [[nodiscard]] double SqNorm_Frobenius() const {
    double s = 0.0;
    for (auto &i : p_) {
      s += i * i;
    }
    return s;
  }
  [[nodiscard]] double Trace() const {
    return p_[0] + p_[4] + p_[8];
  }
  [[nodiscard]] double SecondInvarint() const {
    const CMat3 &m2 = (*this) * (*this);
    const double tr = this->Trace();
    return 0.5 * (tr * tr - m2.Trace());
  }
  void Print() const {
    std::cout << p_[0] << " " << p_[1] << " " << p_[2] << std::endl;
    std::cout << p_[3] << " " << p_[4] << " " << p_[5] << std::endl;
    std::cout << p_[6] << " " << p_[7] << " " << p_[8] << std::endl;
  }
  void PolerDecomp(CMat3 &R, int nitr) const {
    GetRotPolarDecomp(R.p_,
                      p_, nitr);
  }
  // --------------------
  // static functions
  static CMat3 Identity(REAL scale = 1) {
    CMat3 m;
    m.SetIdentity(scale);
    return m;
  }
  static CMat3 Zero() {
    CMat3 m;
    m.setZero();
    return m;
  }
  static CMat3 Spin(const REAL *v) {
    CMat3 m;
    Mat3_Spin(m.p_, v);
    return m;
  }
  static CMat3 OuterProduct(const REAL *v0, const REAL *v1) {
    return CMat3<REAL>(
        v0[0] * v1[0], v0[0] * v1[1], v0[0] * v1[2],
        v0[1] * v1[0], v0[1] * v1[1], v0[1] * v1[2],
        v0[2] * v1[0], v0[2] * v1[1], v0[2] * v1[2]);
  }
  // quaternion order of (x,y,z,w)
  static CMat3 Quat(const REAL *q) {
    const REAL x2 = q[0] * q[0] * 2;
    const REAL y2 = q[1] * q[1] * 2;
    const REAL z2 = q[2] * q[2] * 2;
    const REAL xy = q[0] * q[1] * 2;
    const REAL yz = q[1] * q[2] * 2;
    const REAL zx = q[2] * q[0] * 2;
    const REAL xw = q[0] * q[3] * 2;
    const REAL yw = q[1] * q[3] * 2;
    const REAL zw = q[2] * q[3] * 2;
    CMat3<REAL> m;
    m.p_[0 * 3 + 0] = 1 - y2 - z2;
    m.p_[0 * 3 + 1] = xy - zw;
    m.p_[0 * 3 + 2] = zx + yw;
    m.p_[1 * 3 + 0] = xy + zw;
    m.p_[1 * 3 + 1] = 1 - z2 - x2;
    m.p_[1 * 3 + 2] = yz - xw;
    m.p_[2 * 3 + 0] = zx - yw;
    m.p_[2 * 3 + 1] = yz + xw;
    m.p_[2 * 3 + 2] = 1 - x2 - y2;
    return m;
  }
 public:
  REAL p_[9]; // value with row-major order
};

using CMat3d = CMat3<double>;
using CMat3f = CMat3<float>;

template<typename T>
CMat3<T> Mat3_Identity(T alpha) {
  CMat3<T> m;
  Mat3_Identity(m.p_, alpha);
  return m;
}

}

#ifndef DFM2_STATIC_LIBRARY
#  include "delfem2/mat3.cpp"
#endif

#endif /* DFM2_MAT3_H */
