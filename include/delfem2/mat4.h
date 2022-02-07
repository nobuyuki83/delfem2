/*
 * Copyright (c) 2019 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */

/**
 * @file 4x4 matrix class (CMat4) and functions
 * @detail The order of dependency in delfem2 is "vec2 -> mat2 -> vec3 -> quaternion -> mat3 -> mat4",
 */

#ifndef DFM2_MAT4_H
#define DFM2_MAT4_H

#include <vector>
#include <cassert>
#include <cmath>
#include <iostream>
#include <limits>  // using NaN Check
#include <array>

#include "delfem2/dfm2_inline.h"
#include "delfem2/geo_meta_funcs.h"

#define NEARLY_ZERO 1.e-16

// -----------------------------

namespace delfem2 {

// if matrix is delfem2::CMat4 or Eigen::Mat4
template<typename MAT, typename T = typename MAT::Scalar>
void Print_Mat4(const MAT &m) {
  std::cout << m(0,0) << " " << m(0,1) << " " << m(0,2) << " " << m(0,3) << std::endl;
  std::cout << m(1,0) << " " << m(1,1) << " " << m(1,2) << " " << m(1,3) << std::endl;
  std::cout << m(2,0) << " " << m(2,1) << " " << m(2,2) << " " << m(2,3) << std::endl;
  std::cout << m(3,0) << " " << m(3,1) << " " << m(3,2) << " " << m(3,3) << std::endl;
}

template<typename MAT, typename std::enable_if_t<std::is_array_v<MAT>, std::nullptr_t> = nullptr>
void Print_Mat4(const MAT &m) {
  std::cout << m[0] << " " << m[1] << " " << m[2] << " " << m[3] << std::endl;
  std::cout << m[4] << " " << m[5] << " " << m[6] << " " << m[7] << std::endl;
  std::cout << m[8] << " " << m[9] << " " << m[10] << " " << m[11] << std::endl;
  std::cout << m[12] << " " << m[13] << " " << m[14] << " " << m[15] << std::endl;
}

template<typename REAL>
void Print_Mat4Transp(const REAL m[16]) {
  std::cout << m[0] << " " << m[4] << " " << m[8] << " " << m[12] << std::endl;
  std::cout << m[1] << " " << m[5] << " " << m[9] << " " << m[13] << std::endl;
  std::cout << m[2] << " " << m[6] << " " << m[10] << " " << m[14] << std::endl;
  std::cout << m[3] << " " << m[7] << " " << m[11] << " " << m[15] << std::endl;
}

template<typename T0, typename T1>
void Copy_Mat4(
    T0 M0[16],
    const T1 M1[16]) {
  for (int i = 0; i < 16; ++i) { M0[i] = static_cast<T0>(M1[i]); }
}

template<typename T0, typename T1>
void Mat4_Mat3(
    T0 B[16],
    T1 K[9]) {
  B[0] = K[0];
  B[1] = K[1];
  B[2] = K[2];
  B[3] = 0;
  B[4] = K[3];
  B[5] = K[4];
  B[6] = K[5];
  B[7] = 0;
  B[8] = K[6];
  B[9] = K[7];
  B[10] = K[8];
  B[11] = 0;
  B[12] = 0;
  B[13] = 0;
  B[14] = 0;
  B[15] = 1;
}

template<typename REAL>
void Transpose_Mat4(
    REAL M0[16],
    const REAL M1[16]) {
  for (int i = 0; i < 4; ++i) {
    for (int j = 0; j < 4; ++j) {
      M0[i * 4 + j] = M1[j * 4 + i];
    }
  }
}

template<typename T>
void Mat4_AffineTrans_RotTransl(
    T A[16],
    const T R[9],
    const T t[3]) {
  A[0 * 4 + 0] = R[0 * 3 + 0];
  A[0 * 4 + 1] = R[1 * 3 + 0];
  A[0 * 4 + 2] = R[2 * 3 + 0];
  A[0 * 4 + 3] = 0;
  A[1 * 4 + 0] = R[0 * 3 + 1];
  A[1 * 4 + 1] = R[1 * 3 + 1];
  A[1 * 4 + 2] = R[2 * 3 + 1];
  A[1 * 4 + 3] = 0;
  A[2 * 4 + 0] = R[0 * 3 + 2];
  A[2 * 4 + 1] = R[1 * 3 + 2];
  A[2 * 4 + 2] = R[2 * 3 + 2];
  A[2 * 4 + 3] = 0;
  A[3 * 4 + 0] = t[0];
  A[3 * 4 + 1] = t[1];
  A[3 * 4 + 2] = t[2];
  A[3 * 4 + 3] = 1;
}

template<typename REAL>
DFM2_INLINE void Inverse_Mat4(
    REAL minv[16],
    const REAL m[16]);

template<typename T0, typename T1, typename T2>
void MatMat4(
    T0 *C,
    const T1 *A,
    const T2 *B) {
  for (int i = 0; i < 4; i++) {
    for (int j = 0; j < 4; j++) {
      const T1 c =
          A[i * 4 + 0] * B[0 * 4 + j] +
              A[i * 4 + 1] * B[1 * 4 + j] +
              A[i * 4 + 2] * B[2 * 4 + j] +
              A[i * 4 + 3] * B[3 * 4 + j];
      C[i * 4 + j] = static_cast<T0>(c);
    }
  }
}

// ---------------------
// below: functions to generate Mat4

template<typename REAL>
void Mat4_Identity(
    REAL A[16]);

/**
 * @brief affine matrix for orthogonal projection
 * @details column major order (fortran order)
 */
template<typename T>
void Mat4_AffineProjectionOrtho(
    T mP[16],
    T xmin, T xmax,
    T ymin, T ymax,
    T zmin, T zmax);

template<typename REAL>
void Mat4_AffineLookAt(
    REAL *Mr,
    REAL eyex, REAL eyey, REAL eyez,
    REAL cntx, REAL cnty, REAL cntz,
    REAL upx, REAL upy, REAL upz);

/**
 * construct projection matrix mapping perspective view frustrum to a cube [-1,+1, -1,+1, -1,+1]
 * The view is from the origin to the -Z direction (zmin < zmax < 0).
 * @param[out] matrix affine matrix (column major)
 * @param[in] fovyInRad filed-of-view in radian
 * @param[in] aspectRatio aspect ratio of the window
 * @param[in] zmin minimum Z coordinate for the view frustrum (mapped to the plane Z==-1)
 * @param[in] zmax maximum Z coordinate for the view frustrum (mapped to the plane Z==+1)
 */
template<typename REAL>
void Mat4_AffineProjectionFrustum(
    REAL mP[16],
    REAL fovyInRad,
    REAL aspectRatio,
    REAL zmin,
    REAL zmax);

template<typename REAL>
void Mat4_AffineScale(
    REAL A[16],
    REAL s);

template<typename REAL>
void Mat4_AffineTranslation(
    REAL A[16],
    REAL dx,
    REAL dy,
    REAL dz);

template<typename REAL>
std::array<REAL,16> Mat4_AffineTranslation(
    REAL dx,
    REAL dy,
    REAL dz);

template<typename T>
std::array<T,16> Mat4_AffineRotationRodriguez(
    T dx,
    T dy,
    T dz);

template<typename REAL>
void Mat4_AffineRotationCartesian(
    REAL mat[16],
    const REAL vec[3]);

/**
 * @func multiply rotation affine matrix from left to an affine matrix in 3D
 * @details the ritation is parmeterized with a rodriguez angle
 */
template<typename REAL>
void Rotate_Mat4AffineRodriguez(
    REAL A[16],
    const REAL v[3]);


// ------------------------

template<typename T>
DFM2_INLINE void MatVec4(
    T v[4],
    const T A[16],
    const T x[4]);

template<typename T>
DFM2_INLINE void VecMat4(
    T v[4],
    const T x[4],
    const T A[16]);

// --------------------------------
// below: functions mat4 and vec3

/**
 * @details this is not affine transfromation. Translation part is ignored
 */
template<typename T>
DFM2_INLINE void Mat4Vec3(
    T vo[3],
    const T M[16],
    const T vi[3]);

DFM2_INLINE void Vec3Mat4(
    double vo[3],
    const double vi[3],
    const double M[16]);

/**
 * @brief multiply translation affine matrix from left to an affine matrix in 3D
 */
template<typename REAL>
void Translate_Mat4Affine(
    REAL A[16],
    const REAL v[3]);

template<typename T1, typename T2>
std::array<T1,3> Vec3_Mat4Vec3_Homography(
    const T1 a[16],
    const T2 x0[3]);

template<typename T>
DFM2_INLINE std::array<T, 2> Vec2_Mat4Vec3_Homography(
    const T a[16],
    const T x0[3]);

template<typename T0, typename T1, typename T2>
void Vec3_Vec3Mat4_Homography(
    T0 y0[3],
    const T1 x0[3],
    const T2 a[16]);

template<typename T>
void Vec3_Mat4Vec3_Affine(
    T y0[3],
    const T a[16],
    const T x0[3]);

// ------------------------------------
// below: function with mat4 and quarternion

template<typename REAL>
DFM2_INLINE void Mat4_AffineQuaternion(
    REAL r[],
    const REAL q[]);

/**
 *
 * @tparam REAL
 * @param r
 * @param q quaternion (x,y,z,w)
 */
template<typename REAL>
DFM2_INLINE void Mat4_AffineQuaternionConjugate(
    REAL *r,
    const REAL *q);

/*
template <typename REAL>
DFM2_INLINE void Mat4_AffineTransQuat(
    REAL r[],
    const REAL q[]);
*/

DFM2_INLINE void Mat4_ScaleRotTrans(
    double m[16],
    double scale,
    const double quat[4],
    const double trans[3]);

// =============================================

template<typename T>
class CMat4;

/**
 * matrix-matrix product
 * @tparam T float and double is defined for the static library
 */
template<typename T>
CMat4<T> operator*(const CMat4<T> &lhs, const CMat4<T> &rhs);

template<typename T>
CMat4<T> operator-(const CMat4<T> &lhs, const CMat4<T> &rhs);

template<typename T>
CMat4<T> operator+(const CMat4<T> &lhs, const CMat4<T> &rhs);

/**
 * @brief 4x4 matrix class
 * @class 4x4 matrix class
 * @tparam REAL value type of the matrix. defiend for "double" and "float" for the static library.
 */
template<typename REAL>
class CMat4 {
 public:
  CMat4() : mat{1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1} {}

  template<typename S>
  explicit CMat4(const S *pm) {
    for (int i = 0; i < 16; ++i) { mat[i] = static_cast<S>(pm[i]); }
  }

  template<typename S>
  CMat4(const std::array<S, 16> &&pm) {
    for (int i = 0; i < 16; ++i) { mat[i] = static_cast<S>(pm[i]); }
  }

  CMat4(REAL v00, REAL v01, REAL v02, REAL v03,
        REAL v10, REAL v11, REAL v12, REAL v13,
        REAL v20, REAL v21, REAL v22, REAL v23,
        REAL v30, REAL v31, REAL v32, REAL v33)
      : mat{
      v00, v01, v02, v03,
      v10, v11, v12, v13,
      v20, v21, v22, v23,
      v30, v31, v32, v33} {}
   
  REAL *data() { return mat; }
  const REAL *data() const { return mat; }
 public:
  inline REAL operator()(int i, int j) const {
    assert(i < 4 && j < 4);
    return mat[i * 4 + j];
  }

  inline REAL &operator()(int i, int j) {
    assert(i < 4 && j < 4);
    return mat[i * 4 + j];
  }
  // ------------------------
  // below: "set" functions
  void SetAffineTranslate(REAL x, REAL y, REAL z) {
    Mat4_AffineTranslation(mat,
                           x, y, z);
  }
  void SetQuaternion(const REAL *q) {
    Mat4_Quat(mat,
              q);
  }

  /**
   * @details named same as Eigen
   */
  void setZero() {
    for (auto &v: mat) { v = 0; }
  }
  void SetIdentity() {
    for (auto &v: mat) { v = 0; }
    mat[0 * 4 + 0] = 1;
    mat[1 * 4 + 1] = 1;
    mat[2 * 4 + 2] = 1;
    mat[3 * 4 + 3] = 1;
  }
  void SetScale(REAL x, REAL y, REAL z) {
    this->setZero();
    mat[0 * 4 + 0] = x;
    mat[1 * 4 + 1] = y;
    mat[2 * 4 + 2] = z;
    mat[3 * 4 + 3] = 1;
  }
  // -----------------------
  template<typename S>
  void CopyTo(S *v) const {
    for (int i = 0; i < 16; ++i) { v[i] = mat[i]; }
  }

  std::array<REAL, 9> GetMat3() const {
    return {
        mat[0], mat[1], mat[2],
        mat[4], mat[5], mat[6],
        mat[8], mat[9], mat[10]};
  }

  std::array<REAL, 16> GetStlArray() const {
    return {
        mat[0], mat[1], mat[2], mat[3],
        mat[4], mat[5], mat[6], mat[7],
        mat[8], mat[9], mat[10],mat[11],
        mat[12],mat[13],mat[14],mat[15]};
  }

  std::array<REAL,3> MultVec3_Homography(const REAL* v) const {
    return Vec3_Mat4Vec3_Homography(mat, v);
  }

  std::array<REAL,2> Vec2_MultVec3_Homography(const REAL* v) const {
    return Vec2_Mat4Vec3_Homography(mat, v);
  }
  
  std::array<REAL,3> MultVec3(const REAL* v) const {
    std::array<REAL,3> r;
    Mat4Vec3(r.data(), mat, v);
    return r;
  }

  /**
 * @details named same as Eigen
 */
  CMat4<REAL> transpose() const {
    CMat4<REAL> m1;
    for (int i = 0; i < 4; ++i) {
      for (int j = 0; j < 4; ++j) {
        m1.mat[i * 4 + j] = mat[j * 4 + i];
      }
    }
    return m1;
  }
  CMat4<REAL> Inverse() const;

  template<typename S>
  CMat4<S> cast() const {
    return CMat4<S>(
        static_cast<S>(mat[0]),
        static_cast<S>(mat[1]),
        static_cast<S>(mat[2]),
        static_cast<S>(mat[3]),
        static_cast<S>(mat[4]),
        static_cast<S>(mat[5]),
        static_cast<S>(mat[6]),
        static_cast<S>(mat[7]),
        static_cast<S>(mat[8]),
        static_cast<S>(mat[9]),
        static_cast<S>(mat[10]),
        static_cast<S>(mat[11]),
        static_cast<S>(mat[12]),
        static_cast<S>(mat[13]),
        static_cast<S>(mat[14]),
        static_cast<S>(mat[15]));
  };
  // ----------------------
  // below: static function
  [[nodiscard]] static CMat4<REAL> Identity() {
    return CMat4<REAL>{
        1, 0, 0, 0,
        0, 1, 0, 0,
        0, 0, 1, 0,
        0, 0, 0, 1};
  }
  [[nodiscard]] static CMat4<REAL> AffineScale(REAL s) {
    return CMat4<REAL>{
        s, 0, 0, 0,
        0, s, 0, 0,
        0, 0, s, 0,
        0, 0, 0, 1};
  }
  [[nodiscard]] static CMat4<REAL> ScaleXYZ(REAL sx, REAL sy, REAL sz) {
    return CMat4<REAL>{
        sx, 0, 0, 0,
        0, sy, 0, 0,
        0, 0, sz, 0,
        0, 0, 0, 1};
  }
  [[nodiscard]] static CMat4<REAL> Spin(const REAL *v) {
    return CMat4<REAL>{
        0, -v[2], +v[1], 0,
        +v[2], 0, -v[0], 0,
        -v[1], +v[0], 0, 0,
        0, 0, 0, 1.0};
  }
  [[nodiscard]] static CMat4<REAL> Quat(const REAL *q);

  [[nodiscard]] static CMat4<REAL> Translation(const std::array<REAL,3>&& v) {
    return CMat4<REAL>{
        1, 0, 0, static_cast<REAL>(v[0]),
        0, 1, 0, static_cast<REAL>(v[1]),
        0, 0, 1, static_cast<REAL>(v[2]),
        0, 0, 0, 1 };
  }

  template<typename S>
  [[nodiscard]] static CMat4<REAL> Translation(const S *v) {
    return CMat4<REAL>{
        1, 0, 0, static_cast<REAL>(v[0]),
        0, 1, 0, static_cast<REAL>(v[1]),
        0, 0, 1, static_cast<REAL>(v[2]),
        0, 0, 0, 1
    };
  }

  [[nodiscard]] static CMat4<REAL> Mat3(const REAL *v) {
    return CMat4<REAL>{
        v[0 * 3 + 0], v[0 * 3 + 1], v[0 * 3 + 2], 0,
        v[1 * 3 + 0], v[1 * 3 + 1], v[1 * 3 + 2], 0,
        v[2 * 3 + 0], v[2 * 3 + 1], v[2 * 3 + 2], 0,
        0, 0, 0, 1};
  }

  [[nodiscard]] static CMat4<REAL> AffineOrthogonalProjection(
      REAL xmin, REAL xmax,
      REAL ymin, REAL ymax,
      REAL zmin, REAL zmax){
    CMat4<REAL> m;
    Mat4_AffineProjectionOrtho(m.data(), xmin, xmax, ymin, ymax, zmin, zmax);
    return m;
  }

  [[nodiscard]] static CMat4<REAL> AffineAxisTransform(
      const std::array<REAL,3> &&x_axis,
      const std::array<REAL,3> &&y_axis,
      const std::array<REAL,3> &&z_axis){
    return CMat4<REAL>{
      x_axis[0], y_axis[0], z_axis[0], 0,
      x_axis[1], y_axis[1], z_axis[1], 0,
      x_axis[2], y_axis[2], z_axis[2], 0,
      0, 0, 0, 1
    };
  }

 public:
  REAL mat[16];
  using Scalar = REAL;
};
using CMat4d = CMat4<double>;
using CMat4f = CMat4<float>;

}

#ifndef DFM2_STATIC_LIBRARY
#  include "delfem2/mat4.cpp"
#endif

#endif
