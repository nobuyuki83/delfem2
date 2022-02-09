/*
 * Copyright (c) 2019 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */

/**
 * @file 3x3 matrix class (CMat3) and functions
 * @detail The order of dependency in delfem2 is "vec2 -> mat2 -> vec3 -> quaternion -> mat3 -> mat4",
 */


#ifndef DFM2_MAT3_FUNCS_H
#define DFM2_MAT3_FUNCS_H

#include <vector>
#include <cassert>
#include <cmath>
#include <iostream>
#include <array>
#include <limits> // using NaN Check

#include "delfem2/dfm2_inline.h"
#include "delfem2/geo_meta_funcs.h"

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
    REAL alpha,
    REAL beta);

template<typename REAL>
void Mat3_Identity(
    REAL *mat,
    REAL alpha);

template<typename REAL>
DFM2_INLINE void Mat3_Identity_ScaleAdd(
    REAL *mat,
    REAL alpha = 1,
    REAL beta = 0);

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

// -----------------
// below: axis angle vector

template<typename T>
DFM2_INLINE void AxisAngleVectorCartesian_Mat3(
    T v[3],
    const T m[9]);

template<typename REAL>
void Mat3_Rotation_Cartesian(
    REAL mat[9],
    const REAL vec[3]);

template<typename T>
DFM2_INLINE void AxisAngleVectorCRV_Mat3(
    T crv[3],
    const T mat[9]);

template<typename T>
DFM2_INLINE void EulerAngle_Mat3(
    T ea[3],
    const T m[9],
    const std::array<int, 3> &axis_idxs);

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

template<typename VEC3, typename T = value_type<VEC3>>
std::array<T,3> Mat3Vec3(const T mat[9], const VEC3 &v) {
  return {
    mat[0]*v[0] + mat[1]*v[1] + mat[2]*v[2],
    mat[3]*v[0] + mat[4]*v[1] + mat[5]*v[2],
    mat[6]*v[0] + mat[7]*v[1] + mat[8]*v[2] };
}

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
DFM2_INLINE void Vec2_Mat3Vec2_Homography(
    T y[2],
    const T Z[9],
    const T x[2]);

template<typename T>
std::array<T,2> Vec2_Mat3Vec2_Homography(
    const T Z[9],
    const T x[2]);

template<typename T>
DFM2_INLINE void Vec2_Mat3Vec2_AffineDirection(
    T y[2],
    const T A[9],
    const T x[2]);

/**
 * @brief 3x3 Rotation matrix to rotate V into v with minimum rotation angle
 * this functions is for dfm2::CVec3 or Eigen::Vector3
 * @param[in] V rotation from
 * @param[in] v rotation to
 */
template<typename VEC, typename REAL= typename VEC::Scalar>
std::array<REAL, 9> Mat3_MinimumRotation(
    const VEC &V,
    const VEC &v);

// --------------------------------
// below: mat3 and quat

template<typename REAL>
DFM2_INLINE void Mat3_Quat(
    REAL r[],
    const REAL q[]);

template<typename T>
DFM2_INLINE void Quat_Mat3(
    T quat[4],
    const T p_[9]);

} // namespace delfem2


#ifndef DFM2_STATIC_LIBRARY
#  include "delfem2/mat3_funcs.cpp"
#endif

#endif /* DFM2_MAT3_FUNCS_H */
