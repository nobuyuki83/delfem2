/*
 * Copyright (c) 2019 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */

/**
 * @file 4x4 matrix class (CMat4) and functions
 */


#ifndef DFM2_MAT4_H
#define DFM2_MAT4_H

#include <vector>
#include <cassert>
#include <cmath>
#include <iostream>
#include <limits> // using NaN Check
#include "delfem2/dfm2_inline.h"

#define NEARLY_ZERO 1.e-16

// -----------------------------

namespace delfem2 {

template <typename REAL>
void Print_Mat4(const REAL m[16]){
  std::cout <<  m[ 0] << " " << m[ 1] << " " << m[ 2] << " " << m[ 3] << std::endl;
  std::cout <<  m[ 4] << " " << m[ 5] << " " << m[ 6] << " " << m[ 7] << std::endl;
  std::cout <<  m[ 8] << " " << m[ 9] << " " << m[10] << " " << m[11] << std::endl;
  std::cout <<  m[12] << " " << m[13] << " " << m[14] << " " << m[15] << std::endl;
}

template <typename REAL>
void Print_Mat4Transp(const REAL m[16]){
  std::cout <<  m[ 0] << " " << m[ 4] << " " << m[ 8] << " " << m[12] << std::endl;
  std::cout <<  m[ 1] << " " << m[ 5] << " " << m[ 9] << " " << m[13] << std::endl;
  std::cout <<  m[ 2] << " " << m[ 6] << " " << m[10] << " " << m[14] << std::endl;
  std::cout <<  m[ 3] << " " << m[ 7] << " " << m[11] << " " << m[15] << std::endl;
}

template <typename REAL>
void Copy_Mat4(
    REAL M0[16],
    const REAL M1[16])
{
  for(int i=0;i<16;++i){ M0[i] = M1[i]; }
}


template <typename REAL>
void Transpose_Mat4(
    REAL M0[16],
    const REAL M1[16])
{
  for(int i=0;i<4;++i){
  for(int j=0;j<4;++j){
    M0[i*4+j] = M1[j*4+i];
  }
  }
}

template <typename REAL>
DFM2_INLINE void Inverse_Mat4(
  REAL minv[16],
  const REAL m[16]);


template <typename T>
void MatMat4(
    T* C,
    const T* A,
    const T* B);

/**
 * @brief affine matrix for orthogonal projection
 * @details column major order (fortran order)
 */
template <typename T>
void Mat4_AffineTransProjectionOrtho(
    T mP[16],
    double xmin, double xmax,
    double ymin, double ymax,
    double zmin, double zmax);

void Mat4_AffineTransLookAt(
    float* Mr,
    float eyex, float eyey, float eyez,
    float cntx, float cnty, float cntz,
    float upx, float upy, float upz );

void MultMat4AffineTransTranslateFromRight(
    float *matrix,
    float x,
    float y,
    float z);

/**
 * construct projection matrix mapping perspective view frustrum to a cube [-1,+1, -1,+1, -1,+1]
 * The view is from the origin to the -Z direction (zmin < zmax < 0).
 * @param[out] matrix affine matrix (column major)
 * @param[in] fovyInRad filed-of-view in radian
 * @param[in] aspectRatio aspect ratio of the window
 * @param[in] zmin minimum Z coordinate for the view frustrum (mapped to the plane Z==-1)
 * @param[in] zmax maximum Z coordinate for the view frustrum (mapped to the plane Z==+1)
 */
void Mat4_AffineTransProjectionFrustum(
    float mP[16],
    float fovyInRad,
    float aspectRatio,
    float zmin,
    float zmax);

template <typename REAL>
void Mat4_Identity(
    REAL A[16]);

template <typename REAL>
void Mat4_AffineScale(
    REAL A[16],
    REAL s);

template <typename REAL>
void Mat4_AffineTranslation(
    REAL A[16],
    REAL dx, REAL dy, REAL dz);

template <typename REAL>
void Mat4_AffineTranslation(
    REAL A[16],
    const REAL v[3]);

template <typename REAL>
DFM2_INLINE void Mat4_AffineTransTranslate(
    REAL r[],
    const REAL t[]);

template <typename T>
void Mat4_AffineRotationRodriguez(
    T A[16],
    T dx, T dy, T dz);

/**
 * @func multiply rotation affine matrix from left to an affine matrix in 3D
 * @details the ritation is parmeterized with a rodriguez angle
 */
template <typename REAL>
void Rotate_Mat4AffineRodriguez(
    REAL A[16],
    const REAL v[3]);

// ------------------------

template <typename T>
DFM2_INLINE void MatVec4(
    T v[4],
    const T A[16],
    const T x[4]);

template <typename T>
DFM2_INLINE void VecMat4(
    T v[4],
    const T x[4],
    const T A[16]);

// --------------------------------
// below: functions mat4 and vec3

/**
 * @details this is not affine transfromation. Translation part is ignored
 */
DFM2_INLINE void Mat4Vec3(
     double vo[3],
     const double M[16],
     const double vi[3]);

/**
 * @brief multiply translation affine matrix from left to an affine matrix in 3D
 */
template <typename REAL>
void Translate_Mat4Affine(
    REAL A[16],
    const REAL v[3]);

template <typename T0, typename T1, typename T2>
void Vec3_Mat4Vec3_AffineProjection(
    T0 y0[3],
    const T1 a[16],
    const T2 x0[3]);

template <typename T0, typename T1, typename T2>
void Vec3_Vec3Mat4_AffineProjection(
    T0 y0[3],
    const T1 x0[3],
    const T2 a[16]);

template <typename T>
void Vec3_Mat4Vec3_Affine(
    T y0[3],
    const T a[16],
    const T x0[3]);

// ------------------------------------
// below: function with mat4 and quarternion

template <typename REAL>
DFM2_INLINE void Mat4_Quat(
    REAL r[], const REAL q[]);

DFM2_INLINE void Mat4_QuatConj(
    double r[], const double q[]);

template <typename REAL>
DFM2_INLINE void Mat4_AffineTransQuat(
    REAL r[],
    const REAL q[]);

DFM2_INLINE void Mat4_ScaleRotTrans(
    double m[16],
    double scale,
    const double quat[4],
    const double trans[3]);

// ---------------------------------------

template<typename T>
class CMat4;

template <typename T>
CMat4<T> operator * (const CMat4<T>& lhs, const CMat4<T>& rhs);

template <typename T>
CMat4<T> operator - (const CMat4<T>& lhs, const CMat4<T>& rhs);

template <typename T>
CMat4<T> operator + (const CMat4<T>& lhs, const CMat4<T>& rhs);

/**
 * @brief 4x4 matrix class
 * @class 4x4 matrix class
 * @tparam REAL value type of the matrix. defiend for "double" and "float" for the static library.
 */
template <typename REAL>
class CMat4 {
public:
  CMat4 (){};
  CMat4 (const float* pm){ for(int i=0;i<16;++i){ mat[i] = (REAL)pm[i]; } }
  CMat4 (const double* pm){ for(int i=0;i<16;++i){ mat[i] = (REAL)pm[i]; } }
public:
  // ------------------------
  // below: "set" functions
  void Set_AffineTranslate(REAL x, REAL y, REAL z){
    Mat4_AffineTranslation(mat,
                           x, y, z);
  }
  void Set_Quaternion(const REAL* q){
    Mat4_Quat(mat,
              q);
  }
  void SetZero() {
    for(auto& v : mat){ v = 0; }
  }
  void SetIdentity() {
    for(auto& v : mat){ v = 0; }
    mat[0*4+0] = 1;
    mat[1*4+1] = 1;
    mat[2*4+2] = 1;
    mat[3*4+3] = 1;
  }
  void SetScale(REAL x, REAL y, REAL z){
    this->SetZero();
    mat[0*4+0] = x;
    mat[1*4+1] = y;
    mat[2*4+2] = z;
    mat[3*4+3] = 1;
  }
  // -----------------------
  CMat4<double> Double() const {
    return CMat4<double>(mat);
  }
  CMat4<REAL> MatMat(const CMat4<REAL>& mat0) const;
  CMat4<REAL> Transpose() const{
    CMat4<REAL> m1;
    for(int i=0;i<4;++i){
      for(int j=0;j<4;++j){
        m1.mat[i*4+j] = mat[j*4+i];
      }
    }
    return m1;
  }
  CMat4<REAL> Inverse() const;
  // ----------------------
  // below: static function
  static CMat4<REAL> Identity(){
    CMat4<REAL> m;
    m.SetIdentity();
    return m;
  }
  static CMat4<REAL> Scale(REAL s){
    CMat4<REAL> m;
    m.SetZero();
    m.mat[0*4+0] = s;
    m.mat[1*4+1] = s;
    m.mat[2*4+2] = s;
    m.mat[3*4+3] = 1.0;
    return m;
  }
  static CMat4<REAL> Spin(const REAL* v){
    CMat4<REAL> m;
    m.SetZero();
    m.mat[0*4+0] =  0;     m.mat[0*4+1] = -v[2];   m.mat[0*4+2] = +v[1];
    m.mat[1*4+0] = +v[2];  m.mat[1*4+1] = 0;       m.mat[1*4+2] = -v[0];
    m.mat[2*4+0] = -v[1];  m.mat[2*4+1] = +v[0];   m.mat[2*4+2] = 0;
    m.mat[3*4+3] = 1.0;
    return m;
  }
  static CMat4<REAL> Quat(const REAL* q);
  static CMat4<REAL> Translate(const REAL* v){
    CMat4<REAL> m;
    m.SetIdentity();
    m.mat[0*4+3] = v[0];
    m.mat[1*4+3] = v[1];
    m.mat[2*4+3] = v[2];
    return m;
  }
  static CMat4<REAL> Mat3(const REAL* v){
    CMat4<REAL> m;
    m.SetZero();
    m.mat[0*4+0] = v[0*3+0];  m.mat[0*4+1] = v[0*3+1];  m.mat[0*4+2] = v[0*3+2];
    m.mat[1*4+0] = v[1*3+0];  m.mat[1*4+1] = v[1*3+1];  m.mat[1*4+2] = v[1*3+2];
    m.mat[2*4+0] = v[2*3+0];  m.mat[2*4+1] = v[2*3+1];  m.mat[2*4+2] = v[2*3+2];
    m.mat[3*4+3] = 1.0;
    return m;
  }

public:
  REAL mat[16];
};
using CMat4d = CMat4<double>;
using CMat4f = CMat4<float>;

}

#ifdef DFM2_HEADER_ONLY
#  include "delfem2/mat4.cpp"
#endif

#endif
