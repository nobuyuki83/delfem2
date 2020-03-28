/*
 * Copyright (c) 2019 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */

/**
 * @file 3x3 matrix class (CMat3) and functions
 */


#ifndef DFM2_MAT4_H
#define DFM2_MAT4_H

#include <vector>
#include <cassert>
#include <cmath>
#include <iostream>
#include <limits> // using NaN Check

#define NEARLY_ZERO 1.e-16

// -----------------------------

namespace delfem2 {


void Mat4Vec3(double vo[3],
              const double M[16],
              const double vi[3]);

// -------------------------------

template <typename T>
void MatVec4(T v[4],
    const T A[16],
    const T x[4]);

template <typename T>
void MatMat4(
    T* C,
    const T* A,
    const T* B);


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

/**
 * @func multiply translation affine matrix from left to an affine matrix in 3D
 */
template <typename REAL>
void Translate_Mat4Affine(
    REAL A[16],
    const REAL v[3]);

template <typename T>
void Vec3_Mat4Vec3_AffineProjection(T y0[3],
                                    const T a[16],
                                    const T x0[3]);

template <typename T>
void Vec3_Mat4Vec3_Affine(T y0[3],
                          const T a[16],
                          const T x0[3]);

template <typename REAL>
void Copy_Mat4(REAL M0[16],
               const REAL M1[16])
{
  for(int i=0;i<16;++i){ M0[i] = M1[i]; }
}

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
 * @details Row major data storage. The template is defined for "float" and "double"
 */
template <typename REAL>
class CMat4 {
public:
  CMat4 (){
  };
public:
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
  CMat4<REAL> MatMat(const CMat4<REAL>& mat0) const;
  // -------------
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

#endif
