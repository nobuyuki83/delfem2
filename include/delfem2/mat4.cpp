/*
 * Copyright (c) 2019 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */

#include <cstring>
#include "delfem2/mat4.h"

// ------------------------

namespace delfem2 {
namespace mat4 {

template <typename REAL>
DFM2_INLINE void CalcInvMat
 (REAL *a,
  const unsigned int n,
  int &info)
{
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
    tmp1 = 1.0 / a[i * n + i];
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

}
}


// ------------------------

DFM2_INLINE void delfem2::Mat4Vec3(
    double vo[3],
    const double M[16],
    const double vi[3])
{
  vo[0] = M[0*4+0]*vi[0] + M[0*4+1]*vi[1] + M[0*4+2]*vi[2];
  vo[1] = M[1*4+0]*vi[0] + M[1*4+1]*vi[1] + M[1*4+2]*vi[2];
  vo[2] = M[2*4+0]*vi[0] + M[2*4+1]*vi[1] + M[2*4+2]*vi[2];
}

template <typename T>
DFM2_INLINE void delfem2::MatVec4
(T v[4],
 const T A[16],
 const T x[4])
{
  v[0] = A[0*4+0]*x[0] + A[0*4+1]*x[1] + A[0*4+2]*x[2] + A[0*4+3]*x[3];
  v[1] = A[1*4+0]*x[0] + A[1*4+1]*x[1] + A[1*4+2]*x[2] + A[1*4+3]*x[3];
  v[2] = A[2*4+0]*x[0] + A[2*4+1]*x[1] + A[2*4+2]*x[2] + A[2*4+3]*x[3];
  v[3] = A[3*4+0]*x[0] + A[3*4+1]*x[1] + A[3*4+2]*x[2] + A[3*4+3]*x[3];
}
#ifndef DFM2_HEADER_ONLY
template void delfem2::MatVec4(float v[4], const float A[16], const float x[4]);
template void delfem2::MatVec4(double v[4], const double A[16], const double x[4]);
#endif

template <typename T>
DFM2_INLINE void delfem2::VecMat4(
    T v[4],
    const T x[4],
    const T A[16])
{
  v[0] = A[0*4+0]*x[0] + A[1*4+0]*x[1] + A[2*4+0]*x[2] + A[3*4+0]*x[3];
  v[1] = A[0*4+1]*x[0] + A[1*4+1]*x[1] + A[2*4+1]*x[2] + A[3*4+1]*x[3];
  v[2] = A[0*4+2]*x[0] + A[1*4+2]*x[1] + A[2*4+2]*x[2] + A[3*4+2]*x[3];
  v[3] = A[0*4+3]*x[0] + A[1*4+3]*x[1] + A[2*4+3]*x[2] + A[3*4+3]*x[3];
}
#ifndef DFM2_HEADER_ONLY
template void delfem2::VecMat4(float v[4], const float x[4], const float A[16]);
template void delfem2::VecMat4(double v[4], const double x[4], const double A[16]);
#endif
  
// --------------------------

template <typename T>
DFM2_INLINE void delfem2::MatMat4
(T* C,
 const T* A, const T* B)
{
  for(int i=0;i<4;i++){
    for(int j=0;j<4;j++){
      C[i*4+j] = A[i*4+0]*B[0*4+j] + A[i*4+1]*B[1*4+j] + A[i*4+2]*B[2*4+j] + A[i*4+3]*B[3*4+j];
    }
  }
}
#ifndef DFM2_HEADER_ONLY
template void delfem2::MatMat4(float* C, const float* A, const float* B);
template void delfem2::MatMat4(double* C, const double* A, const double* B);
#endif

  
// ---------------------------

template <typename T>
DFM2_INLINE void delfem2::Vec3_Mat4Vec3_AffineProjection
(T y0[3],
 const T a[16],
 const T x0[3])
{
  const T x1[4] = {x0[0], x0[1], x0[2], 1.0};
  T y1[4]; MatVec4(y1,a,x1);
  y0[0] = y1[0]/y1[3];
  y0[1] = y1[1]/y1[3];
  y0[2] = y1[2]/y1[3];
}
#ifndef DFM2_HEADER_ONLY
template void delfem2::Vec3_Mat4Vec3_AffineProjection(float y0[3], const float a[16], const float x0[3]);
template void delfem2::Vec3_Mat4Vec3_AffineProjection(double y0[3], const double a[16], const double x0[3]);
#endif

// ----------------------

template <typename T>
DFM2_INLINE void delfem2::Vec3_Mat4Vec3_Affine
 (T y0[3],
  const T a[16],
  const T x0[3])
{
  const T x1[4] = {x0[0], x0[1], x0[2], 1.0};
  T y1[4]; MatVec4(y1,a,x1);
  y0[0] = y1[0];
  y0[1] = y1[1];
  y0[2] = y1[2];
}
#ifndef DFM2_HEADER_ONLY
template void delfem2::Vec3_Mat4Vec3_Affine(float y0[3], const float a[16], const float x0[3]);
template void delfem2::Vec3_Mat4Vec3_Affine(double y0[3], const double a[16], const double x0[3]);
#endif

  
// ----------------------

template <typename T>
DFM2_INLINE void delfem2::Mat4_AffineScale
(T A[16],
 T s)
{
  for(int i=0;i<16;++i){ A[i] = 0.0; }
  A[0*4+0] = s;
  A[1*4+1] = s;
  A[2*4+2] = s;
  A[3*4+3] = 1.0;
}
#ifndef DFM2_HEADER_ONLY
template void delfem2::Mat4_AffineScale(float A[16], float s);
template void delfem2::Mat4_AffineScale(double A[16], double s);
#endif
  
// ------------------------

template <typename T>
DFM2_INLINE void delfem2::Mat4_AffineTranslation
(T A[16],
 T dx, T dy, T dz)
{
  for(auto i=0;i<16;++i){ A[i] = 0.0; }
  for(int i=0;i<4;++i){ A[i*4+i] = 1.0; }
  A[0*4+3] = dx;
  A[1*4+3] = dy;
  A[2*4+3] = dz;
}
#ifndef DFM2_HEADER_ONLY
template void delfem2::Mat4_AffineTranslation(
    float A[16],
    float dx, float dy, float dz);
template void delfem2::Mat4_AffineTranslation(
    double A[16],
    double dx, double dy, double dz);
#endif

template <typename REAL>
DFM2_INLINE void delfem2::Mat4_AffineTransTranslate(
    REAL r[],
    const REAL t[])
{
  // column 0
  r[ 0] = 1;
  r[ 1] = 0;
  r[ 2] = 0;
  r[ 3] = 0;
  // column 1
  r[ 4] = 0;
  r[ 5] = 1;
  r[ 6] = 0;
  r[ 7] = 0;
  // column 2
  r[ 8] = 0;
  r[ 9] = 0;
  r[10] = 1;
  r[11] = 0;
  // column 3
  r[12] = t[0];
  r[13] = t[1];
  r[14] = t[2];
  r[15] = 1;
}
#ifndef DFM2_HEADER_ONLY
template void delfem2::Mat4_AffineTransTranslate(float r[], const float t[]);
template void delfem2::Mat4_AffineTransTranslate(double r[], const double t[]);
#endif


// --------------------------

template <typename T>
DFM2_INLINE void delfem2::Mat4_AffineRotationRodriguez
(T A[16],
 T dx, T dy, T dz)
{
  for(int i=0;i<16;++i){ A[i] = 0.0; }
  //
  const double sqlen = dx*dx+dy*dy+dz*dz;
  const double tmp1 = 1.0/(1+0.25*sqlen);
  A[0*4+0] = 1+tmp1*(+0.5*dx*dx-0.5*sqlen);
  A[0*4+1] =  +tmp1*(-dz+0.5*dx*dy);
  A[0*4+2] =  +tmp1*(+dy+0.5*dx*dz);
  A[0*4+3] = 0.0;
  //
  A[1*4+0] =  +tmp1*(+dz+0.5*dy*dx);
  A[1*4+1] = 1+tmp1*(+0.5*dy*dy-0.5*sqlen);
  A[1*4+2] =  +tmp1*(-dx+0.5*dy*dz);
  A[1*4+3] = 0.0;
  //
  A[2*4+0] =  +tmp1*(-dy+0.5*dz*dx);
  A[2*4+1] =  +tmp1*(+dx+0.5*dz*dy);
  A[2*4+2] = 1+tmp1*(+0.5*dz*dz-0.5*sqlen);
  A[2*4+3] = 0.0;
  //
  A[3*4+0] = 0.0;
  A[3*4+1] = 0.0;
  A[3*4+2] = 0.0;
  A[3*4+3] = 1.0;
}
#ifndef DFM2_HEADER_ONLY
template void delfem2::Mat4_AffineRotationRodriguez(
    float A[16],
    float dx, float dy, float dz);
template void delfem2::Mat4_AffineRotationRodriguez(
    double A[16],
    double dx, double dy, double dz);
#endif

// ------------------------------------------------

template <typename REAL>
void delfem2::Mat4_Identity(
    REAL A[16])
{
  for(int i=0;i<16;++i){ A[i] = 0; }
  A[0*4+0] = 1;
  A[1*4+1] = 1;
  A[2*4+2] = 1;
  A[3*4+3] = 1;
}
#ifndef DFM2_HEADER_ONLY
template void delfem2::Mat4_Identity(float A[16]);
template void delfem2::Mat4_Identity(double A[16]);
#endif

  
// ------------------------------------------------
  
template <typename REAL>
void delfem2::Rotate_Mat4AffineRodriguez(
    REAL A[16],
    const REAL V[3])
{
  REAL B[16];
  Mat4_AffineRotationRodriguez(B,
                               V[0],V[1],V[2]);
  REAL C[16];
  MatMat4(C,
      B,A);
  
  for(int i=0;i<16;++i){ A[i] = C[i]; }
}
#ifndef DFM2_HEADER_ONLY
template void delfem2::Rotate_Mat4AffineRodriguez(float A[16], const float V[3]);
template void delfem2::Rotate_Mat4AffineRodriguez(double A[16], const double V[3]);
#endif
  
  
template <typename REAL>
void delfem2::Translate_Mat4Affine(
    REAL A[16],
    const REAL V[3])
{
  A[0*4+3] += V[0];
  A[1*4+3] += V[1];
  A[2*4+3] += V[2];
}
#ifndef DFM2_HEADER_ONLY
template void delfem2::Translate_Mat4Affine(float A[16], const float V[3]);
template void delfem2::Translate_Mat4Affine(double A[16], const double V[3]);
#endif


DFM2_INLINE void delfem2::Mat4_ScaleRotTrans(
     double m[16],
     double scale,
     const double quat[4],
     const double trans[3])
{
  delfem2::Mat4_Quat(m, quat);
  for(int i=0;i<3;++i){
    for(int j=0;j<3;++j){
      m[i*4+j] *= scale;
    }
  }
  m[0*4+3] = trans[0];
  m[1*4+3] = trans[1];
  m[2*4+3] = trans[2];
}

template <typename REAL>
DFM2_INLINE void delfem2::Mat4_Quat(
    REAL r[],
    const REAL q[])
{
  double x2 = q[1] * q[1] * 2.0;
  double y2 = q[2] * q[2] * 2.0;
  double z2 = q[3] * q[3] * 2.0;
  double xy = q[1] * q[2] * 2.0;
  double yz = q[2] * q[3] * 2.0;
  double zx = q[3] * q[1] * 2.0;
  double xw = q[1] * q[0] * 2.0;
  double yw = q[2] * q[0] * 2.0;
  double zw = q[3] * q[0] * 2.0;
  r[ 0] = 1.0 - y2 - z2;
  r[ 1] = xy - zw;
  r[ 2] = zx + yw;
  r[ 4] = xy + zw;
  r[ 5] = 1.0 - z2 - x2;
  r[ 6] = yz - xw;
  r[ 8] = zx - yw;
  r[ 9] = yz + xw;
  r[10] = 1.0 - x2 - y2;
  r[ 3] = r[ 7] = r[11] = r[12] = r[13] = r[14] = 0.0;
  r[15] = 1.0;
}
#ifndef DFM2_HEADER_ONLY
template void delfem2::Mat4_Quat(float r[], const float q[]);
template void delfem2::Mat4_Quat(double r[], const double q[]);
#endif

// return transpose matrix of Mat4_Quat
DFM2_INLINE void delfem2::Mat4_QuatConj(
    double r[],
    const double q[])
{
  double x2 = q[1] * q[1] * 2.0;
  double y2 = q[2] * q[2] * 2.0;
  double z2 = q[3] * q[3] * 2.0;
  double xy = q[1] * q[2] * 2.0;
  double yz = q[2] * q[3] * 2.0;
  double zx = q[3] * q[1] * 2.0;
  double xw = q[1] * q[0] * 2.0;
  double yw = q[2] * q[0] * 2.0;
  double zw = q[3] * q[0] * 2.0;
  r[ 0] = 1.0 - y2 - z2;
  r[ 1] = xy + zw;
  r[ 2] = zx - yw;
  r[ 4] = xy - zw;
  r[ 5] = 1.0 - z2 - x2;
  r[ 6] = yz + xw;
  r[ 8] = zx + yw;
  r[ 9] = yz - xw;
  r[10] = 1.0 - x2 - y2;
  r[ 3] = r[ 7] = r[11] = r[12] = r[13] = r[14] = 0.0;
  r[15] = 1.0;
}

template <typename REAL>
DFM2_INLINE void delfem2::Mat4_AffineTransQuat(
    REAL r[],
    const REAL q[])
{
  REAL x2 = q[1] * q[1] * 2;
  REAL y2 = q[2] * q[2] * 2;
  REAL z2 = q[3] * q[3] * 2;
  REAL xy = q[1] * q[2] * 2;
  REAL yz = q[2] * q[3] * 2;
  REAL zx = q[3] * q[1] * 2;
  REAL xw = q[1] * q[0] * 2;
  REAL yw = q[2] * q[0] * 2;
  REAL zw = q[3] * q[0] * 2;
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
#ifndef DFM2_HEADER_ONLY
template void delfem2::Mat4_AffineTransQuat(float r[], const float q[]);
template void delfem2::Mat4_AffineTransQuat(double r[], const double q[]);
#endif


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


DFM2_INLINE void delfem2::Inverse_Mat4(
    double minv[16],
    const double m[16])
{
  for(int i=0;i<16;++i){ minv[i] = m[i]; }
  int info; mat4::CalcInvMat(minv,4,info);
}

// ------------------------------------------------------------------

template <typename T>
delfem2::CMat4<T> delfem2::CMat4<T>::MatMat(const CMat4<T>& mat0) const{
  CMat4 m;
  ::delfem2::MatMat4(m.mat,
                  this->mat, mat0.mat);
  return m;
}
#ifndef DFM2_HEADER_ONLY
template delfem2::CMat4<float> delfem2::CMat4<float>::MatMat(const CMat4<float>& mat0) const;
template delfem2::CMat4<double> delfem2::CMat4<double>::MatMat(const CMat4<double>& mat0) const;
#endif



template <typename REAL>
delfem2::CMat4<REAL> delfem2::CMat4<REAL>::Quat(const REAL* q)
{
  CMat4<REAL> m;
  REAL x2 = q[1] * q[1] * 2.0;
  REAL y2 = q[2] * q[2] * 2.0;
  REAL z2 = q[3] * q[3] * 2.0;
  REAL xy = q[1] * q[2] * 2.0;
  REAL yz = q[2] * q[3] * 2.0;
  REAL zx = q[3] * q[1] * 2.0;
  REAL xw = q[1] * q[0] * 2.0;
  REAL yw = q[2] * q[0] * 2.0;
  REAL zw = q[3] * q[0] * 2.0;
  m.SetZero();
  m.mat[0*4+0] = 1.0 - y2 - z2; m.mat[0*4+1] = xy - zw;         m.mat[0*4+2] = zx + yw;
  m.mat[1*4+0] = xy + zw;       m.mat[1*4+1] = 1.0 - z2 - x2;   m.mat[1*4+2] = yz - xw;
  m.mat[2*4+0] = zx - yw;       m.mat[2*4+1] = yz + xw;         m.mat[2*4+2] = 1.0 - x2 - y2;
  m.mat[3*4+3] = 1.0;
  return m;
}
#ifndef DFM2_HEADER_ONLY
template delfem2::CMat4<float> delfem2::CMat4<float>::Quat(const float* q);
template delfem2::CMat4<double> delfem2::CMat4<double>::Quat(const double* q);
#endif

// ---------------------------

namespace delfem2{

template <typename T>
CMat4<T> operator * (const CMat4<T>& lhs, const CMat4<T>& rhs)
{
  CMat4<T> q;
  MatMat4(q.mat, lhs.mat, rhs.mat);
  return q;
}
#ifndef DFM2_HEADER_ONLY
template CMat4d operator * (const CMat4d& lhs, const CMat4d& rhs);
template CMat4f operator * (const CMat4f& lhs, const CMat4f& rhs);
#endif


template <typename T>
CMat4<T> operator - (const CMat4<T>& lhs, const CMat4<T>& rhs)
{
  CMat4<T> q;
  for(int i=0;i<16;++i){ q.mat[i] = lhs.mat[i] - rhs.mat[i]; }
  return q;
}
#ifndef DFM2_HEADER_ONLY
template CMat4d operator - (const CMat4d& lhs, const CMat4d& rhs);
template CMat4f operator - (const CMat4f& lhs, const CMat4f& rhs);
#endif


template <typename T>
CMat4<T> operator + (const CMat4<T>& lhs, const CMat4<T>& rhs)
{
  CMat4<T> q;
  for(int i=0;i<16;++i){ q.mat[i] = lhs.mat[i] + rhs.mat[i]; }
  return q;
}
#ifndef DFM2_HEADER_ONLY
template CMat4d operator + (const CMat4d& lhs, const CMat4d& rhs);
template CMat4f operator + (const CMat4f& lhs, const CMat4f& rhs);
#endif

}

// --------------------------------------------

template <typename REAL>
delfem2::CMat4<REAL> delfem2::CMat4<REAL>::Inverse() const
{
  CMat4<REAL> m;
  std::memcpy(m.mat, mat, sizeof(REAL)*16);
  int info;
  mat4::CalcInvMat(m.mat, 4, info);
  return m;
}
#ifndef DFM2_HEADER_ONLY
template delfem2::CMat4d delfem2::CMat4d::Inverse() const;
template delfem2::CMat4f delfem2::CMat4f::Inverse() const;
#endif
