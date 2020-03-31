/*
 * Copyright (c) 2019 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */

#include <random>
#include "delfem2/mat4.h"

namespace dfm2 = delfem2;

// ------------------------

void dfm2::Mat4Vec3(
    double vo[3],
    const double M[16],
    const double vi[3])
{
  vo[0] = M[0*4+0]*vi[0] + M[0*4+1]*vi[1] + M[0*4+2]*vi[2];
  vo[1] = M[1*4+0]*vi[0] + M[1*4+1]*vi[1] + M[1*4+2]*vi[2];
  vo[2] = M[2*4+0]*vi[0] + M[2*4+1]*vi[1] + M[2*4+2]*vi[2];
}

template <typename T>
void dfm2::MatVec4
(T v[4],
 const T A[16],
 const T x[4])
{
  v[0] = A[0*4+0]*x[0] + A[0*4+1]*x[1] + A[0*4+2]*x[2] + A[0*4+3]*x[3];
  v[1] = A[1*4+0]*x[0] + A[1*4+1]*x[1] + A[1*4+2]*x[2] + A[1*4+3]*x[3];
  v[2] = A[2*4+0]*x[0] + A[2*4+1]*x[1] + A[2*4+2]*x[2] + A[2*4+3]*x[3];
  v[3] = A[3*4+0]*x[0] + A[3*4+1]*x[1] + A[3*4+2]*x[2] + A[3*4+3]*x[3];
}
template void dfm2::MatVec4(float v[4], const float A[16], const float x[4]);
template void dfm2::MatVec4(double v[4], const double A[16], const double x[4]);
  

template <typename T>
void dfm2::MatMat4
(T* C,
 const T* A, const T* B)
{
  for(int i=0;i<4;i++){
    for(int j=0;j<4;j++){
      C[i*4+j] = A[i*4+0]*B[0*4+j] + A[i*4+1]*B[1*4+j] + A[i*4+2]*B[2*4+j] + A[i*4+3]*B[3*4+j];
    }
  }
}
template void dfm2::MatMat4(float* C, const float* A, const float* B);
template void dfm2::MatMat4(double* C, const double* A, const double* B);

  
// ---------------------------

template <typename T>
void dfm2::Vec3_Mat4Vec3_AffineProjection
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
template void dfm2::Vec3_Mat4Vec3_AffineProjection(float y0[3], const float a[16], const float x0[3]);
template void dfm2::Vec3_Mat4Vec3_AffineProjection(double y0[3], const double a[16], const double x0[3]);

// ----------------------

template <typename T>
void dfm2::Vec3_Mat4Vec3_Affine
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
template void dfm2::Vec3_Mat4Vec3_Affine(float y0[3], const float a[16], const float x0[3]);
template void dfm2::Vec3_Mat4Vec3_Affine(double y0[3], const double a[16], const double x0[3]);

  
// ----------------------

template <typename T>
void dfm2::Mat4_AffineScale
(T A[16],
 T s)
{
  for(int i=0;i<16;++i){ A[i] = 0.0; }
  A[0*4+0] = s;
  A[1*4+1] = s;
  A[2*4+2] = s;
  A[3*4+3] = 1.0;
}
template void dfm2::Mat4_AffineScale(float A[16], float s);
template void dfm2::Mat4_AffineScale(double A[16], double s);
  
// ------------------------

template <typename T>
void dfm2::Mat4_AffineTranslation
(T A[16],
 T dx, T dy, T dz)
{
  for(auto i=0;i<16;++i){ A[i] = 0.0; }
  for(int i=0;i<4;++i){ A[i*4+i] = 1.0; }
  A[0*4+3] = dx;
  A[1*4+3] = dy;
  A[2*4+3] = dz;
}
template void dfm2::Mat4_AffineTranslation(float A[16],  float dx, float dy, float dz);
template void dfm2::Mat4_AffineTranslation(double A[16],  double dx, double dy, double dz);
  
// --------------------------

template <typename T>
void dfm2::Mat4_AffineRotationRodriguez
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
template void dfm2::Mat4_AffineRotationRodriguez(float A[16],
                                                 float dx, float dy, float dz);
template void dfm2::Mat4_AffineRotationRodriguez(double A[16],
                                                 double dx, double dy, double dz);

// ------------------------------------------------

template <typename REAL>
void dfm2::Mat4_Identity(
    REAL A[16])
{
  for(int i=0;i<16;++i){ A[i] = 0; }
  A[0*4+0] = 1;
  A[1*4+1] = 1;
  A[2*4+2] = 1;
  A[3*4+3] = 1;
}
template void dfm2::Mat4_Identity(float A[16]);
template void dfm2::Mat4_Identity(double A[16]);

  
// ------------------------------------------------
  
template <typename REAL>
void dfm2::Rotate_Mat4AffineRodriguez(
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
template void dfm2::Rotate_Mat4AffineRodriguez(float A[16], const float V[3]);
template void dfm2::Rotate_Mat4AffineRodriguez(double A[16], const double V[3]);
  
  
template <typename REAL>
void dfm2::Translate_Mat4Affine(
    REAL A[16],
    const REAL V[3])
{
  A[0*4+3] += V[0];
  A[1*4+3] += V[1];
  A[2*4+3] += V[2];
}
template void dfm2::Translate_Mat4Affine(float A[16], const float V[3]);
template void dfm2::Translate_Mat4Affine(double A[16], const double V[3]);


void dfm2::Mat4_ScaleRotTrans
    (double m[16],
     double scale,
     const double quat[4],
     const double trans[3])
{
  dfm2::Mat4_Quat(m, quat);
  for(int i=0;i<3;++i){
    for(int j=0;j<3;++j){
      m[i*4+j] *= scale;
    }
  }
  m[0*4+3] = trans[0];
  m[1*4+3] = trans[1];
  m[2*4+3] = trans[2];
}

void dfm2::Mat4_Quat
    (double r[], const double q[])
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

// return transpose matrix of Mat4_Quat
void dfm2::Mat4_QuatConj(
    double r[], const double q[])
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

/*
void dfm2::MatMat4(
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

void dfm2::Copy_Mat4(double m1[16], const double m0[16])
{
  for(int i=0;i<16;++i){ m1[i] = m0[i]; }
}
 */


// ------------------------------------------------------------------

template <typename T>
dfm2::CMat4<T> dfm2::CMat4<T>::MatMat(const CMat4<T>& mat0) const{
  CMat4 m;
  ::dfm2::MatMat4(m.mat,
                  this->mat, mat0.mat);
  return m;
}
template dfm2::CMat4<double> dfm2::CMat4<double>::MatMat(const CMat4<double>& mat0) const;



template <typename REAL>
dfm2::CMat4<REAL> dfm2::CMat4<REAL>::Quat(const REAL* q)
{
  CMat4<REAL> m;
  double x2 = q[1] * q[1] * 2.0;
  double y2 = q[2] * q[2] * 2.0;
  double z2 = q[3] * q[3] * 2.0;
  double xy = q[1] * q[2] * 2.0;
  double yz = q[2] * q[3] * 2.0;
  double zx = q[3] * q[1] * 2.0;
  double xw = q[1] * q[0] * 2.0;
  double yw = q[2] * q[0] * 2.0;
  double zw = q[3] * q[0] * 2.0;
  m.SetZero();
  m.mat[0*4+0] = 1.0 - y2 - z2; m.mat[0*4+1] = xy - zw;         m.mat[0*4+2] = zx + yw;
  m.mat[1*4+0] = xy + zw;       m.mat[1*4+1] = 1.0 - z2 - x2;   m.mat[1*4+2] = yz - xw;
  m.mat[2*4+0] = zx - yw;       m.mat[2*4+1] = yz + xw;         m.mat[2*4+2] = 1.0 - x2 - y2;
  m.mat[3*4+3] = 1.0;
  return m;
}
template dfm2::CMat4<double> dfm2::CMat4<double>::Quat(const double* q);

// ---------------------------

namespace delfem2{

template <typename T>
CMat4<T> operator * (const CMat4<T>& lhs, const CMat4<T>& rhs)
{
  CMat4<T> q;
  MatMat4(q.mat, lhs.mat, rhs.mat);
  return q;
}
template CMat4d operator * (const CMat4d& lhs, const CMat4d& rhs);
template CMat4f operator * (const CMat4f& lhs, const CMat4f& rhs);


template <typename T>
CMat4<T> operator - (const CMat4<T>& lhs, const CMat4<T>& rhs)
{
  CMat4<T> q;
  for(int i=0;i<16;++i){ q.mat[i] = lhs.mat[i] - rhs.mat[i]; }
  return q;
}
template CMat4d operator - (const CMat4d& lhs, const CMat4d& rhs);
template CMat4f operator - (const CMat4f& lhs, const CMat4f& rhs);


template <typename T>
CMat4<T> operator + (const CMat4<T>& lhs, const CMat4<T>& rhs)
{
  CMat4<T> q;
  for(int i=0;i<16;++i){ q.mat[i] = lhs.mat[i] + rhs.mat[i]; }
  return q;
}
template CMat4d operator + (const CMat4d& lhs, const CMat4d& rhs);
template CMat4f operator + (const CMat4f& lhs, const CMat4f& rhs);

}
