/*
 * Copyright (c) 2019 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */


#include <cmath>
#include <cstring>
#include <iostream>

#include "delfem2/camera.h"

// -----------------------------------------------------------------------

namespace delfem2 {
namespace camera {

//! @brief multiply two quaternion
DFM2_INLINE void QuatQuat(double r[], const double p[], const double q[])
{
  r[0] = p[0] * q[0] - p[1] * q[1] - p[2] * q[2] - p[3] * q[3];
  r[1] = p[0] * q[1] + p[1] * q[0] + p[2] * q[3] - p[3] * q[2];
  r[2] = p[0] * q[2] - p[1] * q[3] + p[2] * q[0] + p[3] * q[1];
  r[3] = p[0] * q[3] + p[1] * q[2] - p[2] * q[1] + p[3] * q[0];
}

//! @brief transform vector with quaternion
DFM2_INLINE void QuatVec(double vo[], const double q[], const double vi[])
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
  
  vo[0] = (1.0 - y2 - z2)*vi[0] + (xy + zw      )*vi[1] + (zx - yw      )*vi[2];
  vo[1] = (xy - zw      )*vi[0] + (1.0 - z2 - x2)*vi[1] + (yz + xw      )*vi[2];
  vo[2] = (zx + yw      )*vi[0] + (yz - xw      )*vi[1] + (1.0 - x2 - y2)*vi[2];
}


//! @brief transform vector with conjugate of quaternion
DFM2_INLINE void QuatConjVec(double vo[], const double q[], const double vi[])
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
  
  vo[0] = (1.0 - y2 - z2)*vi[0] + (xy - zw      )*vi[1] + (zx + yw      )*vi[2];
  vo[1] = (xy + zw      )*vi[0] + (1.0 - z2 - x2)*vi[1] + (yz - xw      )*vi[2];
  vo[2] = (zx - yw      )*vi[0] + (yz + xw      )*vi[1] + (1.0 - x2 - y2)*vi[2];
}



//! @brief copy quaternion
DFM2_INLINE void CopyQuat
 (double r[],
  const double p[])
{
  r[0] = p[0];
  r[1] = p[1];
  r[2] = p[2];
  r[3] = p[3];
}

template <typename REAL>
DFM2_INLINE void Mat4_AffineTransQuat(
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

template <typename REAL>
DFM2_INLINE void Mat4_AffineTransTranslate(
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

template <typename REAL>
DFM2_INLINE void Mat4_Identity(REAL r[])
{
  r[ 0] = 1;  r[ 1] = 0;  r[ 2] = 0;  r[ 3] = 0;
  r[ 4] = 0;  r[ 5] = 1;  r[ 6] = 0;  r[ 7] = 0;
  r[ 8] = 0;  r[ 9] = 0;  r[10] = 1;  r[11] = 0;
  r[12] = 0;  r[13] = 0;  r[14] = 0;  r[15] = 1;
}


DFM2_INLINE void Normalize3D(float vec[3])
{
  float len = sqrt(vec[0]*vec[0]+vec[1]*vec[1]+vec[2]*vec[2]);
  float leninv = 1.0/len;
  vec[0] *= leninv;
  vec[1] *= leninv;
  vec[2] *= leninv;
}

DFM2_INLINE void Cross3D(float r[3], const float v1[3], const float v2[3]){
  r[0] = v1[1]*v2[2] - v2[1]*v1[2];
  r[1] = v1[2]*v2[0] - v2[2]*v1[0];
  r[2] = v1[0]*v2[1] - v2[0]*v1[1];
}

DFM2_INLINE void InverseMat3
(float Ainv[],
 const float A[])
{
  const float det =
  + A[0]*A[4]*A[8] + A[3]*A[7]*A[2] + A[6]*A[1]*A[5]
  - A[0]*A[7]*A[5] - A[6]*A[4]*A[2] - A[3]*A[1]*A[8];
  const float inv_det = 1.f/det;
  Ainv[0] = inv_det*(A[4]*A[8]-A[5]*A[7]);
  Ainv[1] = inv_det*(A[2]*A[7]-A[1]*A[8]);
  Ainv[2] = inv_det*(A[1]*A[5]-A[2]*A[4]);
  Ainv[3] = inv_det*(A[5]*A[6]-A[3]*A[8]);
  Ainv[4] = inv_det*(A[0]*A[8]-A[2]*A[6]);
  Ainv[5] = inv_det*(A[2]*A[3]-A[0]*A[5]);
  Ainv[6] = inv_det*(A[3]*A[7]-A[4]*A[6]);
  Ainv[7] = inv_det*(A[1]*A[6]-A[0]*A[7]);
  Ainv[8] = inv_det*(A[0]*A[4]-A[1]*A[3]);
}

DFM2_INLINE void Mat3Vec
(float vo[3],
 const float mat[9],
 const float vi[3])
{
  vo[0] = mat[0]*vi[0] + mat[1]*vi[1] + mat[2]*vi[2];
  vo[1] = mat[3]*vi[0] + mat[4]*vi[1] + mat[5]*vi[2];
  vo[2] = mat[6]*vi[0] + mat[7]*vi[1] + mat[8]*vi[2];
}

/**
 * @brief affine matrix
 * @details column major order
 */
DFM2_INLINE void Mat4_AffineTransProjectionOrtho(
    float mP[16],
    double l, double r,
    double b, double t,
    double n, double f)
{
  // column 0
  mP[0*4+0] = 2.0/(r-l);
  mP[0*4+1] = 0.0;
  mP[0*4+2] = 0.0;
  mP[0*4+3] = 0.0;
  // column 1
  mP[1*4+0] = 0.0;
  mP[1*4+1] = 2.0/(t-b);
  mP[1*4+2] = 0.0;
  mP[1*4+3] = 0.0;
  // column 2
  mP[2*4+0] = 0.0;
  mP[2*4+1] = 0.0;
  mP[2*4+2] = 2.0/(n-f);
  mP[2*4+3] = 0.0;
  // collumn 3
  mP[3*4+0] = -(l+r)/(r-l);
  mP[3*4+1] = -(t+b)/(t-b);
  mP[3*4+2] = -(n+f)/(n-f);
  mP[3*4+3] = 1.0;
}


DFM2_INLINE void solve_GlAffineMatrix(
    float vo[3],
    const float m[16],
    const float p[3])
{
  const float v[3] = {p[0]-m[3*4+0], p[1]-m[3*4+1], p[2]-m[3*4+2]};
  const float M[9] = {
    m[0*4+0],m[1*4+0],m[2*4+0],
    m[0*4+1],m[1*4+1],m[2*4+1],
    m[0*4+2],m[1*4+2],m[2*4+2] };
  float Minv[9];  InverseMat3(Minv, M);
  Mat3Vec(vo,
          Minv,v);
}

DFM2_INLINE void solve_GlAffineMatrixDirection(
    float vo[3],
    const float m[16],
    const float vi[3])
{
  const float M[9] = {
    m[0*4+0],m[1*4+0],m[2*4+0],
    m[0*4+1],m[1*4+1],m[2*4+1],
    m[0*4+2],m[1*4+2],m[2*4+2] };
  float Minv[9];  InverseMat3(Minv, M);
  Mat3Vec(vo,
          Minv,vi);
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


DFM2_INLINE void Mult_MatMat4(
    float *result,
    const float *m1,
    const float *m2)
{
  for(unsigned int i=0;i<4;++i){
    for(unsigned int j=0;j<4;++j){
      result[i * 4 + j] = 0.0;
      for(unsigned int k=0;k<4;++k) {
        result[i * 4 + j] += m1[i * 4 + k] * m2[k * 4 + j];
      }
    }
  }
}

DFM2_INLINE void CalcInvMat(
                            double *a,
                            const unsigned int n,
                            int &info)
{
  double tmp1;
  
  info = 0;
  unsigned int i, j, k;
  for (i = 0; i < n; i++) {
    if (fabs(a[i * n + i]) < 1.0e-30) {
      info = 1;
      return;
    }
    if (a[i * n + i] < 0.0) {
      info--;
    }
    tmp1 = 1.0 / a[i * n + i];
    a[i * n + i] = 1.0;
    for (k = 0; k < n; k++) {
      a[i * n + k] *= tmp1;
    }
    for (j = 0; j < n; j++) {
      if (j != i) {
        tmp1 = a[j * n + i];
        a[j * n + i] = 0.0;
        for (k = 0; k < n; k++) {
          a[j * n + k] -= tmp1 * a[i * n + k];
        }
      }
    }
  }
}


}
}

// static functions ends here
// =====================================================

template <typename REAL>
DFM2_INLINE void delfem2::MatMat4(
    REAL *m0,
    const REAL *m1,
    const REAL *m2)
{
  for (unsigned int i = 0; i < 4; ++i) {
    for (unsigned int j = 0; j < 4; ++j) {
      m0[i * 4 + j]
          = m1[i * 4 + 0] * m2[0 * 4 + j]
          + m1[i * 4 + 1] * m2[1 * 4 + j]
          + m1[i * 4 + 2] * m2[2 * 4 + j]
          + m1[i * 4 + 3] * m2[3 * 4 + j];
    }
  }
}
#ifndef DFM2_HEADER_ONLY
template void delfem2::MatMat4(double*, const double*, const double*);
template void delfem2::MatMat4(float*, const float*, const float*);
#endif

template <typename REAL>
DFM2_INLINE void delfem2::Mult_MatVec4(
    REAL *mv,
    const REAL *m,
    const REAL *v)
{
  for (unsigned int i = 0; i < 4; ++i) {
    mv[i] = 0;
    for (unsigned int j = 0; j < 4; ++j) {
      mv[i] += m[i * 4 + j] * v[j];
    }
  }
}

template <typename REAL>
void delfem2::Mult_VecMat4(
    REAL *mv,
    const REAL *v,
    const REAL *m)
{
  for (unsigned int i = 0; i < 4; ++i) {
    mv[i] = 0;
    for (unsigned int j = 0; j < 4; ++j) {
      mv[i] += v[j] * m[j * 4 + i];
    }
  }
}
#ifndef DFM2_HEADER_ONLY
template void delfem2::Mult_VecMat4(double*, const double*, const double*);
template void delfem2::Mult_VecMat4(float*, const float*, const float*);
#endif

DFM2_INLINE void delfem2::Inverse_Mat4(
    double minv[16],
    const double m[16])
{
  for(int i=0;i<16;++i){ minv[i] = m[i]; }
  int info; camera::CalcInvMat(minv,4,info);
}

DFM2_INLINE void delfem2::Mat4_AffineTransProjectionFrustum(
    float *matrix,
    float left,
    float right,
    float bottom,
    float top,
    float znear,
    float zfar)
{
  float temp, temp2, temp3, temp4;
  temp = 2.f * znear;
  temp2 = right - left;
  temp3 = top - bottom;
  temp4 = zfar - znear;
  // column 0
  matrix[0*4+0] = temp / temp2;
  matrix[0*4+1] = 0.0;
  matrix[0*4+2] = 0.0;
  matrix[0*4+3] = 0.0;
  // column 1
  matrix[1*4+0] = 0.0;
  matrix[1*4+1] = temp / temp3;
  matrix[1*4+2] = 0.0;
  matrix[1*4+3] = 0.0;
  // column 2
  matrix[2*4+0] = (right + left) / temp2;
  matrix[2*4+1] = (top + bottom) / temp3;
  matrix[2*4+2] = (-zfar - znear) / temp4;
  matrix[2*4+3] = -1.0;
  // column 3
  matrix[3*4+0] = 0.0;
  matrix[3*4+1] = 0.0;
  matrix[3*4+2] = (-temp * zfar) / temp4;
  matrix[3*4+3] = 0.0;
}

//matrix will receive the calculated perspective matrix.
//You would have to upload to your shader
// or use glLoadMatrixf if you aren't using shaders.
DFM2_INLINE void delfem2::Mat4_AffineTransProjectionPerspective(
    float *matrix,
    float fovyInDegrees,
    float aspectRatio,
    float znear,
    float zfar)
{
  float ymax, xmax;
  ymax = znear * tanf(fovyInDegrees * 3.14159 / 360.0);
  xmax = ymax * aspectRatio;
  Mat4_AffineTransProjectionFrustum(matrix,
      -xmax, xmax, -ymax, ymax, znear, zfar);
}


DFM2_INLINE void delfem2::MultMat4AffineTransTranslateFromRight(
    float *matrix,
    float x,
    float y,
    float z)
{
  matrix[12]=matrix[0]*x+matrix[4]*y+matrix[8]*z+matrix[12];
  matrix[13]=matrix[1]*x+matrix[5]*y+matrix[9]*z+matrix[13];
  matrix[14]=matrix[2]*x+matrix[6]*y+matrix[10]*z+matrix[14];
  matrix[15]=matrix[3]*x+matrix[7]*y+matrix[11]*z+matrix[15];
}


DFM2_INLINE void delfem2::Mat4_AffineTransLookAt(
    float* Mr,
    float eyex, float eyey, float eyez,
    float cntx, float cnty, float cntz,
    float upx, float upy, float upz )
{
  const float eyePosition3D[3] = {eyex, eyey, eyez};
  const float center3D[3] = {cntx, cnty, cntz};
  const float upVector3D[3] = {upx, upy, upz};
  // ------------------
  float forward[3] = {
    center3D[0] - eyePosition3D[0],
    center3D[1] - eyePosition3D[1],
    center3D[2] - eyePosition3D[2] };
  camera::Normalize3D(forward);
  // ------------------
  // Side = forward x up
  float side[3] = {1,0,0};
  camera::Cross3D(side, forward, upVector3D);
  camera::Normalize3D(side);
  // ------------------
  //Recompute up as: up = side x forward
  float up[3] = {0,1,0};
  camera::Cross3D(up, side, forward);
  // ------------------
  Mr[ 0] = side[0];
  Mr[ 4] = side[1];
  Mr[ 8] = side[2];
  Mr[12] = 0.0;
  // ------------------
  Mr[ 1] = up[0];
  Mr[ 5] = up[1];
  Mr[ 9] = up[2];
  Mr[13] = 0.0;
  // ------------------
  Mr[ 2] = -forward[0];
  Mr[ 6] = -forward[1];
  Mr[10] = -forward[2];
  Mr[14] = 0.0;
  // ------------------
  Mr[ 3] = 0.0;
  Mr[ 7] = 0.0;
  Mr[11] = 0.0;
  Mr[15] = 1.0;
  // ------------------
  delfem2::MultMat4AffineTransTranslateFromRight(Mr,
      -eyePosition3D[0], -eyePosition3D[1], -eyePosition3D[2]);
}

DFM2_INLINE void delfem2::screenUnProjection(
    float vout[3],
    const float v[3],
    const float mMV[16],
    const float mPj[16])
{
  const float D = mPj[11] + mPj[15]; // z is 1 after model view
  const float v0[3] = { D*v[0], D*v[1], 0.0 };
  float v1[3];
  camera::solve_GlAffineMatrix(v1,
      mPj, v0);
  v1[2] = 1;
  camera::solve_GlAffineMatrix(vout,
      mMV, v1);
}


void delfem2::screenUnProjectionDirection(
    float vo[3],
    const float vi[3],
    const float mMV[16],
    const float mPj[16])
{
  float v1[3];
  camera::solve_GlAffineMatrixDirection(v1,
      mPj, vi);
  camera::solve_GlAffineMatrixDirection(vo,
      mMV, v1);
  camera::Normalize3D(vo);
}


// =================================================
// implementation of CCamera class starts here

template <typename REAL>
void delfem2::CCamera<REAL>::Mat4_AffineTransProjection(
    float mP[16],
    double asp,
    double depth) const
{
  if( is_pars ){
    Mat4_AffineTransProjectionPerspective(mP,
                                          fovy, asp, depth*0.01, depth*10);
  }
  else{
    camera::Mat4_AffineTransProjectionOrtho(mP,
               -view_height*asp,
               +view_height*asp,
               -view_height,
               +view_height,
               -depth*10,
               +depth*10);
  }
}
template void delfem2::CCamera<double>::Mat4_AffineTransProjection(
    float mP[16],
    double asp,
    double depth) const;

// ----------------------------

/**
 *
 * @tparam REAL
 * @param mMV model view matrix (column major order)
 * @detail column major
 */
template <typename REAL>
void delfem2::CCamera<REAL>::Mat4_AffineTransModelView(float mMV[16]) const
{
  float Mt[16];
  {
    const float transd[3] = { (float)trans[0], (float)trans[1], (float)trans[2] };
    camera::Mat4_AffineTransTranslate(Mt, transd);
  }
  const float Ms[16] = {
      (float)scale, 0, 0, 0,
      0, (float)scale, 0, 0,
      0, 0, (float)scale, 0,
      0, 0, 0, 1 };
  float Mr[16];
  {
    if (camera_rot_mode == CAMERA_ROT_MODE::YTOP) {
      double x = sin(theta);
      double z = cos(theta);
      double y = sin(psi);
      x *= cos(psi);
      z *= cos(psi);
      Mat4_AffineTransLookAt(Mr, x, y, z, 0, 0, 0, 0, 1, 0);
    } else if (camera_rot_mode == CAMERA_ROT_MODE::ZTOP) {
      double x = sin(theta);
      double y = cos(theta);
      double z = sin(psi);
      x *= cos(psi);
      y *= cos(psi);
      Mat4_AffineTransLookAt(Mr, x, y, z, 0, 0, 0, 0, 0, 1);
    } else if (camera_rot_mode == CAMERA_ROT_MODE::TBALL) {
      const float q[4] = {
          static_cast<float>(Quat_tball[0]),
          static_cast<float>(Quat_tball[1]),
          static_cast<float>(Quat_tball[2]),
          static_cast<float>(Quat_tball[3])};
      camera::Mat4_AffineTransQuat(Mr, q);
    }
  }
  float Mrt[16]; camera::Mult_MatMat4(Mrt, Mr,Mt);
  camera::Mult_MatMat4(mMV, Mrt,Ms);
}
template void delfem2::CCamera<double>::Mat4_AffineTransModelView(float mMV[16]) const;

// -------------------------

template <typename REAL>
void delfem2::CCamera<REAL>::Scale(double s){
  scale *= s;
}
#ifndef DFM2_HEADER_ONLY
template void delfem2::CCamera<double>::Scale(double s);
#endif
  
// --------------------------

template <typename REAL>
void delfem2::CCamera<REAL>::Rot_Camera(double dx, double dy){
  if(      camera_rot_mode == CAMERA_ROT_MODE::YTOP ){
    theta -= dx;
    psi   -= dy;
  }
  else if( camera_rot_mode == CAMERA_ROT_MODE::ZTOP ){
    theta += dx;
    psi   -= dy;
  }
  else if( camera_rot_mode == CAMERA_ROT_MODE::TBALL ){
    double a = sqrt(dx * dx + dy * dy);
    double ar = a*0.5; // angle
    double dq[4] = { cos(ar), -dy*sin(ar)/a, dx*sin(ar)/a, 0.0 };
    if (a != 0.0) {
      double qtmp[4]; camera::QuatQuat(qtmp, dq, Quat_tball);
      camera::CopyQuat(Quat_tball,qtmp);
    }
  }
}
#ifndef DFM2_HEADER_ONLY
template void delfem2::CCamera<double>::Rot_Camera(double dx, double dy);
#endif

// ------------------------

template <typename REAL>
void delfem2::CCamera<REAL>::Pan_Camera(double dx, double dy){
    double s = view_height/scale;
    trans[0] += s*dx;
    trans[1] += s*dy;
    trans[2] += s*0.0;
  
}
#ifndef DFM2_HEADER_ONLY
template void delfem2::CCamera<double>::Pan_Camera(double dx, double dy);
#endif

// ------------------------

namespace delfme2 {

template <typename REAL>
std::ostream &operator<<(std::ostream &output, delfem2::CCamera<REAL>& c)
{
  output.setf(std::ios::scientific);
  output << c.is_pars << std::endl;
  output << c.fovy << std::endl;
  output << c.scale << std::endl;
  output << c.view_height << std::endl;
  output << c.trans[0] << " " << c.trans[1] << " " << c.trans[2] << std::endl;
  output << (int)c.camera_rot_mode << std::endl;
  output << c.theta << " " << c.psi << std::endl;
  output << c.Quat_tball[0] << " " << c.Quat_tball[1] << " " << c.Quat_tball[2] << " " << c.Quat_tball[3] << std::endl;
  return output;
}

template <typename REAL>
std::istream &operator>>(std::istream &input, delfem2::CCamera<REAL>& c)
{
  {
    int is0; input >> is0; c.is_pars = (bool)is0;
  }
  input >> c.fovy;
  input >> c.scale;
  input >> c.view_height;
  input >> c.trans[0] >> c.trans[1] >> c.trans[2];
  {
    int imode;
    input >> imode;
    c.camera_rot_mode = imode;//(delfem2::CCamera<REAL>::CAMERA_ROT_MODE)imode;
  }
  input >> c.theta >> c.psi;
  input >> c.Quat_tball[0] >> c.Quat_tball[1] >> c.Quat_tball[2] >> c.Quat_tball[3];
  return input;
}

}
