/*
 * Copyright (c) 2019 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */


#include <cmath>
#include <cstdio> // memcpy
#include <cstring>
#include <iostream>

#include "delfem2/camera.h"

namespace dfm2 = delfem2;

// -----------------------------------------------------------------------


//! @brief multiply two quaternion
static void QuatQuat(double r[], const double p[], const double q[])
{
  r[0] = p[0] * q[0] - p[1] * q[1] - p[2] * q[2] - p[3] * q[3];
  r[1] = p[0] * q[1] + p[1] * q[0] + p[2] * q[3] - p[3] * q[2];
  r[2] = p[0] * q[2] - p[1] * q[3] + p[2] * q[0] + p[3] * q[1];
  r[3] = p[0] * q[3] + p[1] * q[2] - p[2] * q[1] + p[3] * q[0];
}

//! @brief transform vector with quaternion
inline void QuatVec(double vo[], const double q[], const double vi[])
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
inline void QuatConjVec(double vo[], const double q[], const double vi[])
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
static void CopyQuat
 (double r[],
  const double p[])
{
  r[0] = p[0];
  r[1] = p[1];
  r[2] = p[2];
  r[3] = p[3];
}

static void Mat4f_Quat(float r[], const double q[])
{
  auto x2 = (float)(q[1] * q[1] * 2.0);
  auto y2 = (float)(q[2] * q[2] * 2.0);
  auto z2 = (float)(q[3] * q[3] * 2.0);
  auto xy = (float)(q[1] * q[2] * 2.0);
  auto yz = (float)(q[2] * q[3] * 2.0);
  auto zx = (float)(q[3] * q[1] * 2.0);
  auto xw = (float)(q[1] * q[0] * 2.0);
  auto yw = (float)(q[2] * q[0] * 2.0);
  auto zw = (float)(q[3] * q[0] * 2.0);
  
  r[ 0] = 1.f - y2 - z2;
  r[ 1] = xy + zw;
  r[ 2] = zx - yw;
  r[ 4] = xy - zw;
  r[ 5] = 1.f - z2 - x2;
  r[ 6] = yz + xw;
  r[ 8] = zx + yw;
  r[ 9] = yz - xw;
  r[10] = 1.f - x2 - y2;
  r[ 3] = r[ 7] = r[11] = r[12] = r[13] = r[14] = 0.0;
  r[15] = 1.0;
}


void Mat4f_Identity(float r[])
{
  r[ 0] = 1.0;  r[ 1] = 0.0;  r[ 2] = 0.0;  r[ 3] = 0.0;
  r[ 4] = 0.0;  r[ 5] = 1.0;  r[ 6] = 0.0;  r[ 7] = 0.0;
  r[ 8] = 0.0;  r[ 9] = 0.0;  r[10] = 1.0;  r[11] = 0.0;
  r[12] = 0.0;  r[13] = 0.0;  r[14] = 0.0;  r[15] = 1.0;
}


static void Normalize3D(float vec[3])
{
  float len = sqrt(vec[0]*vec[0]+vec[1]*vec[1]+vec[2]*vec[2]);
  float leninv = 1.0/len;
  vec[0] *= leninv;
  vec[1] *= leninv;
  vec[2] *= leninv;
}

static void Cross3D(float r[3], const float v1[3], const float v2[3]){
  r[0] = v1[1]*v2[2] - v2[1]*v1[2];
  r[1] = v1[2]*v2[0] - v2[2]*v1[0];
  r[2] = v1[0]*v2[1] - v2[0]*v1[1];
}


// static functions ends here
//----------------------------------------------------------------------------------------

void dfm2::glhFrustumf2
(float *matrix,
 float left, float right,
 float bottom, float top,
 float znear, float zfar)
{
  float temp, temp2, temp3, temp4;
  temp = 2.f * znear;
  temp2 = right - left;
  temp3 = top - bottom;
  temp4 = zfar - znear;
  matrix[0] = temp / temp2;
  matrix[1] = 0.0;
  matrix[2] = 0.0;
  matrix[3] = 0.0;
  matrix[4] = 0.0;
  matrix[5] = temp / temp3;
  matrix[6] = 0.0;
  matrix[7] = 0.0;
  matrix[8] = (right + left) / temp2;
  matrix[9] = (top + bottom) / temp3;
  matrix[10] = (-zfar - znear) / temp4;
  matrix[11] = -1.0;
  matrix[12] = 0.0;
  matrix[13] = 0.0;
  matrix[14] = (-temp * zfar) / temp4;
  matrix[15] = 0.0;
}

//matrix will receive the calculated perspective matrix.
//You would have to upload to your shader
// or use glLoadMatrixf if you aren't using shaders.
void dfm2::glhPerspectivef2
(float *matrix,
 float fovyInDegrees, float aspectRatio,
 float znear, float zfar)
{
  float ymax, xmax;
//  float temp, temp2, temp3, temp4;
  ymax = znear * tanf(fovyInDegrees * 3.14159 / 360.0);
  //ymin = -ymax;
  //xmin = -ymax * aspectRatio;
  xmax = ymax * aspectRatio;
  glhFrustumf2(matrix, -xmax, xmax, -ymax, ymax, znear, zfar);
}


void dfm2::glhTranslatef2
(float *matrix,
 float x, float y, float z)
{
  matrix[12]=matrix[0]*x+matrix[4]*y+matrix[8]*z+matrix[12];
  matrix[13]=matrix[1]*x+matrix[5]*y+matrix[9]*z+matrix[13];
  matrix[14]=matrix[2]*x+matrix[6]*y+matrix[10]*z+matrix[14];
  matrix[15]=matrix[3]*x+matrix[7]*y+matrix[11]*z+matrix[15];
}

//PURPOSE:      For square matrices. This is column major for OpenGL
inline void MultiplyMatrices4by4OpenGL_FLOAT
 (float *result,
  const float *m1,
  const float *m2)
{
  result[ 0]=m1[0]*m2[0]+  m1[4]*m2[1]+  m1[8]*m2[2]+  m1[12]*m2[3];
  result[ 4]=m1[0]*m2[4]+  m1[4]*m2[5]+  m1[8]*m2[6]+  m1[12]*m2[7];
  result[ 8]=m1[0]*m2[8]+  m1[4]*m2[9]+  m1[8]*m2[10]+  m1[12]*m2[11];
  result[12]=m1[0]*m2[12]+   m1[4]*m2[13]+  m1[8]*m2[14]+  m1[12]*m2[15];
  
  result[ 1]=m1[1]*m2[0]+  m1[5]*m2[1]+  m1[9]*m2[2]+  m1[13]*m2[3];
  result[ 5]=m1[1]*m2[4]+  m1[5]*m2[5]+  m1[9]*m2[6]+  m1[13]*m2[7];
  result[ 9]=m1[1]*m2[8]+  m1[5]*m2[9]+  m1[9]*m2[10]+  m1[13]*m2[11];
  result[13]=m1[1]*m2[12]+  m1[5]*m2[13]+  m1[9]*m2[14]+  m1[13]*m2[15];
  
  result[ 2]=m1[2]*m2[0]+  m1[6]*m2[1]+  m1[10]*m2[2]+  m1[14]*m2[3];
  result[ 6]=m1[2]*m2[4]+  m1[6]*m2[5]+  m1[10]*m2[6]+  m1[14]*m2[7];
  result[10]=m1[2]*m2[8]+  m1[6]*m2[9]+  m1[10]*m2[10]+  m1[14]*m2[11];
  result[14]=m1[2]*m2[12]+   m1[6]*m2[13]+  m1[10]*m2[14]+  m1[14]*m2[15];
  
  result[ 3]=m1[3]*m2[0]+  m1[7]*m2[1]+  m1[11]*m2[2]+  m1[15]*m2[3];
  result[ 7]=m1[3]*m2[4]+  m1[7]*m2[5]+  m1[11]*m2[6]+  m1[15]*m2[7];
  result[11]=m1[3]*m2[8]+   m1[7]*m2[9]+  m1[11]*m2[10]+  m1[15]*m2[11];
  result[15]=m1[3]*m2[12]+   m1[7]*m2[13]+  m1[11]*m2[14]+  m1[15]*m2[15];
}

void dfm2::glhLookAtf2
( float *matrix,
 float eyex, float eyey, float eyez,
 float cntx, float cnty, float cntz,
 float upx, float upy, float upz )
{
  float eyePosition3D[3] = {eyex, eyey, eyez};
  float center3D[3] = {cntx, cnty, cntz};
  float upVector3D[3] = {upx, upy, upz};
  
  
  //------------------
  float forward[3] = {0,0,1};
  forward[0] = center3D[0] - eyePosition3D[0];
  forward[1] = center3D[1] - eyePosition3D[1];
  forward[2] = center3D[2] - eyePosition3D[2];
  Normalize3D(forward);
  //------------------
  //Side = forward x up
  float side[3] = {1,0,0};
  Cross3D(side, forward, upVector3D);
  Normalize3D(side);
  //------------------
  //Recompute up as: up = side x forward
  float up[3] = {0,1,0};
  Cross3D(up, side, forward);
  //------------------
  float matrix2[16];
  matrix2[0] = side[0];
  matrix2[4] = side[1];
  matrix2[8] = side[2];
  matrix2[12] = 0.0;
  //------------------
  matrix2[1] = up[0];
  matrix2[5] = up[1];
  matrix2[9] = up[2];
  matrix2[13] = 0.0;
  //------------------
  matrix2[2] = -forward[0];
  matrix2[6] = -forward[1];
  matrix2[10] = -forward[2];
  matrix2[14] = 0.0;
  //------------------
  matrix2[3] = matrix2[7] = matrix2[11] = 0.0;
  matrix2[15] = 1.0;
  //------------------
  float resultMatrix[16];
  MultiplyMatrices4by4OpenGL_FLOAT(resultMatrix, matrix, matrix2);
  dfm2::glhTranslatef2(resultMatrix,
                 -eyePosition3D[0], -eyePosition3D[1], -eyePosition3D[2]);
  //------------------
  memcpy(matrix, resultMatrix, 16*sizeof(float));
}


void glhOrthof2
(float mP[16],
 double l, double r,
 double b, double t,
 double n, double f)
{
  mP[0*4+0] = 2.0/(r-l);
  mP[0*4+1] = 0.0;
  mP[0*4+2] = 0.0;
  mP[0*4+3] = -(l+r)/(r-l);
  
  mP[1*4+0] = 0.0;
  mP[1*4+1] = 2.0/(t-b);
  mP[1*4+2] = 0.0;
  mP[1*4+3] = -(t+b)/(t-b);
  
  mP[2*4+0] = 0.0;
  mP[2*4+1] = 0.0;
  mP[2*4+2] = 2.0/(n-f);
  mP[2*4+3] = -(n+f)/(n-f);
  
  mP[3*4+0] = 0.0;
  mP[3*4+1] = 0.0;
  mP[3*4+2] = 0.0;
  mP[3*4+3] = 1.0;
}


void InverseMat3
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

void Mat3Vec
 (float vo[3],
  const float mat[9],
  const float vi[3])
{
  vo[0] = mat[0]*vi[0] + mat[1]*vi[1] + mat[2]*vi[2];
  vo[1] = mat[3]*vi[0] + mat[4]*vi[1] + mat[5]*vi[2];
  vo[2] = mat[6]*vi[0] + mat[7]*vi[1] + mat[8]*vi[2];
}

void solve_GlAffineMatrix
 (float vo[3],
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

void solve_GlAffineMatrixDirection
(float vo[3],
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

void dfm2::screenUnProjection
 (float vout[3],
  const float v[3],
  const float mMV[16],
  const float mPj[16])
{
  const float D = mPj[11] + mPj[15]; // z is 1 after model view
  const float v0[3] = { D*v[0], D*v[1], 0.0 };
  float v1[3];
  solve_GlAffineMatrix(v1,
                       mPj, v0);
  v1[2] = 1;
  solve_GlAffineMatrix(vout,
                       mMV, v1);
}


void dfm2::screenUnProjectionDirection
(float vo[3],
 const float vi[3],
 const float mMV[16],
 const float mPj[16])
{
//  CVector3 v0 = solve_GlAffineMatrixDirection(mPj, v);
//  CVector3 v1 = solve_GlAffineMatrixDirection(mMV, v0);
//  v1.SetNormalizedVector();
//  return v1;
  float v1[3];
  solve_GlAffineMatrixDirection(v1,
                                mPj, vi);
  v1[2] = 1;
  solve_GlAffineMatrixDirection(vo,
                                mMV, v1);
  Normalize3D(vo);
}


// ----------------------------------------------------
// implementation of CCamera class starts here


void dfm2::CCamera::Affine4f_Projection
(float mP[16], double asp, double depth) const 
{
  if( is_pars ){
    glhPerspectivef2(mP, fovy, asp, depth*0.01, depth*10);
  }
  else{
    ////
    glhOrthof2(mP,
               -view_height/scale*asp,
               +view_height/scale*asp,
               -view_height/scale,
               +view_height/scale,
               -depth*10,
               +depth*10);
  }
}

void dfm2::CCamera::Affine4f_ModelView
(float mMV[16]) const
{
  /*
  {
    ::glMatrixMode(GL_MODELVIEW);
    ::glLoadIdentity();
    ::glTranslated(trans[0],trans[1],trans[2]);
  }
   */
  const float B[16] =
  { 1.f, 0.f, 0.f, 0.f,
    0.f, 1.f, 0.f, 0.f,
    0.f, 0.f, 1.f, 0.f,
    (float)trans[0], (float)trans[1], (float)trans[2], 1.f};
  float A[16];
  Mat4f_Identity(A);
  if(      camera_rot_mode == CAMERA_ROT_YTOP  ){
    double x = sin(theta);
    double z = cos(theta);
    double y = sin(psi);
    x *= cos(psi);
    z *= cos(psi);
    glhLookAtf2(A, x,y,z, 0,0,0, 0,1,0);
  }
  else if( camera_rot_mode == CAMERA_ROT_ZTOP  ){
    double x = sin(theta);
    double y = cos(theta);
    double z = sin(psi);
    x *= cos(psi);
    y *= cos(psi);
    glhLookAtf2(A, x,y,z, 0,0,0, 0,0,1);
  }
  else if( camera_rot_mode == CAMERA_ROT_TBALL ){
    Mat4f_Quat(A,Quat_tball);
  }
  MultiplyMatrices4by4OpenGL_FLOAT(mMV, B,A);
}

void dfm2::CCamera::Scale(double s){
  scale *= s;
}
void dfm2::CCamera::Rot_Camera(double dx, double dy){
    if(      camera_rot_mode == CAMERA_ROT_YTOP ){
      theta -= dx;
      psi   -= dy;
    }
    else if( camera_rot_mode == CAMERA_ROT_ZTOP ){
      theta += dx;
      psi   -= dy;
    }
    else if( camera_rot_mode == CAMERA_ROT_TBALL ){
      double a = sqrt(dx * dx + dy * dy);
      double ar = a*0.5; // angle
      double dq[4] = { cos(ar), -dy*sin(ar)/a, dx*sin(ar)/a, 0.0 };
      if (a != 0.0) {
        double qtmp[4]; QuatQuat(qtmp, dq, Quat_tball);
        CopyQuat(Quat_tball,qtmp);
      }
    }
}

void dfm2::CCamera::Pan_Camera(double dx, double dy){
    double s = view_height/scale;
    trans[0] += s*dx;
    trans[1] += s*dy;
    trans[2] += s*0.0;
  
}

std::ostream &operator<<(std::ostream &output, dfm2::CCamera& c)
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

std::istream &operator>>(std::istream &input, dfm2::CCamera& c)
{
  {
    int is0; input >> is0; c.is_pars = (bool)is0;
  }
  input >> c.fovy;
  input >> c.scale;
  input >> c.view_height;
  input >> c.trans[0] >> c.trans[1] >> c.trans[2];
  {
    int imode; input >> imode; c.camera_rot_mode = (dfm2::CAMERA_ROT_MODE)imode;
  }
  input >> c.theta >> c.psi;
  input >> c.Quat_tball[0] >> c.Quat_tball[1] >> c.Quat_tball[2] >> c.Quat_tball[3];
  return input;
}

