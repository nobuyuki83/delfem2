/*
 * Copyright (c) 2019 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */


#include <cmath>
#include <iostream>

#include "delfem2/cam3_m4q.h"

#ifndef M_PI
# define M_PI 3.14159265359
#endif

// =================================================
// implementation of CCvCamera class starts here

template <typename REAL>
void delfem2::CCam3_OnAxisZplusLookOrigin<REAL>::Mat4_AffineTransProjection(
    float mP[16],
    float asp) const
{
  const float mS[16] = {
      (float)scale, 0, 0, 0,
      0, (float)scale, 0, 0,
      0, 0, (float)scale, 0,
      0, 0, 0, 1 };
  float depth = (0.5*view_height)/tan(fovy*0.5*(2*M_PI)/360.0);
  float mP0[16];
  if( is_pars ){
    Mat4_AffineTransProjectionFrustum(mP0,
        fovy*(2*M_PI)/360.0, asp,
        -depth*2, -depth*0.01);
  }
  else{
    Mat4_AffineTransProjectionOrtho(mP0,
        -view_height*asp,
        +view_height*asp,
        -view_height,
        +view_height,
        -2*depth,
        0);
  }
  float mT0[16];
  {
    // the camera is placed at the origin and lookin into the -Z direction in the range [-2*depth,0]
    // to view the object we translate the object at the origin (0,0,-depth)
    const float t0[3] = {0.f,0.f,-(float)depth };
    ::delfem2::Mat4_AffineTransTranslate(mT0, t0 );
  }
  const float mRefZ[16] = { // reflection with the XY plane
      +1.f,  0.f,  0.f,  0.f,
       0.f, +1.f,  0.f,  0.f,
       0.f,  0.f, -1.f,  0.f,
       0.f,  0.f,  0.f, +1.f };
  float mTmp1[16]; ::delfem2::MatMat4(mTmp1,mS,mT0);
  float mTmp0[16]; ::delfem2::MatMat4(mTmp0,mTmp1,mP0);
  ::delfem2::MatMat4(mP,mTmp0,mRefZ);
}
template void delfem2::CCam3_OnAxisZplusLookOrigin<double>::Mat4_AffineTransProjection(
    float mP[16],
    float asp) const;

// ----------------------------

template <typename REAL>
void delfem2::CCam3_OnAxisZplusLookOrigin<REAL>::Mat4_AffineTransModelView(float mMV[16]) const
{
  float Mt[16];
  {
    const float transd[3] = { (float)trans[0], (float)trans[1], (float)trans[2] };
    Mat4_AffineTransTranslate(Mt, transd);
  }
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
      Mat4_AffineTransQuat(Mr, q);
    }
  }
  MatMat4(mMV, Mr,Mt);
}
template void delfem2::CCam3_OnAxisZplusLookOrigin<double>::Mat4_AffineTransModelView(float mMV[16]) const;

// -------------------------

template <typename REAL>
void delfem2::CCam3_OnAxisZplusLookOrigin<REAL>::Scale(double s){
  scale *= s;
}
#ifndef DFM2_HEADER_ONLY
template void delfem2::CCam3_OnAxisZplusLookOrigin<double>::Scale(double s);
#endif
  
// --------------------------

template <typename REAL>
void delfem2::CCam3_OnAxisZplusLookOrigin<REAL>::Rot_Camera(double dx, double dy){
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
      double qtmp[4]; QuatQuat(qtmp, dq, Quat_tball);
      Copy_Quat(Quat_tball,qtmp);
    }
  }
}
#ifndef DFM2_HEADER_ONLY
template void delfem2::CCam3_OnAxisZplusLookOrigin<double>::Rot_Camera(double dx, double dy);
#endif

// ------------------------

template <typename REAL>
void delfem2::CCam3_OnAxisZplusLookOrigin<REAL>::Pan_Camera(double dx, double dy){
    double s = view_height/scale;
    trans[0] += s*dx;
    trans[1] += s*dy;
    trans[2] += s*0.0;
  
}
#ifndef DFM2_HEADER_ONLY
template void delfem2::CCam3_OnAxisZplusLookOrigin<double>::Pan_Camera(double dx, double dy);
#endif

// ------------------------

namespace delfme2 {

template <typename REAL>
std::ostream &operator<<(std::ostream &output, delfem2::CCam3_OnAxisZplusLookOrigin<REAL>& c)
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
std::istream &operator>>(std::istream &input, delfem2::CCam3_OnAxisZplusLookOrigin<REAL>& c)
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
    c.camera_rot_mode = imode;//(delfem2::CCvCamera<REAL>::CAMERA_ROT_MODE)imode;
  }
  input >> c.theta >> c.psi;
  input >> c.Quat_tball[0] >> c.Quat_tball[1] >> c.Quat_tball[2] >> c.Quat_tball[3];
  return input;
}

}
