/*
 * Copyright (c) 2019 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */


#include <cmath>
#include <iostream>

#include "delfem2/cam3_m4q.h"

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
    Mat4_AffineTransProjectionOrtho(mP,
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

template <typename REAL>
void delfem2::CCamera<REAL>::Mat4_AffineTransModelView(float mMV[16]) const
{
  float Mt[16];
  {
    const float transd[3] = { (float)trans[0], (float)trans[1], (float)trans[2] };
    Mat4_AffineTransTranslate(Mt, transd);
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
      Mat4_AffineTransQuat(Mr, q);
    }
  }
  float Mrt[16]; MatMat4(Mrt, Mr,Mt);
  MatMat4(mMV, Mrt,Ms);
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
      double qtmp[4]; QuatQuat(qtmp, dq, Quat_tball);
      Copy_Quat(Quat_tball,qtmp);
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
