/*
 * Copyright (c) 2019 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */

// TODO: rename this file "cam3_m4q.h"
// DONE(2020_11_29): make this file dependent on mat4.h and quat.h

/**
 * @file camera class that define projection and model-view transformation.
 * @details this file should be stand alone
 */

#ifndef DFM2_CAMERA_H
#define DFM2_CAMERA_H

#include "delfem2/quat.h"
#include "delfem2/mat4.h"
#include "delfem2/dfm2_inline.h"
#include <ostream>
#include <cmath>
#include <cstdio> // memcpy
#include <cstring>

// ---------------------------------------------------

namespace delfem2{

// ----------------------------------------------------


template <typename REAL>
class CCamera
{
public:
  enum class CAMERA_ROT_MODE { YTOP, ZTOP, TBALL };

public:
  CCamera() :
  trans{0,0,0}, Quat_tball{1,0,0,0}
  {
    is_pars = false;
    fovy = 60;
    view_height = 1.0;
    scale = 1.0;
    camera_rot_mode =  CAMERA_ROT_MODE::YTOP;
    psi = 0;
    theta = 0;
  }
  
  // -----------------------
  // const methods from here
  
  /**
   * @brief Compute "projection matrix" for OpenGL.
   * @detial OpenGL will draw the object in the cube [-1,+1, -1,+1, -1,+1] looking from -Z
   * The projection matrix will mirror the object Z.
   * @tparam REAL
   * @param mP
   * @param asp
   * @param depth
   */
  void Mat4_AffineTransProjection(float mP[16], double asp, double depth) const;
  
  /**
   *
   * @tparam REAL
   * @param mMV model view matrix (column major order)
   * @detail column major
   */
  void Mat4_AffineTransModelView(float mMV[16]) const;
  
  // ------------------------
  
  void Scale(double s);
  void Rot_Camera(double dx, double dy);
  void Pan_Camera(double dx, double dy);
private:
public:
  bool is_pars;
  double fovy;
  double scale;
  double view_height;
  double trans[3];
  
  CAMERA_ROT_MODE camera_rot_mode;
  
  // ytop
  double theta;
  double psi;
  
  // tball
  double Quat_tball[4];
};

template <typename REAL>
std::ostream &operator<<(std::ostream &output, CCamera<REAL>& c);

template <typename REAL>
std::istream &operator>>(std::istream &input, CCamera<REAL>& c);

} // namespace delfem2

#ifdef DFM2_HEADER_ONLY
#  include "delfem2/camera.cpp"
#endif

#endif
