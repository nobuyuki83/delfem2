/*
 * Copyright (c) 2019 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */

// DONE(2020/12/11): rename the class "CCam3_OnAxisZplus" from CCvCamera
// DONE(2020/12/--): rename this file "cam3_m4q.h"
// DONE(2020/11/29): make this file dependent on mat4.h and quat.h

/**
 * @file camera class that define projection and model-view transformation.
 * @details this file should be stand alone
 */

#ifndef DFM2_CAMERA_H
#define DFM2_CAMERA_H

#include <iostream>
#include <cmath>
#include <cstdio> // memcpy
#include <cstring>

#include "delfem2/quat.h"
#include "delfem2/mat4.h"
#include "delfem2/dfm2_inline.h"

// ---------------------------------------------------

namespace delfem2{

// ----------------------------------------------------

/**
 * the camera is placed on the positive Z axis and looking at the origin.
 * The object will undergo movdel view transformation (rotation->translation)
 * and then perspective transformation (scale->orthognoal/perspective).
 *
 * @tparam REAL either float or double
 */
template <typename REAL>
class CCam3_OnAxisZplusLookOrigin
{

public:
  enum class CAMERA_ROT_MODE { YTOP, ZTOP, TBALL };

public:
  CCam3_OnAxisZplusLookOrigin() :
  trans{0,0,0}, Quat_tball{0,0,0,1}
  {
    is_pars = false;
    fovy = 10;
    view_height = 1.0;
    scale = 1.0;
    camera_rot_mode =  CAMERA_ROT_MODE::TBALL;
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
  void Mat4_AffineTransProjection(float mP[16], float asp) const;
  
  /**
   *
   * @tparam REAL
   * @param mMV model view matrix (column major order)
   * @detail column major
   */
  void Mat4_AffineTransModelView(float mMV[16]) const;

  /**
   * @brief make 4x4 affine matrix for view transformation.
   * @details the matrix strage format is *column-major*. Applying trasformation need to multply vector from *left hand side*.
   * A 3D point is transfromed with this affine matrix and then a cube [-1,+1, -1,+1, -1,+1] is looked from -Z directoin.
   * To look from +Z direction, The transformation needs a mirror transformation in XY plane.
   * We separate mMV and  mP because of the light (light position should not be affected by the modelview transform).
   * @param[out] mMV modelview matrix (column major order)
   * @param[out] mP  projection matrix (column major order)
   * @param asp
   */
  void Mat4_MVP_OpenGL(
      float mMV[16],
      float mP[16],
      float asp) const
  {
    Mat4_AffineTransProjection(mP, asp); // project space into cube [-1,+1,-1,+1,-1,+1] and view from -Z
    Mat4_AffineTransModelView(mMV);
  }
  // ------------------------
  
  void Scale(double s);
  void Rot_Camera(REAL dx, REAL dy);
  void Pan_Camera(REAL dx, REAL dy);
private:
public:
  bool is_pars;
  double fovy;
  double scale;
  double view_height;

  double trans[3];
  
  CAMERA_ROT_MODE camera_rot_mode;
  
  // ytop
  REAL theta;
  REAL psi;
  
  // tball
  double Quat_tball[4];
};

template <typename REAL>
std::ostream &operator<<(std::ostream &output, CCam3_OnAxisZplusLookOrigin<REAL>& c);

template <typename REAL>
std::istream &operator>>(std::istream &input, CCam3_OnAxisZplusLookOrigin<REAL>& c);

} // namespace delfem2

#ifndef DFM2_STATIC_LIBRARY
#  include "delfem2/cam3_m4q.cpp"
#endif

#endif
