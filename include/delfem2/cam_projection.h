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

#ifndef DFM2_CAM_PROJECTION_H
#define DFM2_CAM_PROJECTION_H

#include <iostream>
#include <cmath>
#include <cstdio> // memcpy
#include <cstring>
#include <array>

#include "delfem2/quat.h"
#include "delfem2/mat4.h"
#include "delfem2/cvcamera.h"
#include "delfem2/dfm2_inline.h"

// ---------------------------------------------------

namespace delfem2 {

/**
 * Pure virtual class
 * @tparam REAL
 */
template <typename REAL>
class Projection {
public:
  virtual ~Projection()= default;
  /**
 * @brief Compute "projection matrix" for OpenGL.
 * @detial OpenGL will draw the object in the cube [-1,+1, -1,+1, -1,+1] looking from -Z
 * A 3D point is transfromed with this affine matrix and then a cube [-1,+1, -1,+1, -1,+1] is looked from -Z directoin.
 * The projection matrix will mirror the object Z.
 * To look from +Z direction, The transformation needs a mirror transformation in XY plane.
 * @param mP
 * @param asp
 */
  [[nodiscard]] virtual std::array<float,16> GetMatrix(float asp) const = 0;
};


/**
 * the camera is placed on the positive Z axis and looking at the origin.
 * The object will undergo movdel view transformation (rotation->translation)
 * and then perspective transformation (scale->orthognoal/perspective).
 *
 * @tparam REAL either float or double
 */
template<typename REAL>
class Projection_LookOriginFromZplus : public Projection<REAL>{
 public:
  explicit Projection_LookOriginFromZplus(double view_height = 1,
                                 bool is_pars = false,
                                 double fovy = 10)
  : view_height(view_height), is_pars(is_pars), fovy(fovy){}
  
  ~Projection_LookOriginFromZplus()= default;
  /**
   * @brief Compute "projection matrix" for OpenGL.
   * @detial OpenGL will draw the object in the cube [-1,+1, -1,+1, -1,+1] looking from -Z
   * A 3D point is transfromed with this affine matrix and then a cube [-1,+1, -1,+1, -1,+1] is looked from -Z directoin.
   * The projection matrix will mirror the object Z.
   * To look from +Z direction, The transformation needs a mirror transformation in XY plane.
   * @param mP
   * @param asp
   */
   [[nodiscard]] std::array<float,16> GetMatrix(float asp) const override;
public:
  double view_height = 1;
  bool is_pars = false;
  double fovy = 10;
};


template<typename REAL>
std::array<float,16> delfem2::Projection_LookOriginFromZplus<REAL>::GetMatrix(
    float asp) const {
  REAL fovyInRad = fovy * (2. * M_PI) / 360.f;
  REAL depth = view_height / tan(fovyInRad * 0.5f);
  REAL mP0[16];
  if (is_pars) {
    Mat4_AffineTransProjectionFrustum(
        mP0,
        fovyInRad,
        static_cast<REAL>(asp),
        -depth * 2.,
        -depth * 0.01);
  } else {
    Mat4_AffineTransProjectionOrtho(
        mP0,
        -view_height * asp,
        +view_height * asp,
        -view_height,
        +view_height,
        -2 * depth,
        0);
  }
  REAL mT0[16];
  {
    // the camera is placed at the origin and lookin into the -Z direction in the range [-2*depth,0]
    // to view the object we translate the object at the origin (0,0,-depth)
    const REAL t0[3] = {0.f, 0.f, -depth};
    ::delfem2::Mat4_AffineTransTranslate(mT0, t0);
  }
  const REAL mRefZ[16] = { // reflection with the XY plane
    +1.f, 0.f, 0.f, 0.f,
    0.f, +1.f, 0.f, 0.f,
    0.f, 0.f, -1.f, 0.f,
    0.f, 0.f, 0.f, +1.f};
  REAL mTmp0[16];
  ::delfem2::MatMat4(mTmp0, mT0, mP0);
  float mP1[16];
  ::delfem2::MatMat4(mP1, mTmp0, mRefZ);
  std::array<float,16> mP{};
  ::delfem2::Transpose_Mat4(mP.data(), mP1);
  return mP;
}

} // namespace delfem2

#endif  // DFM2_CAM_PROJECTION_H
