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
#include "delfem2/dfm2_inline.h"

// ---------------------------------------------------

namespace delfem2 {

/**
 * Pure virtual class
 * @tparam REAL
 */
class Projection {
 public:
  virtual ~Projection() = default;
  /**
 * @brief Compute "projection matrix" for OpenGL.
 * @detial OpenGL will draw the object in the cube [-1,+1, -1,+1, -1,+1] looking from -Z
 * A 3D point is transfromed with this affine matrix and then a cube [-1,+1, -1,+1, -1,+1]
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
 */
class Projection_LookOriginFromZplus : public Projection{
 public:
  explicit Projection_LookOriginFromZplus(float view_height = 1,
                                 bool is_pars = false,
                                 float fovy = 10)
  : view_height(view_height), is_pars(is_pars), fovy(fovy){}
  
  ~Projection_LookOriginFromZplus() override = default;
  /**
   * @brief Compute "projection matrix" for OpenGL.
   * @detial OpenGL will draw the object in the cube [-1,+1, -1,+1, -1,+1] looking from -Z
   * A 3D point is transfromed with this affine matrix and then a cube [-1,+1, -1,+1, -1,+1] is looked from -Z directoin.
   * The projection matrix will mirror the object Z.
   * To look from +Z direction, The transformation needs a mirror transformation in XY plane.
   * @param mP
   * @param asp
   */
   [[nodiscard]] DFM2_INLINE std::array<float,16> GetMatrix(float asp) const override{

    const float fovyInRad = fovy * (2.f * static_cast<float>(M_PI)) / 360.f;
    const float depth = view_height / tan(fovyInRad * 0.5f);
    CMat4f mP1;
    if (is_pars) {
      Mat4_AffineProjectionFrustum(
          mP1.data(),
          fovyInRad,
          asp,
          -depth * 2.f,
          -depth * 0.01f);
    } else {
      Mat4_AffineProjectionOrtho(
          mP1.data(),
          -view_height * asp,
          +view_height * asp,
          -view_height,
          +view_height,
          -2.f * depth,
          0.f);
    }
    // the camera is placed at the origin and lookin into the -Z direction in the range [-2*depth,0]
    // to view the object we translate the object at the origin (0,0,-depth)
    CMat4f mT1;
    ::delfem2::Mat4_AffineTranslation(
        mT1.data(),
        std::array<float,3>{0.f, 0.f, -depth}.data());
    return (mP1 * mT1).GetStlArray();
   }
public:
  float view_height = 1;
  bool is_pars = false;
  float fovy = 10;
};

} // namespace delfem2

#endif  // DFM2_CAM_PROJECTION_H
