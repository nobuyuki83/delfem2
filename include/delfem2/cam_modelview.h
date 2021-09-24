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

#ifndef DFM2_CAM_MODELVIEW_H
#define DFM2_CAM_MODELVIEW_H

#include <iostream>
#include <cmath>
#include <cstdio>  // memcpy
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
class ModelView {
 public:
  virtual ~ModelView() = default;
  /**
 * @brief Compute "projection matrix" for OpenGL.
 * @detial OpenGL will draw the object in the cube [-1,+1, -1,+1, -1,+1] looking from -Z
 * A 3D point is transfromed with this affine matrix and then a cube [-1,+1, -1,+1, -1,+1] is looked from -Z directoin.
 * The projection matrix will mirror the object Z.
 * To look from +Z direction, The transformation needs a mirror transformation in XY plane.
 * @param mP
 * @param asp
 */
  [[nodiscard]] virtual std::array<float, 16> GetMatrix() const = 0;
};

// ----------------------------------------------------

class ModelView_Trackball : public ModelView {
 public:
  /**
   *
   * @tparam REAL
   * @param mMV model view matrix (column major order)
   * @detail column major
   */
  [[nodiscard]] std::array<float, 16> GetMatrix() const override {
    const CMat4f mT = CMat4f::Translate(-anchor[0], -anchor[1], -anchor[2]);
    const CMat4f mR = CMat4f::Quat(quaternion);
    std::array<float, 16> mMV{};
    (mR*mT).CopyTo(mMV.data());
    return mMV;
  }
  void Rot_Camera(float dx, float dy) {
    const float a = sqrt(dx * dx + dy * dy);
    const float ar = a * 0.5f; // angle
    const float dq[4] = {
        -dy * std::sin(ar) / a,
        dx * std::sin(ar) / a,
        0.f,
        std::cos(ar)};
    if (a != 0.0) {
      float qtmp[4];
      QuatQuat(qtmp, dq, quaternion);
      Copy_Quat(quaternion, qtmp);
    }
  }
  float anchor[3] = {0, 0, 0};
  float quaternion[4] = {0, 0, 0, 1};
};

// ------------------------------------

template<typename REAL>
class ModelView_Ytop {
 public:
  /**
   *
   * @tparam REAL
   * @param mMV model view matrix (column major order)
   * @detail column major
   */
  void Mat4ColumnMajor(float mMV[16]) const {
    REAL Mt[16];
    {
      const REAL transd[3] = {trans[0], trans[1], trans[2]};
      Mat4_AffineTransTranslate(Mt, transd);
    }
    REAL Mr[16];
    {
      REAL x = std::sin(theta);
      REAL z = std::cos(theta);
      REAL y = std::sin(psi);
      x *= std::cos(psi);
      z *= std::cos(psi);
      Mat4_AffineTransLookAt(Mr, x, y, z, 0., 0., 0., 0., 1., 0.);
    }
    MatMat4(mMV, Mr, Mt);
  }
  void Rot_Camera(REAL dx, REAL dy) {
    theta -= dx;
    psi -= dy;
  }
  void Pan_Camera(REAL dx, REAL dy, REAL s) {
    // double s = view_height / scale;
    trans[0] += s * dx;
    trans[1] += s * dy;
    trans[2] += s * 0.0;
  }
 public:
  double trans[3] = {0, 0, 0};
  REAL theta = 0;
  REAL psi = 0;
};

// ------------------------------------

template<typename REAL>
class ModelView_Ztop {
 public:
  /**
   *
   * @tparam REAL
   * @param mMV model view matrix (column major order)
   * @detail column major
   */
  void Mat4ColumnMajor(float mMV[16]) const {
    REAL Mt[16];
    {
      const REAL transd[3] = {trans[0], trans[1], trans[2]};
      Mat4_AffineTransTranslate(Mt, transd);
    }
    REAL Mr[16];
    {
      REAL x = std::sin(theta);
      REAL y = std::cos(theta);
      REAL z = std::sin(psi);
      x *= std::cos(psi);
      y *= std::cos(psi);
      Mat4_AffineTransLookAt(Mr, x, y, z, 0., 0., 0., 0., 0., 1.);
    }
    MatMat4(mMV, Mr, Mt);
  }
  void Rot_Camera(REAL dx, REAL dy) {
    theta += dx;
    psi -= dy;
  }
  void Pan_Camera(REAL dx, REAL dy, REAL s) {
    // double s = view_height / scale;
    trans[0] += s * dx;
    trans[1] += s * dy;
    trans[2] += s * 0.0;
  }
 public:
  double trans[3] = {0, 0, 0};
  REAL theta = 0;
  REAL psi = 0;
};

} // namespace delfem2

#endif  // DFM2_CAM_MODELVIEW_H
