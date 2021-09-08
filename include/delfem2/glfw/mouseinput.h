#ifndef DFM2_MOUSEINPUT_H
#define DFM2_MOUSEINPUT_H

#include "delfem2/mat4.h"

namespace delfem2 {

class CMouseInput {
public:
  CMouseInput() {
    ibutton = -1;
  }

  void MouseRay(
      float src[3],
      float dir[3],
      [[maybe_unused]] float asp,
      const float mMVP[16]) const {
    float mMVP_inv[16];
    ::delfem2::Inverse_Mat4(mMVP_inv, mMVP);
    const float ps[3] = {(float) mouse_x, (float) mouse_y, -1.0};
    const float pe[3] = {(float) mouse_x, (float) mouse_y, +1.0};
    float qs[3];
    ::delfem2::Vec3_Vec3Mat4_AffineProjection(qs, ps, mMVP_inv);
    float qe[3];
    ::delfem2::Vec3_Vec3Mat4_AffineProjection(qe, pe, mMVP_inv);
    src[0] = qs[0];
    src[1] = qs[1];
    src[2] = qs[2];
    dir[0] = qe[0] - qs[0];
    dir[1] = qe[1] - qs[1];
    dir[2] = qe[2] - qs[2];
  }

  void RayMouseMove(
      float src0[3],
      float src1[3],
      float dir0[3],
      float dir1[3],
      [[maybe_unused]] float asp,
      const float mMVP[16]) {
    float mMVP_inv[16];
    ::delfem2::Inverse_Mat4(mMVP_inv, mMVP);
    const float p0s[3] = {(float) (mouse_x - dx), (float) (mouse_y - dy), -1.f};
    const float p0e[3] = {(float) (mouse_x - dx), (float) (mouse_y - dy), +1.f};
    const float p1s[3] = {(float) mouse_x, (float) mouse_y, -1.f};
    const float p1e[3] = {(float) mouse_x, (float) mouse_y, +1.f};
    float q0s[3];
    ::delfem2::Vec3_Vec3Mat4_AffineProjection(q0s, p0s, mMVP_inv);
    float q0e[3];
    ::delfem2::Vec3_Vec3Mat4_AffineProjection(q0e, p0e, mMVP_inv);
    float q1s[3];
    ::delfem2::Vec3_Vec3Mat4_AffineProjection(q1s, p1s, mMVP_inv);
    float q1e[3];
    ::delfem2::Vec3_Vec3Mat4_AffineProjection(q1e, p1e, mMVP_inv);
    //
    src0[0] = q0s[0];
    src0[1] = q0s[1];
    src0[2] = q0s[2];
    //
    src1[0] = q1s[0];
    src1[1] = q1s[1];
    src1[2] = q1s[2];
    //
    dir0[0] = q0e[0] - q0s[0];
    dir0[1] = q0e[1] - q0s[1];
    dir0[2] = q0e[2] - q0s[2];
    //
    dir1[0] = q1e[0] - q1s[0];
    dir1[1] = q1e[1] - q1s[1];
    dir1[2] = q1e[2] - q1s[2];
  }

public:
  int imodifier;
  int ibutton;
  double mouse_x, mouse_y;
  double dx, dy;
  double mouse_x_down, mouse_y_down;
};


}

#endif