/*
 * Copyright (c) 2020 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */

#include "delfem2/gizmo_geo3.h"

#include <cstdio>

#include "delfem2/vec3.h"
#include "delfem2/geoproximity3_v3.h"


DFM2_INLINE bool delfem2::isPickCircle(
    const CVec3d &axis,
    const CVec3d &org,
    double rad,
    const CVec3d &src,
    const CVec3d &dir,
    double pick_tol) {
  double t = ((org - src).dot(axis)) / (dir.dot(axis));
  CVec3d p0 = src + t * dir;
  double rad0 = (p0 - org).norm();
  return fabs(rad - rad0) < pick_tol;
}

bool delfem2::isPickQuad(
    const CVec3d &p0, const CVec3d &p1, const CVec3d &p2, const CVec3d &p3,
    const delfem2::CVec2d &sp, const CVec3d &pick_dir,
    const float mMV[16], const float mPj[16],
    double eps) {
  const CVec2d sp0 = screenXYProjection(p0, mMV, mPj);
  const CVec2d sp1 = screenXYProjection(p1, mMV, mPj);
  const CVec2d sp2 = screenXYProjection(p2, mMV, mPj);
  const CVec2d sp3 = screenXYProjection(p3, mMV, mPj);
  double a01 = Area_Tri(sp, sp0, sp1);
  double a12 = Area_Tri(sp, sp1, sp2);
  double a23 = Area_Tri(sp, sp2, sp3);
  double a30 = Area_Tri(sp, sp3, sp0);
  double a0123 = a01 + a12 + a23 + a30;
  if (fabs(a0123) < 1.0e-10) return false;
  a01 /= a0123;
  a12 /= a0123;
  a23 /= a0123;
  a30 /= a0123;
  if (a01 < eps || a12 < eps || a23 < eps || a30 < eps) { return false; }
  CVec3d n0123 = Normal(p0, p1, p2) + Normal(p1, p2, p3) + Normal(p2, p3, p0) + Normal(p3, p0, p1);
  return n0123.dot(pick_dir) <= 0;
}

template<typename REAL>
DFM2_INLINE int delfem2::PickHandlerRotation_PosQuat(
    const CVec3<REAL> &src,
    const CVec3<REAL> &dir,
    const CVec3<REAL> &pos,
    const REAL quat[4],
    REAL rad,
    REAL tol) {
  using CV3 = CVec3<REAL>;
  const CV3 ax = QuatVec(quat, CV3(1, 0, 0));
  const CV3 ay = QuatVec(quat, CV3(0, 1, 0));
  const CV3 az = QuatVec(quat, CV3(0, 0, 1));
  CV3 px, qx;
  Nearest_Line_Circle(px, qx, src, dir, pos, ax, rad);
  CV3 py, qy;
  Nearest_Line_Circle(py, qy, src, dir, pos, ay, rad);
  CV3 pz, qz;
  Nearest_Line_Circle(pz, qz, src, dir, pos, az, rad);
  REAL dx = (px - src).dot(dir);
  REAL dy = (py - src).dot(dir);
  REAL dz = (pz - src).dot(dir);
  REAL lx = (px - qx).norm();
  REAL ly = (py - qy).norm();
  REAL lz = (pz - qz).norm();
  REAL dm = (fabs(dx) + fabs(dy) + fabs(dz)) * 1000;
  //  std::cout << lx << " " << ly << " " << lz << " " << dm << std::endl;
  if (lx > tol) { dx = dm; }
  if (ly > tol) { dy = dm; }
  if (lz > tol) { dz = dm; }
  if (dx < dy && dx < dz && dx < 0.9 * dm) { return 0; }
  if (dy < dz && dy < dx && dy < 0.9 * dm) { return 1; }
  if (dz < dx && dz < dy && dz < 0.9 * dm) { return 2; }
  return -1;
}
#ifdef DFM2_STATIC_LIBRARY
template int delfem2::PickHandlerRotation_PosQuat(
    const CVec3f& src, const CVec3f& dir,
    const CVec3f& pos, const float quat[4], float rad,
    float tol);
template int delfem2::PickHandlerRotation_PosQuat(
    const CVec3d& src, const CVec3d& dir,
    const CVec3d& pos, const double quat[4], double rad,
    double tol);
#endif

// -------------------------------------------------

DFM2_INLINE int delfem2::PickHandlerRotation_Mat4(
    const CVec3d &src, const CVec3d &dir,
    const double mat[16], double rad,
    double tol) {
  CVec3d ax = Mat4Vec(mat, CVec3d(1, 0, 0));
  CVec3d ay = Mat4Vec(mat, CVec3d(0, 1, 0));
  CVec3d az = Mat4Vec(mat, CVec3d(0, 0, 1));
  CVec3d pos(mat[3], mat[7], mat[11]);
  CVec3d px, qx;
  Nearest_Line_Circle(px, qx, src, dir, pos, ax, rad);
  CVec3d py, qy;
  Nearest_Line_Circle(py, qy, src, dir, pos, ay, rad);
  CVec3d pz, qz;
  Nearest_Line_Circle(pz, qz, src, dir, pos, az, rad);
  double dx = (px - src).dot(dir);
  double dy = (py - src).dot(dir);
  double dz = (pz - src).dot(dir);
  double lx = (px - qx).norm();
  double ly = (py - qy).norm();
  double lz = (pz - qz).norm();
  double dm = (fabs(dx) + fabs(dy) + fabs(dz)) * 1000;
  if (lx > tol) { dx = dm; }
  if (ly > tol) { dy = dm; }
  if (lz > tol) { dz = dm; }
  if (dx < dy && dx < dz && dx < 0.9 * dm) { return 0; }
  if (dy < dz && dy < dx && dy < 0.9 * dm) { return 1; }
  if (dz < dx && dz < dy && dz < 0.9 * dm) { return 2; }
  return -1;
}

DFM2_INLINE bool delfem2::DragHandlerRot_PosQuat(
    double quat[4],
    int ielem,
    const CVec2d &sp0,
    const CVec2d &sp1,
    const CVec3d &pos,
    const float mMV[16],
    const float mPj[16]) {
  if (ielem >= 0 && ielem < 3) {
    double vi[3] = {0, 0, 0};
    vi[ielem] = 1;
    double vo[3];
    QuatVec(vo, quat, vi);
    CVec3d v0(0, 0, 0);
    v0[ielem] = 1;
    CVec3d v1(vo[0], vo[1], vo[2]);
    v1.normalize();
    double ar = -DragCircle(sp0, sp1, pos, v1, mMV, mPj);
    const double dq[4] = {v0.x * sin(ar * 0.5), v0.y * sin(ar * 0.5), v0.z * sin(ar * 0.5), cos(ar * 0.5)};
    double qtmp[4];
    QuatQuat(qtmp, dq, quat);
    Copy_Quat(quat, qtmp);
    return true;
  }
  return false;
}

bool delfem2::DragHandlerRot_Mat4(
    double quat[4],
    int ielem,
    const delfem2::CVec2d &sp0,
    const delfem2::CVec2d &sp1,
    double mat[16],
    const float mMV[16],
    const float mPj[16]) {
  if (ielem >= 0 && ielem < 3) {
    double vi[3] = {0, 0, 0};
    vi[ielem] = 1;
    double vo[3];
    Mat4Vec3(vo, mat, vi);
    CVec3d v0(0, 0, 0);
    v0[ielem] = 1;
    CVec3d v1(vo[0], vo[1], vo[2]);
    v1.normalize();
    CVec3d pos(mat[3], mat[7], mat[11]);
    const double ar = DragCircle(sp0, sp1, pos, v1, mMV, mPj);
    const double dq[4] = {v0.x * sin(ar * 0.5), v0.y * sin(ar * 0.5), v0.z * sin(ar * 0.5), cos(ar * 0.5)};
    double qtmp[4];
    QuatQuat(qtmp, quat, dq);
    Copy_Quat(quat, qtmp);
    return true;
  }
  return false;
}

// -------------------------------------------------

DFM2_INLINE bool delfem2::isPick_AxisHandler(
    const delfem2::CVec2d &sp,
    const CVec3d &p,
    const CVec3d &axis,
    double len,
    const float *mMV,
    const float *mPj,
    double pick_tol) {
  delfem2::CVec2d sp0 = delfem2::screenXYProjection(p + len * axis, mMV, mPj);
  delfem2::CVec2d sp1 = delfem2::screenXYProjection(p - len * axis, mMV, mPj);
  double sdist = GetDist_LineSeg_Point(sp, sp0, sp1);
  return sdist < pick_tol;
}

DFM2_INLINE delfem2::CVec3d delfem2::drag_AxisHandler(
    const CVec2d &sp0,
    const CVec2d &sp1,
    const CVec3d &p,
    const CVec3d &axis,
    double len,
    const float *mMV,
    const float *mPj) {
  CVec2d spa0 = screenXYProjection(p + len * axis, mMV, mPj);
  CVec2d spa1 = screenXYProjection(p - len * axis, mMV, mPj);
  double r = (spa0 - spa1).dot(sp1 - sp0) / (spa0 - spa1).squaredNorm();
  return r * axis * len;
}

DFM2_INLINE double delfem2::DragCircle(
    const CVec2d &sp0,
    const CVec2d &sp1,
    const CVec3d &p,
    const CVec3d &axis,
    const float *mMV,
    const float *mPj) {
  CVec2d spo0 = screenXYProjection(p, mMV, mPj);
  double area = Area_Tri(sp0, spo0, sp1);
  double angl = area / ((sp0 - spo0).norm() * (sp1 - spo0).norm());
  {
    CVec3d a3 = screenUnProjectionDirection(axis, mMV, mPj);
    if (a3.z < 0) { angl *= -1; }
  }
  return angl;
  //  CMatrix3 R; R.SetRotMatrix_Cartesian(angl*axis);
  //  return R;
}

DFM2_INLINE bool delfem2::isPickPoint(
    const CVec2d &sp,
    const CVec3d &p,
    const float *mMV,
    const float *mPj,
    double pick_tol) {
  CVec2d sp0 = screenXYProjection(p, mMV, mPj);
  return (sp - sp0).norm() < pick_tol;
}

DFM2_INLINE bool delfem2::isPickCircle(
    const CVec2d &sp,
    const CVec3d &p,
    const CVec3d &axis,
    double r,
    const float *mMV,
    const float *mPj,
    double pick_tol) {
  const int ndiv = 32;
  double rdiv = 3.1415 * 2.0 / ndiv;
  CVec3d x, y;
  GetVertical2Vector(axis, x, y);
  for (int idiv = 0; idiv < ndiv + 1; idiv++) {
    int jdiv = idiv + 1;
    CVec3d p0 = p + (r * sin(rdiv * idiv)) * x + (r * cos(rdiv * idiv)) * y;
    CVec3d p1 = p + (r * sin(rdiv * jdiv)) * x + (r * cos(rdiv * jdiv)) * y;
    CVec2d sp0 = screenXYProjection(p0, mMV, mPj);
    CVec2d sp1 = screenXYProjection(p1, mMV, mPj);
    double sdist = GetDist_LineSeg_Point(sp, sp0, sp1);
    if (sdist < pick_tol) { return true; }
  }
  return false;
}


// ================================================

template<typename REAL>
void delfem2::CGizmo_Rotation<REAL>::Pick(
    bool is_down,
    const REAL src[3],
    const REAL dir[3], REAL tol) {
  if (!is_down) {
    ielem_picked = -1;
    return;
  }
  ielem_picked = PickHandlerRotation_PosQuat(CVec3<REAL>(src),
                                             CVec3<REAL>(dir),
                                             pos, quat, size, tol);
}
#ifdef DFM2_STATIC_LIBRARY
template void delfem2::CGizmo_Rotation<float>::Pick
 (bool is_down,
  const float src[3],
  const float dir[3], float tol);
template void delfem2::CGizmo_Rotation<double>::Pick
(bool is_down,
 const double src[3],
 const double dir[3], double tol);
#endif

// --------------------------

template<typename REAL>
void delfem2::CGizmo_Rotation<REAL>::Drag(
    const REAL src0[3],
    const REAL src1[3],
    const REAL dir[3]) {
  using CV3 = CVec3<REAL>;
  using CQ = CQuat<REAL>;
  int ielem = ielem_picked;
  if (ielem >= 0 && ielem < 3) {
    CV3 va = (CQ(quat) * CV3::Axis(ielem)).normalized();
    CV3 pz0, qz0;
    Nearest_Line_Circle(pz0, qz0,
                        CV3(src0), CV3(dir),
                        pos, va, size);
    CV3 pz1, qz1;
    Nearest_Line_Circle(pz1, qz1,
                        CV3(src1), CV3(dir),
                        pos, va, size);
    CV3 a0 = (qz0 - pos) / size;
    CV3 a1 = (qz1 - pos) / size;
    const double ar = atan2((a0 ^ a1).dot(va), a0.dot(a1));
    const REAL dq[4] = {
        (REAL) (va.x * sin(ar * 0.5)),
        (REAL) (va.y * sin(ar * 0.5)),
        (REAL) (va.z * sin(ar * 0.5)),
        (REAL) cos(ar * 0.5)};
    REAL qtmp[4];
    QuatQuat(qtmp, dq, quat);
    Copy_Quat(quat, qtmp);
  }
}
#ifdef DFM2_STATIC_LIBRARY
template void delfem2::CGizmo_Rotation<float>::Drag
(const float src0[3],
 const float src1[3],
 const float dir[3]);
template void delfem2::CGizmo_Rotation<double>::Drag
 (const double src0[3],
  const double src1[3],
  const double dir[3]);
#endif

// ------------------------------------------------

template<typename REAL>
void delfem2::CGizmo_Transl<REAL>::Pick(
    bool is_down,
    const REAL src[3],
    const REAL dir[3],
    REAL tol) {
  if (!is_down) {
    ielem_picked = -1;
    return;
  }
  {
    CVec3<REAL> pls, pl;
    ::delfem2::nearest_LineSeg_Line(pls, pl,
                                    pos - size * CVec3<REAL>::Axis(0),
                                    pos + size * CVec3<REAL>::Axis(0),
                                    CVec3<REAL>(src), CVec3<REAL>(dir));
    if ((pls - pl).norm() < tol) {
      ielem_picked = 0;
      return;
    }
  }
  {
    CVec3<REAL> pls, pl;
    ::delfem2::nearest_LineSeg_Line(pls, pl,
                                    pos - size * CVec3<REAL>::Axis(1),
                                    pos + size * CVec3<REAL>::Axis(1),
                                    CVec3<REAL>(src), CVec3<REAL>(dir));
    if ((pls - pl).norm() < tol) {
      ielem_picked = 1;
      return;
    }
  }
  {
    CVec3<REAL> pls, pl;
    ::delfem2::nearest_LineSeg_Line(pls, pl,
                                    pos - size * CVec3<REAL>::Axis(2),
                                    pos + size * CVec3<REAL>::Axis(2),
                                    CVec3<REAL>(src), CVec3<REAL>(dir));
    if ((pls - pl).norm() < tol) {
      ielem_picked = 2;
      return;
    }
  }
  ielem_picked = -1;
}
#ifdef DFM2_STATIC_LIBRARY
template void delfem2::CGizmo_Transl<float>::Pick
(bool is_down,
 const float src[3],
 const float dir[3],
 float tol);
template void delfem2::CGizmo_Transl<double>::Pick
(bool is_down,
 const double src[3],
 const double dir[3],
 double tol);
#endif

// -----------------

template<typename REAL>
void delfem2::CGizmo_Transl<REAL>::Drag(
    const REAL src0[3],
    const REAL src1[3],
    const REAL dir[3]) {
  if (this->ielem_picked < 0 || this->ielem_picked >= 3) { return; }
  using CV3 = CVec3<REAL>;
  REAL D0;
  CV3 Da0, Db0;
  nearest_Line_Line(D0, Da0, Db0,
                    this->pos, CV3::Axis(ielem_picked),
                    CV3(src0), CV3(dir));
  REAL D1;
  CV3 Da1, Db1;
  nearest_Line_Line(D1, Da1, Db1,
                    this->pos, CV3::Axis(ielem_picked),
                    CV3(src1), CV3(dir));
  if (fabs(D0) > 1.0e-10) {
    pos += (Da1 - Da0) / D0;
  }
}
#ifdef DFM2_STATIC_LIBRARY
template void delfem2::CGizmo_Transl<float>::Drag(const float src0[3],
                                                  const float src1[3],
                                                  const float dir[3]);
template void delfem2::CGizmo_Transl<double>::Drag(const double src0[3],
                                                   const double src1[3],
                                                   const double dir[3]);
#endif
