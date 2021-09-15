//
//  gizmo_geo3.h
//  examples_glfwold_hdronly
//
//  Created by Nobuyuki Umetani on 2020-04-17.
//

#ifndef DFM2_GIZMO_GEO3_H
#define DFM2_GIZMO_GEO3_H

#include "delfem2/dfm2_inline.h"
#include "delfem2/geo3_v23m34q.h"

namespace delfem2 {

DFM2_INLINE bool DragHandlerRot_Mat4(
    double quat[4],
    int ielem,
    const CVec2d& sp0,
    const CVec2d& sp1,
    double mat[16],
    const float mMV[16],
    const float mPj[16]);

DFM2_INLINE int PickHandlerRotation_Mat4(
    const CVec3d& src,
    const CVec3d& dir,
    const double mat[16],
    double rad,
    double tol);


DFM2_INLINE bool isPickPoint(
    const CVec2d& sp,
    const CVec3d& p,
    const float* mMV,
    const float* mPj,
    double pick_tol);

DFM2_INLINE bool isPickCircle(
    const CVec3d& axis,
    const CVec3d& org,
    double rad,
    const CVec3d& src,
    const CVec3d& dir,
    double pick_tol);

DFM2_INLINE bool isPick_AxisHandler(
    const CVec2d& sp,
    const CVec3d& p,
    const CVec3d& axis,
    double len,
    const float* mMV,
    const float* mPj,
    double pick_tol);

DFM2_INLINE double DragCircle(
    const CVec2d& sp0,
    const CVec2d& sp1,
    const CVec3d& p,
    const CVec3d& axis,
    const float* mMV,
    const float* mPj);

DFM2_INLINE bool isPickCircle(
    const CVec2d& sp,
    const CVec3d& p,
    const CVec3d& axis,
    double r,
    const float* mMV,
    const float* mPj,
    double pick_tol);

/**
 * @details defiend for float and double for static library
 */
template <typename REAL>
DFM2_INLINE int PickHandlerRotation_PosQuat(
    const CVec3<REAL>& src,
    const CVec3<REAL>& dir,
    const CVec3<REAL>& pos,
    const REAL quat[4],
    REAL rad,
    REAL tol);

DFM2_INLINE bool DragHandlerRot_PosQuat(
    double quat[4],
    int ielem,
    const CVec2d& sp0,
    const CVec2d& sp1,
    const CVec3d& pos,
    const float mMV[16],
    const float mPj[16]);

DFM2_INLINE CVec3d drag_AxisHandler(
    const CVec2d& sp0,
    const CVec2d& sp1,
    const CVec3d& p,
    const CVec3d& axis,
    double len,
    const float* mMV,
    const float* mPj);

template <typename REAL>
class CGizmo_Rotation{
public:
  CGizmo_Rotation(){
    size = 1.1f;
    quat[0] = 0.0;
    quat[1] = 0.0;
    quat[2] = 0.0;
    quat[3] = 1.0;
    ielem_picked = -1;
    pos = CVec3<REAL>(0,0,0);
  }
  void Pick(bool is_down, const REAL src[3], const REAL dir[3], REAL tol);
  void Drag(const REAL src0[3], const REAL src1[3], const REAL dir[3]);
public:
  REAL size;
  CVec3<REAL> pos;
  REAL quat[4];
  int ielem_picked;
};

template <typename REAL>
class CGizmo_Transl{
public:
  CGizmo_Transl(){
    size = 1.1f;
    ielem_picked = -1;
    pos = CVec3<REAL>(0,0,0);
  }
  void Pick(bool is_down, const REAL src[3], const REAL dir[3], REAL tol);
  void Drag(const REAL src0[3], const REAL src1[3], const REAL dir[3]);
public:
  REAL size;
  CVec3<REAL> pos;
  int ielem_picked;
};

template <typename REAL>
class CGizmo_Affine {
public:
  void SetPivot(
      const CVec3<REAL>& pivot){
    pivot0 = pivot;
  }
  void Pick(const float src[3], const float dir[3]){
    if( this->igizmo_mode == 1 ){
      gizmo_rot.Pick(true, src, dir, 0.1f);
    }
    else if( this->igizmo_mode == 0 ){
      gizmo_trnsl.Pick(true, src, dir, 0.1f);
    }
  }
  void Drag(const float src0[3], const float src1[3], const float dir[3])
  {
    if( this->igizmo_mode == 1 ){
      gizmo_rot.Drag(src0,src1,dir);
    }
    else if( this->igizmo_mode == 0 ){
      gizmo_trnsl.Drag(src0,src1,dir);
      gizmo_rot.pos = gizmo_trnsl.pos;
    }
  }
  CMat4<REAL> Affine() const {
    CMat4<REAL> m0 = delfem2::CMat4<REAL>::Translate((-pivot0).p);
    CMat4<REAL> m1 = delfem2::CMat4<REAL>::Quat(gizmo_rot.quat);
    CMat4<REAL> m2 = delfem2::CMat4<REAL>::Translate(gizmo_trnsl.pos.p);
    return m2*m1*m0;
  }
public:
  int igizmo_mode = -1;
  CGizmo_Transl<REAL> gizmo_trnsl;
  CGizmo_Rotation<REAL> gizmo_rot;
  CVec3<REAL> pivot0;
};

}

#ifndef DFM2_STATIC_LIBRARY
#  include "delfem2/gizmo_geo3.cpp"
#endif


#endif /* gizmo_geo3_h */
