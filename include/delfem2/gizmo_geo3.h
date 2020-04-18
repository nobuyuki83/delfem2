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

DFM2_INLINE bool DragHandlerRot_Mat4
(double quat[4], int ielem,
 const CVec2d& sp0, const CVec2d& sp1, double mat[16],
 const float mMV[16], const float mPj[16]);

DFM2_INLINE int PickHandlerRotation_Mat4
(const CVec3d& src, const CVec3d& dir,
 const double mat[16], double rad,
 double tol);


DFM2_INLINE bool isPickPoint
(const CVec2d& sp,
 const CVec3d& p,
 const float* mMV,
 const float* mPj,
 double pick_tol);

DFM2_INLINE bool isPickCircle
(const CVec3d& axis,
 const CVec3d& org,
 double rad,
 const CVec3d& src,
 const CVec3d& dir,
 double pick_tol);

DFM2_INLINE bool isPick_AxisHandler
(const CVec2d& sp,
 const CVec3d& p,
 const CVec3d& axis,
 double len,
 const float* mMV,
 const float* mPj,
 double pick_tol);


DFM2_INLINE double DragCircle
(const CVec2d& sp0,
 const CVec2d& sp1,
 const CVec3d& p,
 const CVec3d& axis,
 const float* mMV,
 const float* mPj);

DFM2_INLINE bool isPickCircle
(const CVec2d& sp,
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
DFM2_INLINE int PickHandlerRotation_PosQuat
(const CVec3<REAL>& src, const CVec3<REAL>& dir,
 const CVec3<REAL>& pos, const REAL quat[4], REAL rad,
 REAL tol);

DFM2_INLINE bool DragHandlerRot_PosQuat
(double quat[4], int ielem,
 const CVec2d& sp0, const CVec2d& sp1,
 const CVec3d& pos,
 const float mMV[16], const float mPj[16]);

DFM2_INLINE CVec3d drag_AxisHandler
(const CVec2d& sp0,
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
    size = 1.1;
    quat[0] = 1.0;
    quat[1] = 0.0;
    quat[2] = 0.0;
    quat[3] = 0.0;
    ielem_picked = -1;
    pos = CVec3<REAL>(0,0,0);
  }
  void Pick(bool is_down, const REAL src[3], const REAL dir[3], REAL tol){
    if( !is_down ){
      ielem_picked = -1;
      return;
    }
    ielem_picked = PickHandlerRotation_PosQuat(CVec3<REAL>(src),
                                               CVec3<REAL>(dir),
                                               pos, quat,size, tol);
  }
  void Drag(const REAL src0[3], const REAL src1[3], const REAL dir[3]){
    using CV3 = CVec3<REAL>;
    using CQ = CQuat<REAL>;
    int ielem = ielem_picked;
    if( ielem>=0 && ielem<3 ){
      CV3 va = (CQ(quat)*CV3::Axis(ielem)).Normalize();
      CV3 pz0,qz0; Nearest_Line_Circle(pz0,qz0,
                                       CV3(src0), CV3(dir),
                                       pos,va, size);
      CV3 pz1,qz1; Nearest_Line_Circle(pz1,qz1,
                                       CV3(src1), CV3(dir),
                                       pos,va, size);
      CV3 a0 = (qz0-pos)/size;
      CV3 a1 = (qz1-pos)/size;
      double cos01 = a0*a1;
      double sin01 = (a0^a1)*va;
      double ar = atan2(sin01, cos01);
      REAL dq[4] = {
        (REAL)cos(ar*0.5),
        (REAL)(va.x()*sin(ar*0.5)),
        (REAL)(va.y()*sin(ar*0.5)),
        (REAL)(va.z()*sin(ar*0.5)) };
      REAL qtmp[4]; QuatQuat(qtmp, dq, quat);
      Copy_Quat(quat,qtmp);
    }
  }
public:
  REAL size;
  CVec3<REAL> pos;
  REAL quat[4];
  int ielem_picked;
};

}

#ifdef DFM2_HEADER_ONLY
#  include "delfem2/gizmo_geo3.cpp"
#endif


#endif /* gizmo_geo3_h */
