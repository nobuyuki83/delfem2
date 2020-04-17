//
//  gizmo_geo3.h
//  examples_glfwold_hdronly
//
//  Created by Nobuyuki Umetani on 2020-04-17.
//

#ifndef DFM2_GIZMO_GEO3_H
#define DFM2_GIZMO_GEO3_H

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

DFM2_INLINE int PickHandlerRotation_PosQuat
(const CVec3d& src, const CVec3d& dir,
 const CVec3d& pos, const double quat[4], double rad,
 double tol);

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
    pos = CVec3d(0,0,0);
  }
  void Pick(bool is_down, const double src[3], const double dir[3], double tol){
    if( !is_down ){
      ielem_picked = -1;
      return;
    }
    ielem_picked = PickHandlerRotation_PosQuat(src,dir, CVec3d(0,0,0), quat,size, tol);
  }
  void Drag(const float src0[3], const float src1[3], const float dir[3]){
    int ielem = ielem_picked;
    if( ielem>=0 && ielem<3 ){
      CVec3d org_circle(0,0,0);
      CVec3d va = (CQuatd(quat)*CVec3d::Axis(ielem)).Normalize();
      CVec3d pz0,qz0; Nearest_Line_Circle(pz0,qz0,
                                          CVec3d(src0), CVec3d(dir),
                                          org_circle,va, size);
      CVec3d pz1,qz1; Nearest_Line_Circle(pz1,qz1,
                                          CVec3d(src1), CVec3d(dir),
                                          org_circle,va, size);
      CVec3d a0 = (qz0-org_circle)/size;
      CVec3d a1 = (qz1-org_circle)/size;
      double cos01 = a0*a1;
      double sin01 = (a0^a1)*va;
      double ar = atan2(sin01, cos01);
      double dq[4] = { cos(ar*0.5), va.x()*sin(ar*0.5), va.y()*sin(ar*0.5), va.z()*sin(ar*0.5) };
      double qtmp[4]; QuatQuat(qtmp, dq, quat);
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




DFM2_INLINE bool delfem2::isPickCircle
 (const CVec3d& axis,
  const CVec3d& org,
  double rad,
  const CVec3d& src,
  const CVec3d& dir,
  double pick_tol)
{
  double t = ((org-src)*axis)/(dir*axis);
  CVec3d p0 = src+t*dir;
  double rad0 = (p0-org).Length();
  return fabs(rad - rad0) < pick_tol;
}

bool delfem2::isPickQuad
 (const CVec3d& p0,const CVec3d& p1,const CVec3d& p2,const CVec3d& p3,
  const delfem2::CVec2d& sp, const CVec3d& pick_dir,
  const float mMV[16], const float mPj[16],
  double eps)
{
  const CVec2d sp0 = screenXYProjection(p0, mMV, mPj);
  const CVec2d sp1 = screenXYProjection(p1, mMV, mPj);
  const CVec2d sp2 = screenXYProjection(p2, mMV, mPj);
  const CVec2d sp3 = screenXYProjection(p3, mMV, mPj);
  double a01 = Area_Tri(sp,sp0,sp1);
  double a12 = Area_Tri(sp,sp1,sp2);
  double a23 = Area_Tri(sp,sp2,sp3);
  double a30 = Area_Tri(sp,sp3,sp0);
  double a0123 = a01+a12+a23+a30;
  if( fabs(a0123) < 1.0e-10 ) return false;
  a01 /= a0123;
  a12 /= a0123;
  a23 /= a0123;
  a30 /= a0123;
  if( a01<eps || a12<eps || a23<eps || a30<eps ){ return false; }
  CVec3d n0123 = Normal(p0,p1,p2) + Normal(p1,p2,p3) + Normal(p2,p3,p0) + Normal(p3,p0,p1);
  return n0123 * pick_dir <= 0;
}

DFM2_INLINE int delfem2::PickHandlerRotation_PosQuat
 (const CVec3d& src, const CVec3d& dir,
  const CVec3d& pos, const double quat[4], double rad,
  double tol)
{
  const CVec3d ax = QuatVec(quat,CVec3d(1,0,0));
  const CVec3d ay = QuatVec(quat,CVec3d(0,1,0));
  const CVec3d az = QuatVec(quat,CVec3d(0,0,1));
  CVec3d px,qx; Nearest_Line_Circle(px,qx, src,dir, pos,ax, rad);
  CVec3d py,qy; Nearest_Line_Circle(py,qy, src,dir, pos,ay, rad);
  CVec3d pz,qz; Nearest_Line_Circle(pz,qz, src,dir, pos,az, rad);
  double dx = (px-src)*dir;
  double dy = (py-src)*dir;
  double dz = (pz-src)*dir;
  double lx = (px-qx).Length();
  double ly = (py-qy).Length();
  double lz = (pz-qz).Length();
  double dm = (fabs(dx)+fabs(dy)+fabs(dz))*1000;
  std::cout << lx << " " << ly << " " << lz << " " << dm << std::endl;
  if( lx>tol ){ dx = dm; }
  if( ly>tol ){ dy = dm; }
  if( lz>tol ){ dz = dm; }
  if( dx < dy && dx < dz  && dx < 0.9*dm ){ return 0; }
  if( dy < dz && dy < dx  && dy < 0.9*dm ){ return 1; }
  if( dz < dx && dz < dy  && dz < 0.9*dm ){ return 2; }
  return -1;
}

DFM2_INLINE int delfem2::PickHandlerRotation_Mat4
 (const CVec3d& src, const CVec3d& dir,
  const double mat[16], double rad,
  double tol)
{
  CVec3d ax = Mat4Vec(mat,CVec3d(1,0,0));
  CVec3d ay = Mat4Vec(mat,CVec3d(0,1,0));
  CVec3d az = Mat4Vec(mat,CVec3d(0,0,1));
  CVec3d pos(mat[3],mat[7],mat[11]);
  CVec3d px,qx; Nearest_Line_Circle(px,qx, src,dir, pos,ax, rad);
  CVec3d py,qy; Nearest_Line_Circle(py,qy, src,dir, pos,ay, rad);
  CVec3d pz,qz; Nearest_Line_Circle(pz,qz, src,dir, pos,az, rad);
  double dx = (px-src)*dir;
  double dy = (py-src)*dir;
  double dz = (pz-src)*dir;
  double lx = (px-qx).Length();
  double ly = (py-qy).Length();
  double lz = (pz-qz).Length();
  double dm = (fabs(dx)+fabs(dy)+fabs(dz))*1000;
  if( lx>tol ){ dx = dm; }
  if( ly>tol ){ dy = dm; }
  if( lz>tol ){ dz = dm; }
  if( dx < dy && dx < dz  && dx < 0.9*dm ){ return 0; }
  if( dy < dz && dy < dx  && dy < 0.9*dm ){ return 1; }
  if( dz < dx && dz < dy  && dz < 0.9*dm ){ return 2; }
  return -1;
}


DFM2_INLINE bool delfem2::DragHandlerRot_PosQuat
 (double quat[4], int ielem,
  const CVec2d& sp0,
  const CVec2d& sp1,
  const CVec3d& pos,
  const float mMV[16], const float mPj[16])
{
  if( ielem>=0 && ielem<3 ){
    double vi[3] = {0,0,0}; vi[ielem] = 1;
    double vo[3]; QuatVec(vo, quat, vi);
    CVec3d v0(0,0,0); v0[ielem] = 1;
    CVec3d v1(vo[0],vo[1],vo[2]); v1.SetNormalizedVector();
    double ar = -DragCircle(sp0,sp1, pos, v1, mMV, mPj);
    double dq[4] = { cos(ar*0.5), v0.x()*sin(ar*0.5), v0.y()*sin(ar*0.5), v0.z()*sin(ar*0.5) };
    double qtmp[4]; QuatQuat(qtmp, dq, quat);
    Copy_Quat(quat,qtmp);
    return true;
  }
  return false;
}

bool delfem2::DragHandlerRot_Mat4
 (double quat[4], int ielem,
  const delfem2::CVec2d& sp0, const delfem2::CVec2d& sp1, double mat[16],
  const float mMV[16], const float mPj[16])
{
  if( ielem>=0 && ielem<3 ){
    double vi[3] = {0,0,0}; vi[ielem] = 1;
    double vo[3]; Mat4Vec3(vo, mat, vi);
    CVec3d v0(0,0,0); v0[ielem] = 1;
    CVec3d v1(vo[0],vo[1],vo[2]); v1.SetNormalizedVector();
    CVec3d pos(mat[3],mat[7],mat[11]);
    const double ar = DragCircle(sp0,sp1, pos, v1, mMV, mPj);
    const double dq[4] = { cos(ar*0.5), v0.x()*sin(ar*0.5), v0.y()*sin(ar*0.5), v0.z()*sin(ar*0.5) };
    double qtmp[4]; QuatQuat(qtmp, quat, dq);
    Copy_Quat(quat,qtmp);
    return true;
  }
  return false;
}

DFM2_INLINE bool delfem2::isPick_AxisHandler
 (const delfem2::CVec2d& sp,
  const CVec3d& p,
  const CVec3d& axis,
  double len,
  const float* mMV,
  const float* mPj,
  double pick_tol)
{
  delfem2::CVec2d sp0 = delfem2::screenXYProjection(p+len*axis, mMV, mPj);
  delfem2::CVec2d sp1 = delfem2::screenXYProjection(p-len*axis, mMV, mPj);
  double sdist = GetDist_LineSeg_Point(sp, sp0, sp1);
  return sdist < pick_tol;
}

DFM2_INLINE delfem2::CVec3d delfem2::drag_AxisHandler
 (const CVec2d& sp0,
  const CVec2d& sp1,
  const CVec3d& p,
  const CVec3d& axis,
  double len,
  const float* mMV,
  const float* mPj)
{
  CVec2d spa0 = screenXYProjection(p+len*axis, mMV, mPj);
  CVec2d spa1 = screenXYProjection(p-len*axis, mMV, mPj);
  double r = (spa0-spa1)*(sp1-sp0)/(spa0-spa1).SqLength();
  return r*axis*len;
}

DFM2_INLINE double delfem2::DragCircle
 (const CVec2d& sp0,
  const CVec2d& sp1,
  const CVec3d& p,
  const CVec3d& axis,
  const float* mMV,
  const float* mPj)
{
  CVec2d spo0 = screenXYProjection(p, mMV, mPj);
  double area = Area_Tri(sp0, spo0, sp1);
  double angl = area / ( (sp0-spo0).Length() * (sp1-spo0).Length() );
  {
    CVec3d a3 = screenUnProjectionDirection(axis,mMV,mPj);
    if( a3.z() < 0 ){ angl *= -1; }
  }
  return angl;
  //  CMatrix3 R; R.SetRotMatrix_Cartesian(angl*axis);
  //  return R;
}

DFM2_INLINE bool delfem2::isPickPoint
 (const CVec2d& sp,
  const CVec3d& p,
  const float* mMV,
  const float* mPj,
  double pick_tol)
{
  CVec2d sp0 = screenXYProjection(p, mMV, mPj);
  return (sp - sp0).Length() < pick_tol;
}

DFM2_INLINE bool delfem2::isPickCircle
 (const CVec2d& sp,
  const CVec3d& p,
  const CVec3d& axis,
  double r,
  const float* mMV,
  const float* mPj,
  double pick_tol)
{
  const int ndiv = 32;
  double rdiv = 3.1415*2.0/ndiv;
  CVec3d x,y; GetVertical2Vector(axis, x, y);
  for(int idiv=0;idiv<ndiv+1;idiv++){
    int jdiv = idiv+1;
    CVec3d p0 = p+(r*sin(rdiv*idiv))*x+(r*cos(rdiv*idiv))*y;
    CVec3d p1 = p+(r*sin(rdiv*jdiv))*x+(r*cos(rdiv*jdiv))*y;
    CVec2d sp0 = screenXYProjection(p0, mMV, mPj);
    CVec2d sp1 = screenXYProjection(p1, mMV, mPj);
    double sdist = GetDist_LineSeg_Point(sp, sp0, sp1);
    if( sdist < pick_tol ){ return true; }
  }
  return false;
}



#endif /* gizmo_geo3_h */
