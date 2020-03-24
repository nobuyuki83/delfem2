/*
 * Copyright (c) 2019 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */

#include <cmath>
#include "delfem2/vec2.h"
#include "delfem2/vec3.h"
#include "delfem2/mat3.h"
#include "delfem2/mat4.h"
#include "delfem2/quat.h"
//
#include "delfem2/v23m3q.h"

namespace dfm2 = delfem2;

// ----------------------------------------

dfm2::CVec2d dfm2::screenXYProjection
(const CVec3d& v,
 const float* mMV,
 const float* mPj)
{
  CVec3d sp0 = screenProjection(v,mMV,mPj);
  return dfm2::CVec2d(sp0.x(),sp0.y());
}

dfm2::CVec3d dfm2::GetCartesianRotationVector
 (const CMat3d& m)
{
  const double* mat = m.mat;
  CVec3d a;
  a.p[0] = mat[7]-mat[5];
  a.p[1] = mat[2]-mat[6];
  a.p[2] = mat[3]-mat[1];
  double act = (m.Trace()-1)*0.5;
  if( act > +1 ){ act = +1; }
  if( act < -1 ){ act = -1; }
  double theta = acos(act);
  if( myIsNAN_Matrix3(theta) ){ return a; }
  if( fabs(theta) < 1.0e-5 ){ return a*0.5; }
  double mag = 0.5*theta/sin(theta);
  a *= mag;
  return a;
}

dfm2::CVec3d dfm2::GetSpinVector(const CMat3d& m)
{
  const double* mat = m.mat;
  CVec3d r;
  r.p[0] = (mat[7]-mat[5])*0.5;
  r.p[1] = (mat[2]-mat[6])*0.5;
  r.p[2] = (mat[3]-mat[1])*0.5;
  return r;
}

dfm2::CVec3d dfm2::MatVec(const CMat3d& m, const CVec3d& vec0)
{
  CVec3d vec1;
  dfm2::MatVec3(vec1.p, m.mat,vec0.p);
  return vec1;
}

dfm2::CVec3d dfm2::MatVecTrans
 (const CMat3d& m, const CVec3d& vec0)
{
  CVec3d vec1;
  MatTVec3(vec1.p, m.mat,vec0.p);
  return vec1;
}

// ---------------------------------------------------------------------

void dfm2::SetDiag(CMat3d& m, const CVec3d& d)
{
  double* mat = m.mat;
  mat[0*3+0] = d.x();
  mat[1*3+1] = d.y();
  mat[2*3+2] = d.z();
}

void dfm2::SetRotMatrix_Cartesian(CMat3d& m, const CVec3d& v)
{
//  const double vec[3] = { v.x, v.y, v.z };
  m.SetRotMatrix_Cartesian(v.p);
}

void dfm2::SetSpinTensor(CMat3d& m, const CVec3d& vec0)
{
  Mat3_Spin(m.mat, vec0.p);
  /*
  double* mat = m.mat;
  mat[0] =  0;         mat[1] = -vec0.z();   mat[2] = +vec0.y();
  mat[3] = +vec0.z();  mat[4] = 0;           mat[5] = -vec0.x();
  mat[6] = -vec0.y();  mat[7] = +vec0.x();   mat[8] = 0;
   */
}

void dfm2::SetOuterProduct
 (CMat3d& m,
  const CVec3d& vec0,
  const CVec3d& vec1 )
{
  double* mat = m.mat;
  mat[0] = vec0.x()*vec1.x(); mat[1] = vec0.x()*vec1.y(); mat[2] = vec0.x()*vec1.z();
  mat[3] = vec0.y()*vec1.x(); mat[4] = vec0.y()*vec1.y(); mat[5] = vec0.y()*vec1.z();
  mat[6] = vec0.z()*vec1.x(); mat[7] = vec0.z()*vec1.y(); mat[8] = vec0.z()*vec1.z();
}

void dfm2::SetProjection(CMat3d& m, const CVec3d& vec0)
{
  double* mat = m.mat;
  const CVec3d& u = vec0.Normalize();
  mat[0] = 1-u.x()*u.x(); mat[1] = 0-u.x()*u.y(); mat[2] = 0-u.x()*u.z();
  mat[3] = 0-u.y()*u.x(); mat[4] = 1-u.y()*u.y(); mat[5] = 0-u.y()*u.z();
  mat[6] = 0-u.z()*u.x(); mat[7] = 0-u.z()*u.y(); mat[8] = 1-u.z()*u.z();
}

// ----------------------------

dfm2::CMat3d dfm2::Mirror(const CVec3d& n)
{
  CVec3d N = n;
  N.SetNormalizedVector();
  return CMat3d::Identity() - 2*dfm2::Mat3_OuterProduct(N,N);
}

dfm2::CMat3d dfm2::Mat3_CrossCross(const CVec3d& v)
{
  return Mat3(v)*Mat3(v);
}

dfm2::CMat3d dfm2::RotMatrix_Cartesian(const CVec3d& v){
 CMat3d m;
 SetRotMatrix_Cartesian(m,v);
 return m;
}

dfm2::CMat3d dfm2::Mat3(const CVec3d& vec0){
  CMat3d m;
  SetSpinTensor(m,vec0);
  return m;
}

dfm2::CMat3d dfm2::Mat3(const CVec3d& vec0, const CVec3d& vec1){
  CMat3d m;
  SetOuterProduct(m,vec0, vec1);
  return m;
}

dfm2::CMat3d dfm2::Mat3(const CVec3d& vec0, const CVec3d& vec1, const CVec3d& vec2)
{
  CMat3d m;
  double* mat = m.mat;
  mat[0*3+0]=vec0.x(); mat[0*3+1]=vec1.x(); mat[0*3+2]=vec2.x();
  mat[1*3+0]=vec0.y(); mat[1*3+1]=vec1.y(); mat[1*3+2]=vec2.y();
  mat[2*3+0]=vec0.z(); mat[2*3+1]=vec1.z(); mat[2*3+2]=vec2.z();
  return m;
}

dfm2::CMat3d dfm2::Mat3_Spin(const CVec3d& vec0){
  CMat3d m;
  ::dfm2::Mat3_Spin(m.mat,vec0.p);
  return m;
}

dfm2::CMat3d dfm2::Mat3_OuterProduct(const CVec3d& vec0, const CVec3d& vec1 )
{
  CMat3d m;
  SetOuterProduct(m,vec0,vec1);
  return m;
}

dfm2::CMat3d dfm2::Mat3_RotCartesian(const CVec3d& vec0)
{
  CMat3d m;
  m.SetRotMatrix_Cartesian(vec0.x(), vec0.y(), vec0.z());
  return m;
}

// ------------------

namespace delfem2 {
  
CVec3d operator* (const CVec3d& v, const CMat3d& m){
  return MatVecTrans(m,v);
}
  
CVec3d operator* (const CMat3d& m, const CVec3d& v)
{
  return MatVec(m,v);
}
  
}

// ------------------------------

template <typename REAL>
dfm2::CMat3<REAL> dfm2::Mat3_MinimumRotation
(const CVec3<REAL>& V,
 const CVec3<REAL>& v)
{
  CVec3<REAL> ep = V.Normalize();
  CVec3<REAL> eq = v.Normalize();
  CVec3<REAL> n = ep^eq;
  const double st2 = n*n;
  CMat3<REAL> m;
  if( st2 < 1.0e-4f ){
    m.mat[0] = 1.f      +0.5f*(n.x()*n.x()-st2);
    m.mat[1] =    -n.z()+0.5f*(n.x()*n.y());
    m.mat[2] =    +n.y()+0.5f*(n.x()*n.z());
    m.mat[3] =    +n.z()+0.5f*(n.y()*n.x());
    m.mat[4] = 1.f      +0.5f*(n.y()*n.y()-st2);
    m.mat[5] =    -n.x()+0.5f*(n.y()*n.z());
    m.mat[6] =    -n.y()+0.5f*(n.z()*n.x());
    m.mat[7] =    +n.x()+0.5f*(n.z()*n.y());
    m.mat[8] = 1.f      +0.5f*(n.z()*n.z()-st2);
    return m;
  }
  const double st = sqrt(st2);
  const double ct = ep*eq;
  n.SetNormalizedVector();
  m.mat[0] = ct         +(1.f-ct)*n.x()*n.x();
  m.mat[1] =   -n.z()*st+(1.f-ct)*n.x()*n.y();
  m.mat[2] =   +n.y()*st+(1.f-ct)*n.x()*n.z();
  m.mat[3] =   +n.z()*st+(1.f-ct)*n.y()*n.x();
  m.mat[4] = ct         +(1.f-ct)*n.y()*n.y();
  m.mat[5] =   -n.x()*st+(1.f-ct)*n.y()*n.z();
  m.mat[6] =   -n.y()*st+(1.f-ct)*n.z()*n.x();
  m.mat[7] =   +n.x()*st+(1.f-ct)*n.z()*n.y();
  m.mat[8] = ct         +(1.f-ct)*n.z()*n.z();
  return m;
}
template dfm2::CMat3d dfm2::Mat3_MinimumRotation(const CVec3d& V, const CVec3d& v);


// --------------------------

dfm2::CMat3d dfm2::Mat3_ParallelTransport
(const CVec3d& p0,
 const CVec3d& p1,
 const CVec3d& q0,
 const CVec3d& q1)
{
  return Mat3_MinimumRotation(p1-p0, q1-q0);
}

// -----------------------------------------------------
// below: rotational inertia

// moment of inertia around origin triangle vtx (origin,d0,d1,d2) the area_density=1
// see http://www.dcs.warwick.ac.uk/~rahil/files/RigidBodySimulation.pdf
dfm2::CMat3d dfm2::Mat3_IrotTri
(const CVec3d& d0,
 const CVec3d& d1,
 const CVec3d& d2)
{
  
  CVec3d dv = d0+d1+d2;
  CMat3d I0 = Mat3_OuterProduct(d0,d0) + Mat3_OuterProduct(d1,d1) + Mat3_OuterProduct(d2,d2) + Mat3_OuterProduct(dv,dv);
  double tr0 = I0.Trace();
  CMat3d I = tr0*CMat3d::Identity()-I0;
  
  double darea = ((d1-d0)^(d2-d0)).Length();
  I *= darea/24.0;
  return I;
}

// moment of inertia triangle pyramid with vtx (origin,d0,d1,d2) volume_density = 1
// see http://www.dcs.warwick.ac.uk/~rahil/files/RigidBodySimulation.pdf
dfm2::CMat3d dfm2::Mat3_IrotTriSolid
(const CVec3d& d0,
 const CVec3d& d1,
 const CVec3d& d2)
{
  CVec3d dv = d0+d1+d2;
  CMat3d I0 = Mat3_OuterProduct(d0,d0) + Mat3_OuterProduct(d1,d1) + Mat3_OuterProduct(d2,d2) + Mat3_OuterProduct(dv,dv);
  double tr0 = I0.Trace();
  CMat3d I = tr0*CMat3d::Identity()-I0;
  
  double darea = (d0*(d1^d2));
  I *= darea/120.0;
  return I;
}

dfm2::CMat3d dfm2::Mat3_IrotLineSeg
(const CVec3d& d0,
 const CVec3d& d1)
{
  CVec3d dv = d1-d0;
  double l = dv.Length();
  CMat3d I;
  {
    I = dv.DLength()*CMat3d::Identity()-Mat3_OuterProduct(dv,dv);
    I *= l/12.0;
  }
  CVec3d p = (d0+d1)*0.5;
  I += l*(p.DLength()*CMat3d::Identity()-Mat3_OuterProduct(p,p));
  return I;
}

dfm2::CMat3d dfm2::Mat3_IrotPoint
(const CVec3d& d0)
{
  return (d0.DLength()*CMat3d::Identity()-Mat3_OuterProduct(d0,d0));
}


// above: rotational inertia
// ---------------------------------------------------------------------


void dfm2::Mat4_MatTransl(double m[16], const CMat3d& mat, const CVec3d& trans)
{
  mat.AffineMatrixTrans(m);
  m[3*4+0] = trans.x();
  m[3*4+1] = trans.y();
  m[3*4+2] = trans.z();
}


void dfm2::Mat4_ScaleMatTransl(double m[16], double scale, const CMat3d& mat, const CVec3d& trans)
{
  mat.AffineMatrixTrans(m);
  for(int i=0;i<3;++i){
  for(int j=0;j<3;++j){
    m[i*4+j] *= scale;
  }
  }
  m[3*4+0] = trans.x();
  m[3*4+1] = trans.y();
  m[3*4+2] = trans.z();
}

// ---------------------------------------------------------------------------------------

bool dfm2::isPickCircle
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

bool dfm2::isPickQuad
(const CVec3d& p0,const CVec3d& p1,const CVec3d& p2,const CVec3d& p3,
 const dfm2::CVec2d& sp, const CVec3d& pick_dir,
 const float mMV[16], const float mPj[16],
 double eps)
{
  const dfm2::CVec2d sp0 = dfm2::screenXYProjection(p0, mMV, mPj);
  const dfm2::CVec2d sp1 = dfm2::screenXYProjection(p1, mMV, mPj);
  const dfm2::CVec2d sp2 = dfm2::screenXYProjection(p2, mMV, mPj);
  const dfm2::CVec2d sp3 = dfm2::screenXYProjection(p3, mMV, mPj);
  double a01 = dfm2::Area_Tri(sp,sp0,sp1);
  double a12 = dfm2::Area_Tri(sp,sp1,sp2);
  double a23 = dfm2::Area_Tri(sp,sp2,sp3);
  double a30 = dfm2::Area_Tri(sp,sp3,sp0);
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

int dfm2::PickHandlerRotation_PosQuat
(const CVec3d& src, const CVec3d& dir,
 const CVec3d& pos, const double quat[4], double rad,
 double tol)
{
  CVec3d ax = QuatVec(quat,CVec3d(1,0,0));
  CVec3d ay = QuatVec(quat,CVec3d(0,1,0));
  CVec3d az = QuatVec(quat,CVec3d(0,0,1));
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

int dfm2::PickHandlerRotation_Mat4
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


bool dfm2::DragHandlerRot_PosQuat
(double quat[4], int ielem,
 const dfm2::CVec2d& sp0,
 const dfm2::CVec2d& sp1,
 const CVec3d& pos,
 const float mMV[16], const float mPj[16])
{
  if( ielem>=0 && ielem<3 ){
    double vi[3] = {0,0,0}; vi[ielem] = 1;
    double vo[3]; dfm2::QuatVec(vo, quat, vi);
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

bool dfm2::DragHandlerRot_Mat4
(double quat[4], int ielem,
 const dfm2::CVec2d& sp0, const dfm2::CVec2d& sp1, double mat[16],
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

bool dfm2::isPick_AxisHandler
(const dfm2::CVec2d& sp,
 const CVec3d& p,
 const CVec3d& axis,
 double len,
 const float* mMV,
 const float* mPj,
 double pick_tol)
{
  dfm2::CVec2d sp0 = dfm2::screenXYProjection(p+len*axis, mMV, mPj);
  dfm2::CVec2d sp1 = dfm2::screenXYProjection(p-len*axis, mMV, mPj);
  double sdist = GetDist_LineSeg_Point(sp, sp0, sp1);
  return sdist < pick_tol;
}

dfm2::CVec3d dfm2::drag_AxisHandler
(const dfm2::CVec2d& sp0,
 const dfm2::CVec2d& sp1,
 const CVec3d& p,
 const CVec3d& axis,
 double len,
 const float* mMV,
 const float* mPj)
{
  dfm2::CVec2d spa0 = dfm2::screenXYProjection(p+len*axis, mMV, mPj);
  dfm2::CVec2d spa1 = dfm2::screenXYProjection(p-len*axis, mMV, mPj);
  double r = (spa0-spa1)*(sp1-sp0)/(spa0-spa1).SqLength();
  return r*axis*len;
}

double dfm2::DragCircle
(const dfm2::CVec2d& sp0,
 const dfm2::CVec2d& sp1,
 const CVec3d& p,
 const CVec3d& axis,
 const float* mMV,
 const float* mPj)
{
  dfm2::CVec2d spo0 = dfm2::screenXYProjection(p, mMV, mPj);
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

bool dfm2::isPickPoint
(const dfm2::CVec2d& sp,
 const CVec3d& p,
 const float* mMV,
 const float* mPj,
 double pick_tol)
{
  dfm2::CVec2d sp0 = dfm2::screenXYProjection(p, mMV, mPj);
  return (sp - sp0).Length() < pick_tol;
}

bool dfm2::isPickCircle
(const dfm2::CVec2d& sp,
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
    dfm2::CVec2d sp0 = dfm2::screenXYProjection(p0, mMV, mPj);
    dfm2::CVec2d sp1 = dfm2::screenXYProjection(p1, mMV, mPj);
    double sdist = GetDist_LineSeg_Point(sp, sp0, sp1);
    if( sdist < pick_tol ){ return true; }
  }
  return false;
}

// ----------------------------------------------------
// quaternion

namespace delfem2 {

CVec3d operator* (const CQuatd& q, const CVec3d& v)
{
  CVec3d p;
  QuatVec(p.p, q.q, v.p);
  return p;
}

}

// ------------

dfm2::CQuatd dfm2::Quat_CartesianAngle(const CVec3d& p)
{
  CQuatd q;
  Quat_CartesianAngle(q.q, p.p);
  return q;
}
