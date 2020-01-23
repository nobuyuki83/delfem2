/*
 * Copyright (c) 2019 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */

#include "delfem2/vec2.h"
#include "delfem2/vec3.h"
#include "delfem2/mat3.h"
#include "delfem2/quat.h"
//
#include "delfem2/v23m3q.h"

namespace dfm2 = delfem2;

// ----------------------------------------

dfm2::CVec2 dfm2::screenXYProjection
(const CVec3& v,
 const float* mMV,
 const float* mPj)
{
  CVec3 sp0 = screenProjection(v,mMV,mPj);
  return dfm2::CVec2(sp0.x(),sp0.y());
}

dfm2::CVec3 dfm2::GetCartesianRotationVector
 (const CMat3d& m)
{
  const double* mat = m.mat;
  CVec3 a;
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

dfm2::CVec3 dfm2::GetSpinVector(const CMat3d& m)
{
  const double* mat = m.mat;
  CVec3 r;
  r.p[0] = (mat[7]-mat[5])*0.5;
  r.p[1] = (mat[2]-mat[6])*0.5;
  r.p[2] = (mat[3]-mat[1])*0.5;
  return r;
}

dfm2::CVec3 dfm2::MatVec(const CMat3d& m, const CVec3& vec0)
{
  CVec3 vec1;
  dfm2::MatVec3(vec1.p, m.mat,vec0.p);
  return vec1;
    //  const double* mat = m.mat;
    //  vec1.p[0] = mat[0]*vec0.x + mat[1]*vec0.y + mat[2]*vec0.z;
    //  vec1.p[1] = mat[3]*vec0.x + mat[4]*vec0.y + mat[5]*vec0.z;
    //  vec1.p[2] = mat[6]*vec0.x + mat[7]*vec0.y + mat[8]*vec0.z;
}

dfm2::CVec3 dfm2::MatVecTrans
 (const CMat3d& m, const CVec3& vec0)
{
  CVec3 vec1;
  MatTransVec3(vec1.p, m.mat,vec0.p);
//  const double* mat = m.mat;
//  vec1.x = mat[0]*vec0.x + mat[3]*vec0.y + mat[6]*vec0.z;
//  vec1.y = mat[1]*vec0.x + mat[4]*vec0.y + mat[7]*vec0.z;
//  vec1.z = mat[2]*vec0.x + mat[5]*vec0.y + mat[8]*vec0.z;
  return vec1;
}

// ---------------------------------------------------------------------

void dfm2::SetDiag(CMat3d& m, const CVec3& d)
{
  double* mat = m.mat;
  mat[0*3+0] = d.x();
  mat[1*3+1] = d.y();
  mat[2*3+2] = d.z();
}

void dfm2::SetRotMatrix_Cartesian(CMat3d& m, const CVec3& v)
{
//  const double vec[3] = { v.x, v.y, v.z };
  m.SetRotMatrix_Cartesian(v.p);
}

void dfm2::SetSpinTensor(CMat3d& m, const CVec3& vec0)
{
  double* mat = m.mat;
  mat[0] =  0;       mat[1] = -vec0.z();   mat[2] = +vec0.y();
  mat[3] = +vec0.z();  mat[4] = 0;         mat[5] = -vec0.x();
  mat[6] = -vec0.y();  mat[7] = +vec0.x();   mat[8] = 0;
}

void dfm2::SetOuterProduct
 (CMat3d& m,
  const CVec3& vec0,
  const CVec3& vec1 )
{
  double* mat = m.mat;
  mat[0] = vec0.x()*vec1.x(); mat[1] = vec0.x()*vec1.y(); mat[2] = vec0.x()*vec1.z();
  mat[3] = vec0.y()*vec1.x(); mat[4] = vec0.y()*vec1.y(); mat[5] = vec0.y()*vec1.z();
  mat[6] = vec0.z()*vec1.x(); mat[7] = vec0.z()*vec1.y(); mat[8] = vec0.z()*vec1.z();
}

void dfm2::SetProjection(CMat3d& m, const CVec3& vec0)
{
  double* mat = m.mat;
  const CVec3& u = vec0.Normalize();
  mat[0] = 1-u.x()*u.x(); mat[1] = 0-u.x()*u.y(); mat[2] = 0-u.x()*u.z();
  mat[3] = 0-u.y()*u.x(); mat[4] = 1-u.y()*u.y(); mat[5] = 0-u.y()*u.z();
  mat[6] = 0-u.z()*u.x(); mat[7] = 0-u.z()*u.y(); mat[8] = 1-u.z()*u.z();
}

// ----------------------------

dfm2::CMat3d dfm2::Mirror(const CVec3& n)
{
  CVec3 N = n;
  N.SetNormalizedVector();
  return CMat3d::Identity() - 2*dfm2::Mat3_OuterProduct(N,N);
}

dfm2::CMat3d dfm2::RotMatrix_Cartesian(const CVec3& v){
 CMat3d m;
 SetRotMatrix_Cartesian(m,v);
 return m;
}

dfm2::CMat3d dfm2::Mat3(const CVec3& vec0){
  CMat3d m;
  SetSpinTensor(m,vec0);
  return m;
}

dfm2::CMat3d dfm2::Mat3(const CVec3& vec0, const CVec3& vec1){
  CMat3d m;
  SetOuterProduct(m,vec0, vec1);
  return m;
}

dfm2::CMat3d dfm2::Mat3(const CVec3& vec0, const CVec3& vec1, const CVec3& vec2)
{
  CMat3d m;
  double* mat = m.mat;
  mat[0*3+0]=vec0.x(); mat[0*3+1]=vec1.x(); mat[0*3+2]=vec2.x();
  mat[1*3+0]=vec0.y(); mat[1*3+1]=vec1.y(); mat[1*3+2]=vec2.y();
  mat[2*3+0]=vec0.z(); mat[2*3+1]=vec1.z(); mat[2*3+2]=vec2.z();
  return m;
}

dfm2::CMat3d dfm2::Mat3_Spin(const CVec3& vec0){
  CMat3d m;
  SetSpinTensor(m,vec0);
  return m;
}

dfm2::CMat3d dfm2::Mat3_OuterProduct(const CVec3& vec0, const CVec3& vec1 )
{
  CMat3d m;
  SetOuterProduct(m,vec0,vec1);
  return m;
}

dfm2::CMat3d dfm2::Mat3_RotCartesian(const CVec3& vec0)
{
  CMat3d m;
  m.SetRotMatrix_Cartesian(vec0.x(), vec0.y(), vec0.z());
  return m;
}

// ------------------

namespace delfem2 {
  
CVec3 operator* (const CVec3& v, const CMat3d& m){
  return MatVecTrans(m,v);
}
  
CVec3 operator* (const CMat3d& m, const CVec3& v)
{
  return MatVec(m,v);
}
  
}

// ------------------------------

dfm2::CMat3d dfm2::Mat3_MinimumRotation
(const CVec3& V,
 const CVec3& v)
{
  CVec3 ep = V.Normalize();
  CVec3 eq = v.Normalize();
  CVec3 n = ep^eq;
  const double st2 = n*n;
  CMat3d m;
  if( st2 < 1.0e-4f ){
    m.mat[0] = 1.f    +0.5f*(n.x()*n.x()-st2);
    m.mat[1] =    -n.z()+0.5f*(n.x()*n.y());
    m.mat[2] =    +n.y()+0.5f*(n.x()*n.z());
    m.mat[3] =    +n.z()+0.5f*(n.y()*n.x());
    m.mat[4] = 1.f    +0.5f*(n.y()*n.y()-st2);
    m.mat[5] =    -n.x()+0.5f*(n.y()*n.z());
    m.mat[6] =    -n.y()+0.5f*(n.z()*n.x());
    m.mat[7] =    +n.x()+0.5f*(n.z()*n.y());
    m.mat[8] = 1.f    +0.5f*(n.z()*n.z()-st2);
    return m;
  }
  const double st = sqrt(st2);
  const double ct = ep*eq;
  n.SetNormalizedVector();
  m.mat[0] = ct       +(1.f-ct)*n.x()*n.x();
  m.mat[1] =   -n.z()*st+(1.f-ct)*n.x()*n.y();
  m.mat[2] =   +n.y()*st+(1.f-ct)*n.x()*n.z();
  m.mat[3] =   +n.z()*st+(1.f-ct)*n.y()*n.x();
  m.mat[4] = ct       +(1.f-ct)*n.y()*n.y();
  m.mat[5] =   -n.x()*st+(1.f-ct)*n.y()*n.z();
  m.mat[6] =   -n.y()*st+(1.f-ct)*n.z()*n.x();
  m.mat[7] =   +n.x()*st+(1.f-ct)*n.z()*n.y();
  m.mat[8] = ct       +(1.f-ct)*n.z()*n.z();
  return m;
}




void dfm2::Energy_MIPS
(double& E, double dE[3][3], double ddE[3][3][3][3],
 const double c[3][3],
 const double C[3][3])
{
  /*
   double area = TriArea3D(c[0], c[1], c[2]);
   double Area = TriArea3D(C[0], C[1], C[2]);
   double la = SquareDistance3D(c[1], c[2]);
   double lb = SquareDistance3D(c[2], c[0]);
   double lc = SquareDistance3D(c[0], c[1]);
   double lA = SquareDistance3D(C[1], C[2]);
   double lB = SquareDistance3D(C[2], C[0]);
   double lC = SquareDistance3D(C[0], C[1]);
   double cot0 = (-la+lb+lc);
   double cot1 = (+la-lb+lc);
   double cot2 = (+la+lb-lc);
   //  double E_angle = (cot0*lA+cot1*lB+cot2*lC)/(8*Area*area);
   //  double E_area = Area/area + area/Area;
   //  E = E_angle*E_area;
   */
  CVec3 p0(c[0][0],c[0][1],c[0][2]);
  CVec3 p1(c[1][0],c[1][1],c[1][2]);
  CVec3 p2(c[2][0],c[2][1],c[2][2]);
  CVec3 P0(C[0][0],C[0][1],C[0][2]);
  CVec3 P1(C[1][0],C[1][1],C[1][2]);
  CVec3 P2(C[2][0],C[2][1],C[2][2]);
  CVec3 v01 = p1-p0;
  CVec3 v12 = p2-p1;
  CVec3 v20 = p0-p2;
  CVec3 n = v01^v20;
  double area = n.Length()*0.5;
  double Area = ((P1-P0)^(P2-P0)).Length()*0.5;
  double la = (p1-p2)*(p1-p2);
  double lb = (p2-p0)*(p2-p0);
  double lc = (p0-p1)*(p0-p1);
  double lA = (P1-P2)*(P1-P2);
  double lB = (P2-P0)*(P2-P0);
  double lC = (P0-P1)*(P0-P1);
  double cot0 = (-la+lb+lc);
  double cot1 = (+la-lb+lc);
  double cot2 = (+la+lb-lc);
  double tmp0 = 1.0/(8*Area);
  double EC = (cot0*lA+cot1*lB+cot2*lC)*tmp0;
  //  CVector3 dECd0 = (2*p0-p2-p1)*lA*2 + (p2-p1)*lB*2 + (p1-p2)*lC*2;
  //  CVector3 dECd1 = (p2-p0)*lA*2 + (2*p1-p2-p0)*lB*2 + (p0-p2)*lC*2;
  //  CVector3 dECd2 = (p1-p0)*lA*2 + (p0-p1)*lB*2 + (2*p2-p0-p1)*lC*2;
  double t00 = 4*lA*tmp0;
  double t11 = 4*lB*tmp0;
  double t22 = 4*lC*tmp0;
  double t01 = (-2*lA-2*lB+2*lC)*tmp0;
  double t02 = (-2*lA+2*lB-2*lC)*tmp0;
  double t12 = (+2*lA-2*lB-2*lC)*tmp0;
  CVec3 dECd0 = t00*p0 + t01*p1 + t02*p2;
  CVec3 dECd1 = t01*p0 + t11*p1 + t12*p2;
  CVec3 dECd2 = t02*p0 + t12*p1 + t22*p2;
  ////
  dE[0][0]=dECd0.x(); dE[0][1]=dECd0.y(); dE[0][2]=dECd0.z();
  dE[1][0]=dECd1.x(); dE[1][1]=dECd1.y(); dE[1][2]=dECd1.z();
  dE[2][0]=dECd2.x(); dE[2][1]=dECd2.y(); dE[2][2]=dECd2.z();
  
  CMat3d (*op)(const CVec3&, const CVec3&) = Mat3_OuterProduct;
  
  double tmp1 = 0.25/area;
  CVec3 dad0 = ((v20*v12)*v01-(v01*v12)*v20)*tmp1;
  CVec3 dad1 = ((v01*v20)*v12-(v12*v20)*v01)*tmp1;
  CVec3 dad2 = ((v12*v01)*v20-(v20*v01)*v12)*tmp1;
  CMat3d ddad0d0 = (CMat3d::Identity(v12*v12) - op(v12,v12) - 4*op(dad0,dad0))*tmp1;
  CMat3d ddad0d1 = (CMat3d::Identity(v20*v12) - op(v20,v12-v01) - op(v01,v20) - 4*op(dad0,dad1))*tmp1;
  CMat3d ddad0d2 = (CMat3d::Identity(v01*v12) - op(v01,v12-v20) - op(v20,v01) - 4*op(dad0,dad2))*tmp1;
  CMat3d ddad1d0 = (CMat3d::Identity(v12*v20) - op(v12,v20-v01) - op(v01,v12) - 4*op(dad1,dad0))*tmp1;
  CMat3d ddad1d1 = (CMat3d::Identity(v20*v20) - op(v20,v20)                             - 4*op(dad1,dad1))*tmp1;
  CMat3d ddad1d2 = (CMat3d::Identity(v01*v20) - op(v01,v20-v12) - op(v12,v01) - 4*op(dad1,dad2))*tmp1;
  CMat3d ddad2d0 = (CMat3d::Identity(v12*v01) - op(v12,v01-v20) - op(v20,v12) - 4*op(dad2,dad0))*tmp1;
  CMat3d ddad2d1 = (CMat3d::Identity(v20*v01) - op(v20,v01-v12) - op(v12,v20) - 4*op(dad2,dad1))*tmp1;
  CMat3d ddad2d2 = (CMat3d::Identity(v01*v01) - op(v01,v01)                             - 4*op(dad2,dad2))*tmp1;
  
  double ADR = Area/area+area/Area;
  double EA = ADR;
  double dADR = 1.0/Area-Area/(area*area);
  double dEA = dADR;
  double ddADR = 2*Area/(area*area*area);
  double ddEA = ddADR;
  E = EC*EA;
  for(int idim=0;idim<3;++idim){
    dE[0][idim]=EC*dEA*dad0[idim]+EA*dECd0[idim];
    dE[1][idim]=EC*dEA*dad1[idim]+EA*dECd1[idim];
    dE[2][idim]=EC*dEA*dad2[idim]+EA*dECd2[idim];
  }
  CMat3d ddEd0d0 = EC*dEA*ddad0d0 + EC*ddEA*op(dad0,dad0) + EA*CMat3d::Identity(t00) + dEA*op(dad0,dECd0)*2;
  CMat3d ddEd0d1 = EC*dEA*ddad0d1 + EC*ddEA*op(dad0,dad1) + EA*CMat3d::Identity(t01) + dEA*op(dad0,dECd1) + dEA*op(dad1,dECd0);
  CMat3d ddEd0d2 = EC*dEA*ddad0d2 + EC*ddEA*op(dad0,dad2) + EA*CMat3d::Identity(t02) + dEA*op(dad0,dECd2) + dEA*op(dad2,dECd0);
  CMat3d ddEd1d0 = EC*dEA*ddad1d0 + EC*ddEA*op(dad1,dad0) + EA*CMat3d::Identity(t01) + dEA*op(dad1,dECd0) + dEA*op(dad0,dECd1);
  CMat3d ddEd1d1 = EC*dEA*ddad1d1 + EC*ddEA*op(dad1,dad1) + EA*CMat3d::Identity(t11) + dEA*op(dad1,dECd1)*2;
  CMat3d ddEd1d2 = EC*dEA*ddad1d2 + EC*ddEA*op(dad1,dad2) + EA*CMat3d::Identity(t12) + dEA*op(dad1,dECd2) + dEA*op(dad2,dECd1);
  CMat3d ddEd2d0 = EC*dEA*ddad2d0 + EC*ddEA*op(dad2,dad0) + EA*CMat3d::Identity(t02) + dEA*op(dad2,dECd0) + dEA*op(dad0,dECd2);
  CMat3d ddEd2d1 = EC*dEA*ddad2d1 + EC*ddEA*op(dad2,dad1) + EA*CMat3d::Identity(t12) + dEA*op(dad2,dECd1) + dEA*op(dad1,dECd2);
  CMat3d ddEd2d2 = EC*dEA*ddad2d2 + EC*ddEA*op(dad2,dad2) + EA*CMat3d::Identity(t22) + dEA*op(dad2,dECd2)*2;
  
  for(int idim=0;idim<3;++idim){
    for(int jdim=0;jdim<3;++jdim){
      ddE[0][0][idim][jdim] = ddEd0d0(idim,jdim);
      ddE[0][1][idim][jdim] = ddEd0d1(idim,jdim);
      ddE[0][2][idim][jdim] = ddEd0d2(idim,jdim);
      ddE[1][0][idim][jdim] = ddEd1d0(idim,jdim);
      ddE[1][1][idim][jdim] = ddEd1d1(idim,jdim);
      ddE[1][2][idim][jdim] = ddEd1d2(idim,jdim);
      ddE[2][0][idim][jdim] = ddEd2d0(idim,jdim);
      ddE[2][1][idim][jdim] = ddEd2d1(idim,jdim);
      ddE[2][2][idim][jdim] = ddEd2d2(idim,jdim);
    }
  }
}

void CheckEnergyMIPS(){
  double C[3][3];
  for(auto & ino : C){
    ino[0] = (double)rand()/(RAND_MAX+1.0);
    ino[1] = (double)rand()/(RAND_MAX+1.0);
    ino[2] = (double)rand()/(RAND_MAX+1.0);
  }
  double c[3][3];
  dfm2::CMat3d m;
  m.SetRotMatrix_Cartesian(0.3, 1.0, 0.5);
  for(int ino=0;ino<3;++ino){
    m.MatVec(C[ino], c[ino]);
    //        c[ino][0] = C[ino][0];
    //        c[ino][1] = C[ino][1];
    //        c[ino][2] = C[ino][2];
    //        c[ino][0] *= 1.5;
    //        c[ino][1] *= 1.5;
    //        c[ino][2] *= 1.5;
  }
  double E, dE[3][3], ddE[3][3][3][3];
  dfm2::Energy_MIPS(E, dE, ddE,
              c, C);
  std::cout << E << std::endl;
  for(int ino=0;ino<3;++ino){
    for(int idim=0;idim<3;++idim){
      double c1[3][3] = {
        {c[0][0],c[0][1],c[0][2]},
        {c[1][0],c[1][1],c[1][2]},
        {c[2][0],c[2][1],c[2][2]} };
      double eps = 1.0e-4;
      c1[ino][idim] += eps;
      double E1, dE1[3][3], ddE1[3][3][3][3];
      dfm2::Energy_MIPS(E1, dE1, ddE1,
                  c1, C);
      std::cout << " " << (E1-E)/eps << " " << dE[ino][idim] << std::endl;
      for(int jno=0;jno<3;++jno){
        for(int jdim=0;jdim<3;++jdim){
          std::cout << "   " << jno << " " << jdim << "   -->  " << (dE1[jno][jdim]-dE[jno][jdim])/eps << " " << ddE[jno][ino][jdim][idim] << std::endl;
        }
      }
    }
  }
}

dfm2::CMat3d dfm2::Mat3_ParallelTransport
(const CVec3& p0,
 const CVec3& p1,
 const CVec3& q0,
 const CVec3& q1)
{
  return Mat3_MinimumRotation(p1-p0, q1-q0);
}

// moment of inertia around origin triangle vtx (origin,d0,d1,d2) the area_density=1
dfm2::CMat3d dfm2::Mat3_IrotTri
(const CVec3& d0,
 const CVec3& d1,
 const CVec3& d2)
{
  // see http://www.dcs.warwick.ac.uk/~rahil/files/RigidBodySimulation.pdf
  
  CVec3 dv = d0+d1+d2;
  CMat3d I0 = Mat3_OuterProduct(d0,d0) + Mat3_OuterProduct(d1,d1) + Mat3_OuterProduct(d2,d2) + Mat3_OuterProduct(dv,dv);
  double tr0 = I0.Trace();
  CMat3d I = tr0*CMat3d::Identity()-I0;
  
  double darea = ((d1-d0)^(d2-d0)).Length();
  I *= darea/24.0;
  return I;
}

// moment of inertia triangle pyramid with vtx (origin,d0,d1,d2) volume_density = 1
dfm2::CMat3d dfm2::Mat3_IrotTriSolid
(const CVec3& d0,
 const CVec3& d1,
 const CVec3& d2)
{
  // see http://www.dcs.warwick.ac.uk/~rahil/files/RigidBodySimulation.pdf
  
  CVec3 dv = d0+d1+d2;
  CMat3d I0 = Mat3_OuterProduct(d0,d0) + Mat3_OuterProduct(d1,d1) + Mat3_OuterProduct(d2,d2) + Mat3_OuterProduct(dv,dv);
  double tr0 = I0.Trace();
  CMat3d I = tr0*CMat3d::Identity()-I0;
  
  double darea = (d0*(d1^d2));
  I *= darea/120.0;
  return I;
}

dfm2::CMat3d dfm2::Mat3_IrotLineSeg
(const CVec3& d0,
 const CVec3& d1)
{
  CVec3 dv = d1-d0;
  double l = dv.Length();
  CMat3d I;
  {
    I = dv.DLength()*CMat3d::Identity()-Mat3_OuterProduct(dv,dv);
    I *= l/12.0;
  }
  CVec3 p = (d0+d1)*0.5;
  I += l*(p.DLength()*CMat3d::Identity()-Mat3_OuterProduct(p,p));
  return I;
}

dfm2::CMat3d dfm2::Mat3_IrotPoint
(const CVec3& d0)
{
  return (d0.DLength()*CMat3d::Identity()-Mat3_OuterProduct(d0,d0));
}



void dfm2::Mat4_MatTransl(double m[16], const CMat3d& mat, const CVec3& trans)
{
  mat.AffineMatrixTrans(m);
  m[3*4+0] = trans.x();
  m[3*4+1] = trans.y();
  m[3*4+2] = trans.z();
}


void dfm2::Mat4_ScaleMatTransl(double m[16], double scale, const CMat3d& mat, const CVec3& trans)
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
(const CVec3& axis,
 const CVec3& org,
 double rad,
 const CVec3& src,
 const CVec3& dir,
 double pick_tol)
{
  double t = ((org-src)*axis)/(dir*axis);
  CVec3 p0 = src+t*dir;
  double rad0 = (p0-org).Length();
  return fabs(rad - rad0) < pick_tol;
}

bool dfm2::isPickQuad
(const CVec3& p0,const CVec3& p1,const CVec3& p2,const CVec3& p3,
 const dfm2::CVec2& sp, const CVec3& pick_dir,
 const float mMV[16], const float mPj[16],
 double eps)
{
  const dfm2::CVec2 sp0 = dfm2::screenXYProjection(p0, mMV, mPj);
  const dfm2::CVec2 sp1 = dfm2::screenXYProjection(p1, mMV, mPj);
  const dfm2::CVec2 sp2 = dfm2::screenXYProjection(p2, mMV, mPj);
  const dfm2::CVec2 sp3 = dfm2::screenXYProjection(p3, mMV, mPj);
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
  CVec3 n0123 = Normal(p0,p1,p2) + Normal(p1,p2,p3) + Normal(p2,p3,p0) + Normal(p3,p0,p1);
  return n0123 * pick_dir <= 0;
}

int dfm2::PickHandlerRotation_PosQuat
(const CVec3& src, const CVec3& dir,
 const CVec3& pos, const double quat[4], double rad,
 double tol)
{
  CVec3 ax = QuatVec(quat,CVec3(1,0,0));
  CVec3 ay = QuatVec(quat,CVec3(0,1,0));
  CVec3 az = QuatVec(quat,CVec3(0,0,1));
  CVec3 px,qx; Nearest_Line_Circle(px,qx, src,dir, pos,ax, rad);
  CVec3 py,qy; Nearest_Line_Circle(py,qy, src,dir, pos,ay, rad);
  CVec3 pz,qz; Nearest_Line_Circle(pz,qz, src,dir, pos,az, rad);
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
(const CVec3& src, const CVec3& dir,
 const double mat[16], double rad,
 double tol)
{
  CVec3 ax = Mat4Vec(mat,CVec3(1,0,0));
  CVec3 ay = Mat4Vec(mat,CVec3(0,1,0));
  CVec3 az = Mat4Vec(mat,CVec3(0,0,1));
  CVec3 pos(mat[3],mat[7],mat[11]);
  CVec3 px,qx; Nearest_Line_Circle(px,qx, src,dir, pos,ax, rad);
  CVec3 py,qy; Nearest_Line_Circle(py,qy, src,dir, pos,ay, rad);
  CVec3 pz,qz; Nearest_Line_Circle(pz,qz, src,dir, pos,az, rad);
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
 const dfm2::CVec2& sp0,
 const dfm2::CVec2& sp1,
 const CVec3& pos,
 const float mMV[16], const float mPj[16])
{
  if( ielem>=0 && ielem<3 ){
    double vi[3] = {0,0,0}; vi[ielem] = 1;
    double vo[3]; dfm2::QuatVec(vo, quat, vi);
    CVec3 v0(0,0,0); v0[ielem] = 1;
    CVec3 v1(vo[0],vo[1],vo[2]); v1.SetNormalizedVector();
    double ar = -DragCircle(sp0,sp1, pos, v1, mMV, mPj);
    double dq[4] = { cos(ar*0.5), v0.x()*sin(ar*0.5), v0.y()*sin(ar*0.5), v0.z()*sin(ar*0.5) };
    double qtmp[4]; QuatQuat(qtmp, dq, quat);
    QuatCopy(quat,qtmp);
    return true;
  }
  return false;
}

bool dfm2::DragHandlerRot_Mat4
(double quat[4], int ielem,
 const dfm2::CVec2& sp0, const dfm2::CVec2& sp1, double mat[16],
 const float mMV[16], const float mPj[16])
{
  if( ielem>=0 && ielem<3 ){
    double vi[3] = {0,0,0}; vi[ielem] = 1;
    double vo[3]; Mat4Vec3(vo, mat, vi);
    CVec3 v0(0,0,0); v0[ielem] = 1;
    CVec3 v1(vo[0],vo[1],vo[2]); v1.SetNormalizedVector();
    CVec3 pos(mat[3],mat[7],mat[11]);
    const double ar = DragCircle(sp0,sp1, pos, v1, mMV, mPj);
    const double dq[4] = { cos(ar*0.5), v0.x()*sin(ar*0.5), v0.y()*sin(ar*0.5), v0.z()*sin(ar*0.5) };
    double qtmp[4]; QuatQuat(qtmp, quat, dq);
    QuatCopy(quat,qtmp);
    return true;
  }
  return false;
}

bool dfm2::isPick_AxisHandler
(const dfm2::CVec2& sp,
 const CVec3& p,
 const CVec3& axis,
 double len,
 const float* mMV,
 const float* mPj,
 double pick_tol)
{
  dfm2::CVec2 sp0 = dfm2::screenXYProjection(p+len*axis, mMV, mPj);
  dfm2::CVec2 sp1 = dfm2::screenXYProjection(p-len*axis, mMV, mPj);
  double sdist = GetDist_LineSeg_Point(sp, sp0, sp1);
  return sdist < pick_tol;
}

dfm2::CVec3 dfm2::drag_AxisHandler
(const dfm2::CVec2& sp0,
 const dfm2::CVec2& sp1,
 const CVec3& p,
 const CVec3& axis,
 double len,
 const float* mMV,
 const float* mPj)
{
  dfm2::CVec2 spa0 = dfm2::screenXYProjection(p+len*axis, mMV, mPj);
  dfm2::CVec2 spa1 = dfm2::screenXYProjection(p-len*axis, mMV, mPj);
  double r = (spa0-spa1)*(sp1-sp0)/(spa0-spa1).SqLength();
  return r*axis*len;
}

double dfm2::DragCircle
(const dfm2::CVec2& sp0,
 const dfm2::CVec2& sp1,
 const CVec3& p,
 const CVec3& axis,
 const float* mMV,
 const float* mPj)
{
  dfm2::CVec2 spo0 = dfm2::screenXYProjection(p, mMV, mPj);
  double area = Area_Tri(sp0, spo0, sp1);
  double angl = area / ( (sp0-spo0).Length() * (sp1-spo0).Length() );
  {
    CVec3 a3 = screenUnProjectionDirection(axis,mMV,mPj);
    if( a3.z() < 0 ){ angl *= -1; }
  }
  return angl;
  //  CMatrix3 R; R.SetRotMatrix_Cartesian(angl*axis);
  //  return R;
}

bool dfm2::isPickPoint
(const dfm2::CVec2& sp,
 const CVec3& p,
 const float* mMV,
 const float* mPj,
 double pick_tol)
{
  dfm2::CVec2 sp0 = dfm2::screenXYProjection(p, mMV, mPj);
  return (sp - sp0).Length() < pick_tol;
}

bool dfm2::isPickCircle
(const dfm2::CVec2& sp,
 const CVec3& p,
 const CVec3& axis,
 double r,
 const float* mMV,
 const float* mPj,
 double pick_tol)
{
  const int ndiv = 32;
  double rdiv = 3.1415*2.0/ndiv;
  CVec3 x,y; GetVertical2Vector(axis, x, y);
  for(int idiv=0;idiv<ndiv+1;idiv++){
    int jdiv = idiv+1;
    CVec3 p0 = p+(r*sin(rdiv*idiv))*x+(r*cos(rdiv*idiv))*y;
    CVec3 p1 = p+(r*sin(rdiv*jdiv))*x+(r*cos(rdiv*jdiv))*y;
    dfm2::CVec2 sp0 = dfm2::screenXYProjection(p0, mMV, mPj);
    dfm2::CVec2 sp1 = dfm2::screenXYProjection(p1, mMV, mPj);
    double sdist = GetDist_LineSeg_Point(sp, sp0, sp1);
    if( sdist < pick_tol ){ return true; }
  }
  return false;
}

