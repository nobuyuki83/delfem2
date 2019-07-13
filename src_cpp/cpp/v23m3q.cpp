/*
 * Copyright (c) 2019 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */


#include "delfem2/v23m3q.h"

CVector3 GetCartesianRotationVector(const CMatrix3& m)
{
  const double* mat = m.mat;
  CVector3 a;
  a.x = mat[7]-mat[5];
  a.y = mat[2]-mat[6];
  a.z = mat[3]-mat[1];
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

CVector3 GetSpinVector(const CMatrix3& m)
{
  const double* mat = m.mat;
  CVector3 r;
  r.x = (mat[7]-mat[5])*0.5;
  r.y = (mat[2]-mat[6])*0.5;
  r.z = (mat[3]-mat[1])*0.5;
  return r;
}

CVector3 MatVec(const CMatrix3& m, const CVector3& vec0)
{
  const double* mat = m.mat;
  CVector3 vec1;
  vec1.x = mat[0]*vec0.x + mat[1]*vec0.y + mat[2]*vec0.z;
  vec1.y = mat[3]*vec0.x + mat[4]*vec0.y + mat[5]*vec0.z;
  vec1.z = mat[6]*vec0.x + mat[7]*vec0.y + mat[8]*vec0.z;
  return vec1;
}

CVector3 MatVecTrans(const CMatrix3& m, const CVector3& vec0)
{
  CVector3 vec1;
  const double* mat = m.mat;
  vec1.x = mat[0]*vec0.x + mat[3]*vec0.y + mat[6]*vec0.z;
  vec1.y = mat[1]*vec0.x + mat[4]*vec0.y + mat[7]*vec0.z;
  vec1.z = mat[2]*vec0.x + mat[5]*vec0.y + mat[8]*vec0.z;
  return vec1;
}

////////////

void SetDiag(CMatrix3& m, const CVector3& d)
{
  double* mat = m.mat;
  mat[0*3+0] = d.x;
  mat[1*3+1] = d.y;
  mat[2*3+2] = d.z;
}

void SetRotMatrix_Cartesian(CMatrix3& m, const CVector3& v)
{
  const double vec[3] = { v.x, v.y, v.z };
  m.SetRotMatrix_Cartesian(vec);
}

void SetSpinTensor(CMatrix3& m, const CVector3& vec0)
{
  double* mat = m.mat;
  mat[0] =  0;       mat[1] = -vec0.z;   mat[2] = +vec0.y;
  mat[3] = +vec0.z;  mat[4] = 0;         mat[5] = -vec0.x;
  mat[6] = -vec0.y;  mat[7] = +vec0.x;   mat[8] = 0;
}

void SetOuterProduct(CMatrix3& m, const CVector3& vec0, const CVector3& vec1 )
{
  double* mat = m.mat;
  mat[0] = vec0.x*vec1.x; mat[1] = vec0.x*vec1.y; mat[2] = vec0.x*vec1.z;
  mat[3] = vec0.y*vec1.x; mat[4] = vec0.y*vec1.y; mat[5] = vec0.y*vec1.z;
  mat[6] = vec0.z*vec1.x; mat[7] = vec0.z*vec1.y; mat[8] = vec0.z*vec1.z;
}

void SetProjection(CMatrix3& m, const CVector3& vec0)
{
  double* mat = m.mat;
  const CVector3& u = vec0.Normalize();
  mat[0] = 1-u.x*u.x; mat[1] = 0-u.x*u.y; mat[2] = 0-u.x*u.z;
  mat[3] = 0-u.y*u.x; mat[4] = 1-u.y*u.y; mat[5] = 0-u.y*u.z;
  mat[6] = 0-u.z*u.x; mat[7] = 0-u.z*u.y; mat[8] = 1-u.z*u.z;
}

//////////

CMatrix3 Mirror(const CVector3& n){
  CVector3 N = n;
  N.SetNormalizedVector();
  return CMatrix3::Identity() - 2*Mat3_OuterProduct(N,N);
}

CMatrix3 RotMatrix_Cartesian(const CVector3& v){
 CMatrix3 m;
 SetRotMatrix_Cartesian(m,v);
 return m;
}

CMatrix3 Mat3(const CVector3& vec0){
  CMatrix3 m;
  SetSpinTensor(m,vec0);
  return m;
}

CMatrix3 Mat3(const CVector3& vec0, const CVector3& vec1){
  CMatrix3 m;
  SetOuterProduct(m,vec0, vec1);
  return m;
}

CMatrix3 Mat3(const CVector3& vec0, const CVector3& vec1, const CVector3& vec2)
{
  CMatrix3 m;
  double* mat = m.mat;
  mat[0*3+0]=vec0.x; mat[0*3+1]=vec1.x; mat[0*3+2]=vec2.x;
  mat[1*3+0]=vec0.y; mat[1*3+1]=vec1.y; mat[1*3+2]=vec2.y;
  mat[2*3+0]=vec0.z; mat[2*3+1]=vec1.z; mat[2*3+2]=vec2.z;
  return m;
}

CMatrix3 Mat3_Spin(const CVector3& vec0){
  CMatrix3 m;
  SetSpinTensor(m,vec0);
  return m;
}

CMatrix3 Mat3_OuterProduct(const CVector3& vec0, const CVector3& vec1 )
{
  CMatrix3 m;
  SetOuterProduct(m,vec0,vec1);
  return m;
}

CMatrix3 Mat3_RotCartesian(const CVector3& vec0)
{
  CMatrix3 m;
  m.SetRotMatrix_Cartesian(vec0.x, vec0.y, vec0.z);
  return m;
}

////////////

CVector3 operator* (const CVector3& v, const CMatrix3& m){
  return MatVecTrans(m,v);
}
CVector3 operator* (const CMatrix3& m, const CVector3& v)
{
  return MatVec(m,v);
}

//////////////////////////////////////////////////////////////////////

CMatrix3 Mat3_MinimumRotation
(const CVector3& V,
 const CVector3& v)
{
  CVector3 ep = V.Normalize();
  CVector3 eq = v.Normalize();
  CVector3 n = ep^eq;
  const double st2 = n*n;
  CMatrix3 m;
  if( st2 < 1.0e-4f ){
    m.mat[0] = 1.f    +0.5f*(n.x*n.x-st2);
    m.mat[1] =    -n.z+0.5f*(n.x*n.y);
    m.mat[2] =    +n.y+0.5f*(n.x*n.z);
    m.mat[3] =    +n.z+0.5f*(n.y*n.x);
    m.mat[4] = 1.f    +0.5f*(n.y*n.y-st2);
    m.mat[5] =    -n.x+0.5f*(n.y*n.z);
    m.mat[6] =    -n.y+0.5f*(n.z*n.x);
    m.mat[7] =    +n.x+0.5f*(n.z*n.y);
    m.mat[8] = 1.f    +0.5f*(n.z*n.z-st2);
    return m;
  }
  const double st = sqrt(st2);
  const double ct = ep*eq;
  n.SetNormalizedVector();
  m.mat[0] = ct       +(1.f-ct)*n.x*n.x;
  m.mat[1] =   -n.z*st+(1.f-ct)*n.x*n.y;
  m.mat[2] =   +n.y*st+(1.f-ct)*n.x*n.z;
  m.mat[3] =   +n.z*st+(1.f-ct)*n.y*n.x;
  m.mat[4] = ct       +(1.f-ct)*n.y*n.y;
  m.mat[5] =   -n.x*st+(1.f-ct)*n.y*n.z;
  m.mat[6] =   -n.y*st+(1.f-ct)*n.z*n.x;
  m.mat[7] =   +n.x*st+(1.f-ct)*n.z*n.y;
  m.mat[8] = ct       +(1.f-ct)*n.z*n.z;
  return m;
}




void Energy_MIPS
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
  CVector3 p0(c[0][0],c[0][1],c[0][2]);
  CVector3 p1(c[1][0],c[1][1],c[1][2]);
  CVector3 p2(c[2][0],c[2][1],c[2][2]);
  CVector3 P0(C[0][0],C[0][1],C[0][2]);
  CVector3 P1(C[1][0],C[1][1],C[1][2]);
  CVector3 P2(C[2][0],C[2][1],C[2][2]);
  CVector3 v01 = p1-p0;
  CVector3 v12 = p2-p1;
  CVector3 v20 = p0-p2;
  CVector3 n = v01^v20;
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
  CVector3 dECd0 = t00*p0 + t01*p1 + t02*p2;
  CVector3 dECd1 = t01*p0 + t11*p1 + t12*p2;
  CVector3 dECd2 = t02*p0 + t12*p1 + t22*p2;
  ////
  dE[0][0]=dECd0.x; dE[0][1]=dECd0.y; dE[0][2]=dECd0.z;
  dE[1][0]=dECd1.x; dE[1][1]=dECd1.y; dE[1][2]=dECd1.z;
  dE[2][0]=dECd2.x; dE[2][1]=dECd2.y; dE[2][2]=dECd2.z;
  
  CMatrix3 (*op)(const CVector3&, const CVector3&) = Mat3_OuterProduct;
  
  double tmp1 = 0.25/area;
  CVector3 dad0 = ((v20*v12)*v01-(v01*v12)*v20)*tmp1;
  CVector3 dad1 = ((v01*v20)*v12-(v12*v20)*v01)*tmp1;
  CVector3 dad2 = ((v12*v01)*v20-(v20*v01)*v12)*tmp1;
  CMatrix3 ddad0d0 = (CMatrix3::Identity(v12*v12) - op(v12,v12) - 4*op(dad0,dad0))*tmp1;
  CMatrix3 ddad0d1 = (CMatrix3::Identity(v20*v12) - op(v20,v12-v01) - op(v01,v20) - 4*op(dad0,dad1))*tmp1;
  CMatrix3 ddad0d2 = (CMatrix3::Identity(v01*v12) - op(v01,v12-v20) - op(v20,v01) - 4*op(dad0,dad2))*tmp1;
  CMatrix3 ddad1d0 = (CMatrix3::Identity(v12*v20) - op(v12,v20-v01) - op(v01,v12) - 4*op(dad1,dad0))*tmp1;
  CMatrix3 ddad1d1 = (CMatrix3::Identity(v20*v20) - op(v20,v20)                             - 4*op(dad1,dad1))*tmp1;
  CMatrix3 ddad1d2 = (CMatrix3::Identity(v01*v20) - op(v01,v20-v12) - op(v12,v01) - 4*op(dad1,dad2))*tmp1;
  CMatrix3 ddad2d0 = (CMatrix3::Identity(v12*v01) - op(v12,v01-v20) - op(v20,v12) - 4*op(dad2,dad0))*tmp1;
  CMatrix3 ddad2d1 = (CMatrix3::Identity(v20*v01) - op(v20,v01-v12) - op(v12,v20) - 4*op(dad2,dad1))*tmp1;
  CMatrix3 ddad2d2 = (CMatrix3::Identity(v01*v01) - op(v01,v01)                             - 4*op(dad2,dad2))*tmp1;
  
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
  CMatrix3 ddEd0d0 = EC*dEA*ddad0d0 + EC*ddEA*op(dad0,dad0) + EA*CMatrix3::Identity(t00) + dEA*op(dad0,dECd0)*2;
  CMatrix3 ddEd0d1 = EC*dEA*ddad0d1 + EC*ddEA*op(dad0,dad1) + EA*CMatrix3::Identity(t01) + dEA*op(dad0,dECd1) + dEA*op(dad1,dECd0);
  CMatrix3 ddEd0d2 = EC*dEA*ddad0d2 + EC*ddEA*op(dad0,dad2) + EA*CMatrix3::Identity(t02) + dEA*op(dad0,dECd2) + dEA*op(dad2,dECd0);
  CMatrix3 ddEd1d0 = EC*dEA*ddad1d0 + EC*ddEA*op(dad1,dad0) + EA*CMatrix3::Identity(t01) + dEA*op(dad1,dECd0) + dEA*op(dad0,dECd1);
  CMatrix3 ddEd1d1 = EC*dEA*ddad1d1 + EC*ddEA*op(dad1,dad1) + EA*CMatrix3::Identity(t11) + dEA*op(dad1,dECd1)*2;
  CMatrix3 ddEd1d2 = EC*dEA*ddad1d2 + EC*ddEA*op(dad1,dad2) + EA*CMatrix3::Identity(t12) + dEA*op(dad1,dECd2) + dEA*op(dad2,dECd1);
  CMatrix3 ddEd2d0 = EC*dEA*ddad2d0 + EC*ddEA*op(dad2,dad0) + EA*CMatrix3::Identity(t02) + dEA*op(dad2,dECd0) + dEA*op(dad0,dECd2);
  CMatrix3 ddEd2d1 = EC*dEA*ddad2d1 + EC*ddEA*op(dad2,dad1) + EA*CMatrix3::Identity(t12) + dEA*op(dad2,dECd1) + dEA*op(dad1,dECd2);
  CMatrix3 ddEd2d2 = EC*dEA*ddad2d2 + EC*ddEA*op(dad2,dad2) + EA*CMatrix3::Identity(t22) + dEA*op(dad2,dECd2)*2;
  
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
  for(int ino=0;ino<3;++ino){
    C[ino][0] = (double)rand()/(RAND_MAX+1.0);
    C[ino][1] = (double)rand()/(RAND_MAX+1.0);
    C[ino][2] = (double)rand()/(RAND_MAX+1.0);
  }
  double c[3][3];
  CMatrix3 m;
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
  Energy_MIPS(E, dE, ddE, c, C);
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
      Energy_MIPS(E1, dE1, ddE1, c1, C);
      std::cout << " " << (E1-E)/eps << " " << dE[ino][idim] << std::endl;
      for(int jno=0;jno<3;++jno){
        for(int jdim=0;jdim<3;++jdim){
          std::cout << "   " << jno << " " << jdim << "   -->  " << (dE1[jno][jdim]-dE[jno][jdim])/eps << " " << ddE[jno][ino][jdim][idim] << std::endl;
        }
      }
    }
  }
}

CMatrix3 Mat3_ParallelTransport
(const CVector3& p0,
 const CVector3& p1,
 const CVector3& q0,
 const CVector3& q1)
{
  return Mat3_MinimumRotation(p1-p0, q1-q0);
}

// moment of inertia around origin triangle vtx (origin,d0,d1,d2) the area_density=1
CMatrix3 Mat3_IrotTri
(const CVector3& d0,
 const CVector3& d1,
 const CVector3& d2)
{
  // see http://www.dcs.warwick.ac.uk/~rahil/files/RigidBodySimulation.pdf
  
  CVector3 dv = d0+d1+d2;
  CMatrix3 I0 = Mat3_OuterProduct(d0,d0) + Mat3_OuterProduct(d1,d1) + Mat3_OuterProduct(d2,d2) + Mat3_OuterProduct(dv,dv);
  double tr0 = I0.Trace();
  CMatrix3 I = tr0*CMatrix3::Identity()-I0;
  
  double darea = ((d1-d0)^(d2-d0)).Length();
  I *= darea/24.0;
  return I;
}

// moment of inertia triangle pyramid with vtx (origin,d0,d1,d2) volume_density = 1
CMatrix3 Mat3_IrotTriSolid
(const CVector3& d0,
 const CVector3& d1,
 const CVector3& d2)
{
  // see http://www.dcs.warwick.ac.uk/~rahil/files/RigidBodySimulation.pdf
  
  CVector3 dv = d0+d1+d2;
  CMatrix3 I0 = Mat3_OuterProduct(d0,d0) + Mat3_OuterProduct(d1,d1) + Mat3_OuterProduct(d2,d2) + Mat3_OuterProduct(dv,dv);
  double tr0 = I0.Trace();
  CMatrix3 I = tr0*CMatrix3::Identity()-I0;
  
  double darea = (d0*(d1^d2));
  I *= darea/120.0;
  return I;
}

CMatrix3 Mat3_IrotLineSeg
(const CVector3& d0,
 const CVector3& d1)
{
  CVector3 dv = d1-d0;
  double l = dv.Length();
  CMatrix3 I;
  {
    I = dv.DLength()*CMatrix3::Identity()-Mat3_OuterProduct(dv,dv);
    I *= l/12.0;
  }
  CVector3 p = (d0+d1)*0.5;
  I += l*(p.DLength()*CMatrix3::Identity()-Mat3_OuterProduct(p,p));
  return I;
}

CMatrix3 Mat3_IrotPoint
(const CVector3& d0)
{
  return (d0.DLength()*CMatrix3::Identity()-Mat3_OuterProduct(d0,d0));
}



void AffineMatrixTrans(double m[16], const CMatrix3& mat, const CVector3& trans)
{
  mat.AffineMatrixTrans(m);
  m[3*4+0] = trans.x;
  m[3*4+1] = trans.y;
  m[3*4+2] = trans.z;
}
