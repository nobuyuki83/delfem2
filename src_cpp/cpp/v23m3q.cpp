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
  return CMatrix3::Identity() - 2*OuterProduct(N,N);
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

CMatrix3 Spin(const CVector3& vec0){
  CMatrix3 m;
  SetSpinTensor(m,vec0);
  return m;
}

CMatrix3 OuterProduct(const CVector3& vec0, const CVector3& vec1 )
{
  CMatrix3 m;
  SetOuterProduct(m,vec0,vec1);
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

CMatrix3 MinimumRotation
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
  
  
  double tmp1 = 0.25/area;
  CVector3 dad0 = ((v20*v12)*v01-(v01*v12)*v20)*tmp1;
  CVector3 dad1 = ((v01*v20)*v12-(v12*v20)*v01)*tmp1;
  CVector3 dad2 = ((v12*v01)*v20-(v20*v01)*v12)*tmp1;
  CMatrix3 ddad0d0 = (CMatrix3::Identity(v12*v12) - OuterProduct(v12,v12)                             - 4*OuterProduct(dad0,dad0))*tmp1;
  CMatrix3 ddad0d1 = (CMatrix3::Identity(v20*v12) - OuterProduct(v20,v12-v01) - OuterProduct(v01,v20) - 4*OuterProduct(dad0,dad1))*tmp1;
  CMatrix3 ddad0d2 = (CMatrix3::Identity(v01*v12) - OuterProduct(v01,v12-v20) - OuterProduct(v20,v01) - 4*OuterProduct(dad0,dad2))*tmp1;
  CMatrix3 ddad1d0 = (CMatrix3::Identity(v12*v20) - OuterProduct(v12,v20-v01) - OuterProduct(v01,v12) - 4*OuterProduct(dad1,dad0))*tmp1;
  CMatrix3 ddad1d1 = (CMatrix3::Identity(v20*v20) - OuterProduct(v20,v20)                             - 4*OuterProduct(dad1,dad1))*tmp1;
  CMatrix3 ddad1d2 = (CMatrix3::Identity(v01*v20) - OuterProduct(v01,v20-v12) - OuterProduct(v12,v01) - 4*OuterProduct(dad1,dad2))*tmp1;
  CMatrix3 ddad2d0 = (CMatrix3::Identity(v12*v01) - OuterProduct(v12,v01-v20) - OuterProduct(v20,v12) - 4*OuterProduct(dad2,dad0))*tmp1;
  CMatrix3 ddad2d1 = (CMatrix3::Identity(v20*v01) - OuterProduct(v20,v01-v12) - OuterProduct(v12,v20) - 4*OuterProduct(dad2,dad1))*tmp1;
  CMatrix3 ddad2d2 = (CMatrix3::Identity(v01*v01) - OuterProduct(v01,v01)                             - 4*OuterProduct(dad2,dad2))*tmp1;
  
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
  CMatrix3 ddEd0d0 = EC*dEA*ddad0d0 + EC*ddEA*OuterProduct(dad0,dad0) + EA*CMatrix3::Identity(t00) + dEA*OuterProduct(dad0,dECd0)*2;
  CMatrix3 ddEd0d1 = EC*dEA*ddad0d1 + EC*ddEA*OuterProduct(dad0,dad1) + EA*CMatrix3::Identity(t01) + dEA*OuterProduct(dad0,dECd1) + dEA*OuterProduct(dad1,dECd0);
  CMatrix3 ddEd0d2 = EC*dEA*ddad0d2 + EC*ddEA*OuterProduct(dad0,dad2) + EA*CMatrix3::Identity(t02) + dEA*OuterProduct(dad0,dECd2) + dEA*OuterProduct(dad2,dECd0);
  CMatrix3 ddEd1d0 = EC*dEA*ddad1d0 + EC*ddEA*OuterProduct(dad1,dad0) + EA*CMatrix3::Identity(t01) + dEA*OuterProduct(dad1,dECd0) + dEA*OuterProduct(dad0,dECd1);
  CMatrix3 ddEd1d1 = EC*dEA*ddad1d1 + EC*ddEA*OuterProduct(dad1,dad1) + EA*CMatrix3::Identity(t11) + dEA*OuterProduct(dad1,dECd1)*2;
  CMatrix3 ddEd1d2 = EC*dEA*ddad1d2 + EC*ddEA*OuterProduct(dad1,dad2) + EA*CMatrix3::Identity(t12) + dEA*OuterProduct(dad1,dECd2) + dEA*OuterProduct(dad2,dECd1);
  CMatrix3 ddEd2d0 = EC*dEA*ddad2d0 + EC*ddEA*OuterProduct(dad2,dad0) + EA*CMatrix3::Identity(t02) + dEA*OuterProduct(dad2,dECd0) + dEA*OuterProduct(dad0,dECd2);
  CMatrix3 ddEd2d1 = EC*dEA*ddad2d1 + EC*ddEA*OuterProduct(dad2,dad1) + EA*CMatrix3::Identity(t12) + dEA*OuterProduct(dad2,dECd1) + dEA*OuterProduct(dad1,dECd2);
  CMatrix3 ddEd2d2 = EC*dEA*ddad2d2 + EC*ddEA*OuterProduct(dad2,dad2) + EA*CMatrix3::Identity(t22) + dEA*OuterProduct(dad2,dECd2)*2;
  
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

CMatrix3 ParallelTransport
(const CVector3& p0,
 const CVector3& p1,
 const CVector3& q0,
 const CVector3& q1)
{
  return MinimumRotation(p1-p0, q1-q0);
}

// moment of inertia around origin triangle vtx (origin,d0,d1,d2) the area_density=1
CMatrix3 Irot_Tri
(const CVector3& d0,
 const CVector3& d1,
 const CVector3& d2)
{
  // see http://www.dcs.warwick.ac.uk/~rahil/files/RigidBodySimulation.pdf
  
  CVector3 dv = d0+d1+d2;
  CMatrix3 I0 = OuterProduct(d0,d0) + OuterProduct(d1,d1) + OuterProduct(d2,d2) + OuterProduct(dv,dv);
  double tr0 = I0.Trace();
  CMatrix3 I = tr0*CMatrix3::Identity()-I0;
  
  double darea = ((d1-d0)^(d2-d0)).Length();
  I *= darea/24.0;
  return I;
}

// moment of inertia triangle pyramid with vtx (origin,d0,d1,d2) volume_density = 1
CMatrix3 Irot_TriSolid
(const CVector3& d0,
 const CVector3& d1,
 const CVector3& d2)
{
  // see http://www.dcs.warwick.ac.uk/~rahil/files/RigidBodySimulation.pdf
  
  CVector3 dv = d0+d1+d2;
  CMatrix3 I0 = OuterProduct(d0,d0) + OuterProduct(d1,d1) + OuterProduct(d2,d2) + OuterProduct(dv,dv);
  double tr0 = I0.Trace();
  CMatrix3 I = tr0*CMatrix3::Identity()-I0;
  
  double darea = (d0*(d1^d2));
  I *= darea/120.0;
  return I;
}

CMatrix3 Irot_LineSeg
(const CVector3& d0,
 const CVector3& d1)
{
  CVector3 dv = d1-d0;
  double l = dv.Length();
  CMatrix3 I;
  {
    I = dv.DLength()*CMatrix3::Identity()-OuterProduct(dv,dv);
    I *= l/12.0;
  }
  CVector3 p = (d0+d1)*0.5;
  I += l*(p.DLength()*CMatrix3::Identity()-OuterProduct(p,p));
  return I;
}

CMatrix3 Irot_Point
(const CVector3& d0)
{
  return (d0.DLength()*CMatrix3::Identity()-OuterProduct(d0,d0));
}


///////////////////////////////////////////////////////////////////


void ConstraintProjection_Rigid3D
(double* aXYZt,
 double stiffness,
 const int* clstr_ind, int nclstr_ind,
 const int* clstr,     int nclstr0,
 const double* aXYZ0,   int nXYZ0)
{
  const int nclstr = nclstr_ind-1;
  for(int iclstr=0;iclstr<nclstr;++iclstr){
    CVector3 pc(0, 0, 0), qc(0, 0, 0);
    for (int iip=clstr_ind[iclstr];iip<clstr_ind[iclstr+1]; iip++){
      const int ip = clstr[iip];
      qc += CVector3(aXYZ0[ip*3+0],aXYZ0[ip*3+1],aXYZ0[ip*3+2]);
      pc += CVector3(aXYZt[ip*3+0],aXYZt[ip*3+1],aXYZt[ip*3+2]);
    }
    qc /= (clstr_ind[iclstr+1]-clstr_ind[iclstr]);
    pc /= (clstr_ind[iclstr+1]-clstr_ind[iclstr]);
    
    double A[9] = { 0,0,0,  0,0,0, 0,0,0 };
    for (int iip=clstr_ind[iclstr];iip<clstr_ind[iclstr+1]; iip++){
      const int ip = clstr[iip];
      const CVector3 dq = CVector3(aXYZ0[ip*3+0],aXYZ0[ip*3+1],aXYZ0[ip*3+2])-qc; // undeform
      const CVector3 dp = CVector3(aXYZt[ip*3+0],aXYZt[ip*3+1],aXYZ0[ip*3+2])-pc; // deform
      A[0*3+0] += dp[0]*dq[0];
      A[0*3+1] += dp[0]*dq[1];
      A[0*3+2] += dp[0]*dq[2];
      A[1*3+0] += dp[1]*dq[0];
      A[1*3+1] += dp[1]*dq[1];
      A[1*3+2] += dp[1]*dq[2];
      A[2*3+0] += dp[2]*dq[0];
      A[2*3+1] += dp[2]*dq[1];
      A[2*3+2] += dp[2]*dq[2];
    }
    double R[9]; GetRotPolarDecomp(R,A, 20);
    //    std::cout << R[0] << " " << R[1] << " " << R[2] << " " << R[3] << std::endl;
    
    for (int iip=clstr_ind[iclstr];iip<clstr_ind[iclstr+1]; iip++){
      const int ip = clstr[iip];
      CVector3 dq = CVector3(aXYZ0[ip*3+0],aXYZ0[ip*3+1],aXYZ0[ip*3+2])-qc;
      CVector3 pg = pc+MatVec(R, dq); // goal position
      CVector3 pg2 = stiffness*pg+(1-stiffness)*CVector3(aXYZt[ip*3+0],aXYZt[ip*3+1],aXYZt[ip*3+2]);
      aXYZt[ip*3+0] = pg2.x;
      aXYZt[ip*3+1] = pg2.y;
      aXYZt[ip*3+2] = pg2.z;
    }
  }
}

void ConstraintProjection_Rigid2D
(double* aXYt,
 double stiffness,
 const int* clstr_ind, int nclstr_ind,
 const int* clstr,     int nclstr0,
 const double* aXY0,   int nXY0)
{
  const int nclstr = nclstr_ind-1;
  for(int iclstr=0;iclstr<nclstr;++iclstr){
    CVector2 pc(0, 0), qc(0, 0);
    for (int iip=clstr_ind[iclstr];iip<clstr_ind[iclstr+1]; iip++){
      const int ip = clstr[iip];
      qc += CVector2(aXY0[ip*2+0],aXY0[ip*2+1]);
      pc += CVector2(aXYt[ip*2+0],aXYt[ip*2+1]);
    }
    qc /= (clstr_ind[iclstr+1]-clstr_ind[iclstr]);
    pc /= (clstr_ind[iclstr+1]-clstr_ind[iclstr]);
    
    double A[4] = { 0, 0, 0, 0 };
    for (int iip=clstr_ind[iclstr];iip<clstr_ind[iclstr+1]; iip++){
      const int ip = clstr[iip];
      const CVector2 dq = CVector2(aXY0[ip*2+0],aXY0[ip*2+1])-qc; // undeform
      const CVector2 dp = CVector2(aXYt[ip*2+0],aXYt[ip*2+1])-pc; // deform
      A[0*2+0] += dp[0]*dq[0];
      A[0*2+1] += dp[0]*dq[1];
      A[1*2+0] += dp[1]*dq[0];
      A[1*2+1] += dp[1]*dq[1];
    }
    double R[4]; RotationalComponentOfMatrix2(R,A);
    //    std::cout << R[0] << " " << R[1] << " " << R[2] << " " << R[3] << std::endl;
    
    for (int iip=clstr_ind[iclstr];iip<clstr_ind[iclstr+1]; iip++){
      const int ip = clstr[iip];
      CVector2 dq = CVector2(aXY0[ip*2+0],aXY0[ip*2+1])-qc;
      CVector2 pg = pc+matVec(R, dq); // goal position
      CVector2 pg2 = stiffness*pg+(1-stiffness)*CVector2(aXYt[ip*2+0],aXYt[ip*2+1]);
      aXYt[ip*2+0] = pg2.x;
      aXYt[ip*2+1] = pg2.y;
    }
  }
}

void PBD_ConstraintProjection_Strain
(double C[3],
 double dCdp[3][9],
 const double P[3][2], // (in) undeformed triangle vertex positions
 const double p[3][3] // (in) deformed triangle vertex positions
)
{
  const CVector3 Gd0( P[1][0]-P[0][0], P[1][1]-P[0][1], 0.0 );
  const CVector3 Gd1( P[2][0]-P[0][0], P[2][1]-P[0][1], 0.0 );
  CVector3 Gd2 = Cross(Gd0, Gd1);
  const double Area = Gd2.Length()*0.5;
  Gd2 /= (Area*2.0);
  
  CVector3 Gu0 = Cross( Gd1, Gd2 ); Gu0 /= Dot(Gu0,Gd0);
  CVector3 Gu1 = Cross( Gd2, Gd0 ); Gu1 /= Dot(Gu1,Gd1);
  
  const CVector3 gd0( p[1][0]-p[0][0], p[1][1]-p[0][1], p[1][2]-p[0][2] );
  const CVector3 gd1( p[2][0]-p[0][0], p[2][1]-p[0][1], p[2][2]-p[0][2] );
  
  const double E2[3] = {  // green lagrange strain (with engineer's notation)
    0.5*( Dot(gd0,gd0) - Dot(Gd0,Gd0) ),
    0.5*( Dot(gd1,gd1) - Dot(Gd1,Gd1) ),
    1.0*( Dot(gd0,gd1) - Dot(Gd0,Gd1) ) };
  
  const double GuGu_xx[3] = { Gu0.x*Gu0.x, Gu1.x*Gu1.x, Gu0.x*Gu1.x };
  const double GuGu_yy[3] = { Gu0.y*Gu0.y, Gu1.y*Gu1.y, Gu0.y*Gu1.y };
  const double GuGu_xy[3] = { 2.0*Gu0.x*Gu0.y, 2.0*Gu1.x*Gu1.y, Gu0.x*Gu1.y+Gu0.y*Gu1.x };
  C[0] = E2[0]*GuGu_xx[0] + E2[1]*GuGu_xx[1] + E2[2]*GuGu_xx[2];
  C[1] = E2[0]*GuGu_yy[0] + E2[1]*GuGu_yy[1] + E2[2]*GuGu_yy[2];
  C[2] = E2[0]*GuGu_xy[0] + E2[1]*GuGu_xy[1] + E2[2]*GuGu_xy[2];
  
  const CVector3 dC0dp0 = -(GuGu_xx[0]+GuGu_xx[2])*gd0 - (GuGu_xx[1]+GuGu_xx[2])*gd1;
  const CVector3 dC0dp1 = GuGu_xx[0]*gd0 + GuGu_xx[2]*gd1;
  const CVector3 dC0dp2 = GuGu_xx[1]*gd1 + GuGu_xx[2]*gd0;
  const CVector3 dC1dp0 = -(GuGu_yy[0]+GuGu_yy[2])*gd0 - (GuGu_yy[1]+GuGu_yy[2])*gd1;
  const CVector3 dC1dp1 = GuGu_yy[0]*gd0 + GuGu_yy[2]*gd1;
  const CVector3 dC1dp2 = GuGu_yy[1]*gd1 + GuGu_yy[2]*gd0;
  const CVector3 dC2dp0 = -(GuGu_xy[0]+GuGu_xy[2])*gd0 - (GuGu_xy[1]+GuGu_xy[2])*gd1;
  const CVector3 dC2dp1 = GuGu_xy[0]*gd0 + GuGu_xy[2]*gd1;
  const CVector3 dC2dp2 = GuGu_xy[1]*gd1 + GuGu_xy[2]*gd0;
  
  dCdp[0][0*3+0] = dC0dp0.x; dCdp[0][0*3+1] = dC0dp0.y;  dCdp[0][0*3+2] = dC0dp0.z;
  dCdp[0][1*3+0] = dC0dp1.x; dCdp[0][1*3+1] = dC0dp1.y;  dCdp[0][1*3+2] = dC0dp1.z;
  dCdp[0][2*3+0] = dC0dp2.x; dCdp[0][2*3+1] = dC0dp2.y;  dCdp[0][2*3+2] = dC0dp2.z;
  dCdp[1][0*3+0] = dC1dp0.x; dCdp[1][0*3+1] = dC1dp0.y;  dCdp[1][0*3+2] = dC1dp0.z;
  dCdp[1][1*3+0] = dC1dp1.x; dCdp[1][1*3+1] = dC1dp1.y;  dCdp[1][1*3+2] = dC1dp1.z;
  dCdp[1][2*3+0] = dC1dp2.x; dCdp[1][2*3+1] = dC1dp2.y;  dCdp[1][2*3+2] = dC1dp2.z;
  dCdp[2][0*3+0] = dC2dp0.x; dCdp[2][0*3+1] = dC2dp0.y;  dCdp[2][0*3+2] = dC2dp0.z;
  dCdp[2][1*3+0] = dC2dp1.x; dCdp[2][1*3+1] = dC2dp1.y;  dCdp[2][1*3+2] = dC2dp1.z;
  dCdp[2][2*3+0] = dC2dp2.x; dCdp[2][2*3+1] = dC2dp2.y;  dCdp[2][2*3+2] = dC2dp2.z;
}

void Check_ConstraintProjection_Strain
(const double P[3][2], // (in) undeformed triangle vertex positions
 const double p[3][3] // (in) deformed triangle vertex positions)
)
{
  double C[3], dCdp[3][9];
  PBD_ConstraintProjection_Strain(C, dCdp, P, p);
  for(int ine=0;ine<3;++ine){
    for(int idim=0;idim<3;++idim){
      double eps = 1.0e-10;
      double p1[3][3]; for(int i=0;i<9;++i){ (&p1[0][0])[i] = (&p[0][0])[i]; }
      p1[ine][idim] += eps;
      double C1[3], dCdp1[3][9];
      PBD_ConstraintProjection_Strain(C1, dCdp1, P, p1);
      double diff0 = (C1[0]-C[0])/eps-dCdp[0][ine*3+idim];
      double diff1 = (C1[1]-C[1])/eps-dCdp[1][ine*3+idim];
      double diff2 = (C1[2]-C[2])/eps-dCdp[2][ine*3+idim];
      diff0 = abs(diff0);
      diff1 = abs(diff1);
      diff2 = abs(diff2);
      std::cout << diff0 << " " << diff1 << " " << diff2 << std::endl;
    }
  }
}



void PBD_ConstraintProjection_DistanceTri2D3D
(double C[3],
 double dCdp[3][9],
 const double P[3][2], // (in) undeformed triangle vertex positions
 const double p[3][3] // (in) deformed triangle vertex positions
)
{
  const double L12 = Distance2D(P[1],P[2]);
  const double L20 = Distance2D(P[2],P[0]);
  const double L01 = Distance2D(P[0],P[1]);
  CVector3 v12(p[1][0]-p[2][0], p[1][1]-p[2][1], p[1][2]-p[2][2]);
  CVector3 v20(p[2][0]-p[0][0], p[2][1]-p[0][1], p[2][2]-p[0][2]);
  CVector3 v01(p[0][0]-p[1][0], p[0][1]-p[1][1], p[0][2]-p[1][2]);
  const double l12 = v12.Length();
  const double l20 = v20.Length();
  const double l01 = v01.Length();
  C[0] = l12-L12;
  C[1] = l20-L20;
  C[2] = l01-L01;
  v12 /= l12;
  v20 /= l20;
  v01 /= l01;
  for(int i=0;i<27;++i){ (&dCdp[0][0])[i] = 0.0; }
  v12.CopyValueTo(dCdp[0]+3*1);  v12.CopyValueToScale(dCdp[0]+3*2,-1.0);
  v20.CopyValueTo(dCdp[1]+3*2);  v20.CopyValueToScale(dCdp[1]+3*0,-1.0);
  v01.CopyValueTo(dCdp[2]+3*0);  v01.CopyValueToScale(dCdp[2]+3*1,-1.0);
}


void Check_ConstraintProjection_DistanceTri2D3D
(const double P[3][2], // (in) undeformed triangle vertex positions
 const double p[3][3] // (in) deformed triangle vertex positions)
)
{
  double C[3], dCdp[3][9];
  PBD_ConstraintProjection_DistanceTri2D3D(C, dCdp, P, p);
  for(int ine=0;ine<3;++ine){
    for(int idim=0;idim<3;++idim){
      double eps = 1.0e-6;
      double p1[3][3]; for(int i=0;i<9;++i){ (&p1[0][0])[i] = (&p[0][0])[i]; }
      p1[ine][idim] += eps;
      double C1[3], dCdp1[3][9];
      PBD_ConstraintProjection_DistanceTri2D3D(C1, dCdp1, P, p1);
      double diff0 = (C1[0]-C[0])/eps-dCdp[0][ine*3+idim];
      double diff1 = (C1[1]-C[1])/eps-dCdp[1][ine*3+idim];
      double diff2 = (C1[2]-C[2])/eps-dCdp[2][ine*3+idim];
      diff0 = abs(diff0);
      diff1 = abs(diff1);
      diff2 = abs(diff2);
      std::cout << ine << " " << idim << "   " << diff0 << " " << diff1 << " " << diff2 << std::endl;
    }
  }
}

void PBD_ConstraintProjection_EnergyStVK
(double& C, // (out) energy
 double dCdp[9], // (out) 1st derivative of energy
 ////
 const double P[3][2], // (in) undeformed triangle vertex positions
 const double p[3][3], // (in) deformed triangle vertex positions
 const double lambda, // (in) Lame's 1st parameter
 const double myu)     // (in) Lame's 2nd parameter
{
  
  const CVector3 Gd0( P[1][0]-P[0][0], P[1][1]-P[0][1], 0.0 );
  const CVector3 Gd1( P[2][0]-P[0][0], P[2][1]-P[0][1], 0.0 );
  CVector3 Gd2 = Cross(Gd0, Gd1);
  const double Area = Gd2.Length()*0.5;
  Gd2 /= (Area*2.0);
  
  CVector3 Gu0 = Cross( Gd1, Gd2 ); Gu0 /= Dot(Gu0,Gd0);
  CVector3 Gu1 = Cross( Gd2, Gd0 ); Gu1 /= Dot(Gu1,Gd1);
  
  const CVector3 gd0( p[1][0]-p[0][0], p[1][1]-p[0][1], p[1][2]-p[0][2] );
  const CVector3 gd1( p[2][0]-p[0][0], p[2][1]-p[0][1], p[2][2]-p[0][2] );
  
  const double E2[3] = {  // green lagrange strain (with engineer's notation)
    0.5*( Dot(gd0,gd0) - Dot(Gd0,Gd0) ),
    0.5*( Dot(gd1,gd1) - Dot(Gd1,Gd1) ),
    1.0*( Dot(gd0,gd1) - Dot(Gd0,Gd1) ) };
  
  /*
  double Gd[3][3] = { // undeformed edge vector
    { C[1][0]-C[0][0], C[1][1]-C[0][1], C[1][2]-C[0][2] },
    { C[2][0]-C[0][0], C[2][1]-C[0][1], C[2][2]-C[0][2] }, { 0,0,0 } };
  double Area;
  UnitNormalAreaTri3D(Gd[2], Area, C[0], C[1], C[2]);
  
  double Gu[2][3]; // inverse of Gd
  {
    Cross3D(Gu[0], Gd[1], Gd[2]);
    const double invtmp1 = 1.0/Dot3D(Gu[0],Gd[0]);
    Gu[0][0] *= invtmp1;  Gu[0][1] *= invtmp1;  Gu[0][2] *= invtmp1;
    ////
    Cross3D(Gu[1], Gd[2], Gd[0]);
    const double invtmp2 = 1.0/Dot3D(Gu[1],Gd[1]);
    Gu[1][0] *= invtmp2;  Gu[1][1] *= invtmp2;  Gu[1][2] *= invtmp2;
  }
  
  const double gd[2][3] = { // deformed edge vector
    { c[1][0]-c[0][0], c[1][1]-c[0][1], c[1][2]-c[0][2] },
    { c[2][0]-c[0][0], c[2][1]-c[0][1], c[2][2]-c[0][2] } };
  
  const double E2[3] = {  // green lagrange strain (with engineer's notation)
    0.5*( Dot3D(gd[0],gd[0]) - Dot3D(Gd[0],Gd[0]) ),
    0.5*( Dot3D(gd[1],gd[1]) - Dot3D(Gd[1],Gd[1]) ),
    1.0*( Dot3D(gd[0],gd[1]) - Dot3D(Gd[0],Gd[1]) ) };
   */
  const double GuGu2[3] = { Gu0*Gu0, Gu1*Gu1, Gu1*Gu0 };
  const double Cons2[3][3] = { // constitutive tensor
    { lambda*GuGu2[0]*GuGu2[0] + 2*myu*(GuGu2[0]*GuGu2[0]),
      lambda*GuGu2[0]*GuGu2[1] + 2*myu*(GuGu2[2]*GuGu2[2]),
      lambda*GuGu2[0]*GuGu2[2] + 2*myu*(GuGu2[0]*GuGu2[2]) },
    { lambda*GuGu2[1]*GuGu2[0] + 2*myu*(GuGu2[2]*GuGu2[2]),
      lambda*GuGu2[1]*GuGu2[1] + 2*myu*(GuGu2[1]*GuGu2[1]),
      lambda*GuGu2[1]*GuGu2[2] + 2*myu*(GuGu2[2]*GuGu2[1]) },
    { lambda*GuGu2[2]*GuGu2[0] + 2*myu*(GuGu2[0]*GuGu2[2]),
      lambda*GuGu2[2]*GuGu2[1] + 2*myu*(GuGu2[2]*GuGu2[1]),
      lambda*GuGu2[2]*GuGu2[2] + 1*myu*(GuGu2[0]*GuGu2[1] + GuGu2[2]*GuGu2[2]) } };
  const double S2[3] = {  // 2nd Piola-Kirchhoff stress
    Cons2[0][0]*E2[0] + Cons2[0][1]*E2[1] + Cons2[0][2]*E2[2],
    Cons2[1][0]*E2[0] + Cons2[1][1]*E2[1] + Cons2[1][2]*E2[2],
    Cons2[2][0]*E2[0] + Cons2[2][1]*E2[1] + Cons2[2][2]*E2[2] };
  
  // compute energy
  C = 0.5*Area*(E2[0]*S2[0] + E2[1]*S2[1] + E2[2]*S2[2]);
  
  // compute 1st derivative
  const double dNdr[3][2] = { {-1.0, -1.0}, {+1.0, +0.0}, {+0.0, +1.0} };
  for(int ino=0;ino<3;ino++){
    dCdp[ino*3+0] = Area*(+S2[0]*gd0[0]*dNdr[ino][0] + S2[2]*gd0[0]*dNdr[ino][1] + S2[2]*gd1[0]*dNdr[ino][0] + S2[1]*gd1[0]*dNdr[ino][1]);
    dCdp[ino*3+1] = Area*(+S2[0]*gd0[1]*dNdr[ino][0] + S2[2]*gd0[1]*dNdr[ino][1] + S2[2]*gd1[1]*dNdr[ino][0] + S2[1]*gd1[1]*dNdr[ino][1]);
    dCdp[ino*3+2] = Area*(+S2[0]*gd0[2]*dNdr[ino][0] + S2[2]*gd0[2]*dNdr[ino][1] + S2[2]*gd1[2]*dNdr[ino][0] + S2[1]*gd1[2]*dNdr[ino][1]);
  }
}


void Check_ConstraintProjection_EnergyStVK
(const double P[3][2], // (in) undeformed triangle vertex positions
 const double p[3][3], // (in) deformed triangle vertex positions)
 const double lambda,
 const double myu
)
{
  double C, dCdp[9];
  PBD_ConstraintProjection_EnergyStVK(C, dCdp, P, p, lambda, myu);
  for(int ine=0;ine<3;++ine){
    for(int idim=0;idim<3;++idim){
      double eps = 1.0e-8;
      double p1[3][3]; for(int i=0;i<9;++i){ (&p1[0][0])[i] = (&p[0][0])[i]; }
      p1[ine][idim] += eps;
      double C1, dCdp1[9];
      PBD_ConstraintProjection_EnergyStVK(C1, dCdp1, P, p1, lambda, myu);
      double diff0 = (C1-C)/eps-dCdp[ine*3+idim];
      diff0 = abs(diff0);
      std::cout << ine << " " << idim << "   " << diff0 << std::endl;
    }
  }
}


void PBD_ConstraintProjection_DistanceTet
(double C[6],
 double dCdp[6][12],
 const double P[4][3], // (in) undeformed triangle vertex positions
 const double p[4][3] // (in) deformed triangle vertex positions
)
{
  const double L01 = Distance3D(P[0],P[1]);
  const double L02 = Distance3D(P[0],P[2]);
  const double L03 = Distance3D(P[0],P[3]);
  const double L12 = Distance3D(P[1],P[2]);
  const double L13 = Distance3D(P[1],P[3]);
  const double L23 = Distance3D(P[2],P[3]);
  CVector3 v01(p[0][0]-p[1][0], p[0][1]-p[1][1], p[0][2]-p[1][2]);
  CVector3 v02(p[0][0]-p[2][0], p[0][1]-p[2][1], p[0][2]-p[2][2]);
  CVector3 v03(p[0][0]-p[3][0], p[0][1]-p[3][1], p[0][2]-p[3][2]);
  CVector3 v12(p[1][0]-p[2][0], p[1][1]-p[2][1], p[1][2]-p[2][2]);
  CVector3 v13(p[1][0]-p[3][0], p[1][1]-p[3][1], p[1][2]-p[3][2]);
  CVector3 v23(p[2][0]-p[3][0], p[2][1]-p[3][1], p[2][2]-p[3][2]);
  const double l01 = v01.Length();
  const double l02 = v02.Length();
  const double l03 = v03.Length();
  const double l12 = v12.Length();
  const double l13 = v13.Length();
  const double l23 = v23.Length();
  C[0] = l01-L01;
  C[1] = l02-L02;
  C[2] = l03-L03;
  C[3] = l12-L12;
  C[4] = l13-L13;
  C[5] = l23-L23;
  ////
  v01 /= l01;
  v02 /= l02;
  v03 /= l03;
  v12 /= l12;
  v13 /= l13;
  v23 /= l23;
  ////
  for(int i=0;i<6*3*4;++i){ (&dCdp[0][0])[i] = 0.0; }
  v01.CopyValueTo(dCdp[0]+0*3);  v01.CopyValueToScale(dCdp[0]+1*3,-1.0);
  v02.CopyValueTo(dCdp[1]+0*3);  v02.CopyValueToScale(dCdp[1]+2*3,-1.0);
  v03.CopyValueTo(dCdp[2]+0*3);  v03.CopyValueToScale(dCdp[2]+3*3,-1.0);
  v12.CopyValueTo(dCdp[3]+1*3);  v12.CopyValueToScale(dCdp[3]+2*3,-1.0);
  v13.CopyValueTo(dCdp[4]+1*3);  v13.CopyValueToScale(dCdp[4]+3*3,-1.0);
  v23.CopyValueTo(dCdp[5]+2*3);  v23.CopyValueToScale(dCdp[5]+3*3,-1.0);
}
