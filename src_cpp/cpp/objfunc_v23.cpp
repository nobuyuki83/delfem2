/*
 * Copyright (c) 2019 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */

#include <stdio.h>

#include "delfem2/vec2.h"
#include "delfem2/vec3.h"
#include "delfem2/mat3.h"

#include "delfem2/objfunc_v23.h"


void PBD_Post
(std::vector<double>& aXYZ,
 std::vector<double>& aUVW,
 double dt,
 const std::vector<double>& aXYZt,
 const std::vector<int>& aBCFlag)
{
  const unsigned int ndof = aXYZt.size();
  assert( aBCFlag.size() == ndof );
  assert( aXYZ.size() == ndof );
  assert( aUVW.size() == ndof );
  for(unsigned int idof=0;idof<ndof;++idof){
    if( aBCFlag[idof] != 0 ) continue;
    aUVW[idof] = (aXYZt[idof]-aXYZ[idof])/dt;
    aXYZ[idof] = aXYZt[idof];
  }
}

void PBD_Pre3D
(std::vector<double>& aXYZt,
 double dt,
 const double gravity[3],
 const std::vector<double>& aXYZ,
 const std::vector<double>& aUVW,
 const std::vector<int>& aBCFlag)
{
  const unsigned int np = aXYZ.size()/3;
  assert( aUVW.size() == np*3 );
  assert( aBCFlag.size() == np*3 );
  aXYZt.resize(np*3);
  for(unsigned int ip=0;ip<np;++ip){
    for(int idim=0;idim<3;++idim){
      const int idof = ip*3+idim;
      if( aBCFlag[idof] != 0 ){ aXYZt[idof] = aXYZ[idof]; continue; }
      aXYZt[idof] = aXYZ[idof]+dt*aUVW[idof]+dt*dt*gravity[idim];
    }
  }
}


void PBD_Update_Const3
(double* aXYZt,
 const int np,
 const int ndim,
 const double* m,
 const double* C,
 const double* dCdp,
 const int* aIP)
{
  double mi[np];
  for(int ip=0;ip<np;++ip){ mi[ip] = 1.0/m[ip]; }
  /////
  const int nc = 3;
  double MinvC[nc*np*ndim];
  for(int ic=0;ic<nc;++ic){
    for(int ine=0;ine<np;++ine){
      MinvC[ic*np*ndim+ine*3+0] = dCdp[ic*np*ndim+ine*ndim+0]*mi[ine];
      MinvC[ic*np*ndim+ine*3+1] = dCdp[ic*np*ndim+ine*ndim+1]*mi[ine];
      MinvC[ic*np*ndim+ine*3+2] = dCdp[ic*np*ndim+ine*ndim+2]*mi[ine];
    }
  }
  double A[nc*nc] = {0,0,0, 0,0,0, 0,0,0};
  for(int i=0;i<np*ndim;++i){
    A[0*3+0] += MinvC[0*np*ndim+i]*dCdp[0*np*ndim+i];
    A[0*3+1] += MinvC[0*np*ndim+i]*dCdp[1*np*ndim+i];
    A[0*3+2] += MinvC[0*np*ndim+i]*dCdp[2*np*ndim+i];
    A[1*3+0] += MinvC[1*np*ndim+i]*dCdp[0*np*ndim+i];
    A[1*3+1] += MinvC[1*np*ndim+i]*dCdp[1*np*ndim+i];
    A[1*3+2] += MinvC[1*np*ndim+i]*dCdp[2*np*ndim+i];
    A[2*3+0] += MinvC[2*np*ndim+i]*dCdp[0*np*ndim+i];
    A[2*3+1] += MinvC[2*np*ndim+i]*dCdp[1*np*ndim+i];
    A[2*3+2] += MinvC[2*np*ndim+i]*dCdp[2*np*ndim+i];
  }
  double Ainv[nc*nc]; InverseMat3(Ainv, A);
  double lmd[nc]; MatVec3(lmd, Ainv, C);
  for(int ine=0;ine<np;++ine){
    const int ip0 = aIP[ine];
    for(int ic=0;ic<nc;++ic){
      for(int idim=0;idim<ndim;++idim){
        aXYZt[ip0*3+idim] -= MinvC[ic*np*ndim+ine*ndim+idim]*lmd[ic];
      }
    }

  }
}

void PBD_ConstProj_Rigid3D
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

void PBD_ConstProj_Rigid2D
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

void PBD_CdC_TriStrain2D3D
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
  
  dC0dp0.CopyValueTo(dCdp[0]+0*3);
  dC0dp1.CopyValueTo(dCdp[0]+1*3);
  dC0dp2.CopyValueTo(dCdp[0]+2*3);
  dC1dp0.CopyValueTo(dCdp[1]+0*3);
  dC1dp1.CopyValueTo(dCdp[1]+1*3);
  dC1dp2.CopyValueTo(dCdp[1]+2*3);
  dC2dp0.CopyValueTo(dCdp[2]+0*3);
  dC2dp1.CopyValueTo(dCdp[2]+1*3);
  dC2dp2.CopyValueTo(dCdp[2]+2*3);
}

void Check_CdC_TriStrain
(const double P[3][2], // (in) undeformed triangle vertex positions
 const double p[3][3] // (in) deformed triangle vertex positions)
)
{
  double C[3], dCdp[3][9];
  PBD_CdC_TriStrain2D3D(C, dCdp, P, p);
  for(int ine=0;ine<3;++ine){
    for(int idim=0;idim<3;++idim){
      double eps = 1.0e-10;
      double p1[3][3]; for(int i=0;i<9;++i){ (&p1[0][0])[i] = (&p[0][0])[i]; }
      p1[ine][idim] += eps;
      double C1[3], dCdp1[3][9];
      PBD_CdC_TriStrain2D3D(C1, dCdp1, P, p1);
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


void PBD_CdC_QuadBend
(double C[3],
 double dCdp[3][12],
 const double P[4][3],
 const double p[4][3])
{
  const double A0 = TriArea3D(P[0],P[2],P[3]);
  const double A1 = TriArea3D(P[1],P[3],P[2]);
  const double L = Distance3D(P[2],P[3]);
  const double H0 = 2.0*A0/L;
  const double H1 = 2.0*A1/L;
  const CVector3 e23(P[3][0]-P[2][0], P[3][1]-P[2][1], P[3][2]-P[2][2]);
  const CVector3 e02(P[2][0]-P[0][0], P[2][1]-P[0][1], P[2][2]-P[0][2]);
  const CVector3 e03(P[3][0]-P[0][0], P[3][1]-P[0][1], P[3][2]-P[0][2]);
  const CVector3 e12(P[2][0]-P[1][0], P[2][1]-P[1][1], P[2][2]-P[1][2]);
  const CVector3 e13(P[3][0]-P[1][0], P[3][1]-P[1][1], P[3][2]-P[1][2]);
  const double cot023 = -(e02*e23)/H0;
  const double cot032 = +(e03*e23)/H0;
  const double cot123 = -(e12*e23)/H1;
  const double cot132 = +(e13*e23)/H1;
  const double tmp0 = sqrt(3.0)/(sqrt(A0+A1)*L);
  const double K[4] = {
    (-cot023-cot032)*tmp0,
    (-cot123-cot132)*tmp0,
    (+cot032+cot132)*tmp0,
    (+cot023+cot123)*tmp0 };
  C[0] = K[0]*p[0][0] + K[1]*p[1][0] + K[2]*p[2][0] + K[3]*p[3][0];
  C[1] = K[0]*p[0][1] + K[1]*p[1][1] + K[2]*p[2][1] + K[3]*p[3][1];
  C[2] = K[0]*p[0][2] + K[1]*p[1][2] + K[2]*p[2][2] + K[3]*p[3][2];
  for(int i=0;i<36;++i){ (&dCdp[0][0])[i] = 0.0; }
  for(int idim=0;idim<3;++idim){
    dCdp[idim][0*3+idim] = K[0];
    dCdp[idim][1*3+idim] = K[1];
    dCdp[idim][2*3+idim] = K[2];
    dCdp[idim][3*3+idim] = K[3];
  }
}
