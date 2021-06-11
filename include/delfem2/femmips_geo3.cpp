/*
 * Copyright (c) 2019 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */

#include "delfem2/geo3_v23m34q.h"
#include "delfem2/femmips_geo3.h"

// =======================================

DFM2_INLINE void delfem2::WdWddW_MIPS
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
  CVec3d p0(c[0][0],c[0][1],c[0][2]);
  CVec3d p1(c[1][0],c[1][1],c[1][2]);
  CVec3d p2(c[2][0],c[2][1],c[2][2]);
  CVec3d P0(C[0][0],C[0][1],C[0][2]);
  CVec3d P1(C[1][0],C[1][1],C[1][2]);
  CVec3d P2(C[2][0],C[2][1],C[2][2]);
  CVec3d v01 = p1-p0;
  CVec3d v12 = p2-p1;
  CVec3d v20 = p0-p2;
  CVec3d n = v01^v20;
  double area = n.norm()*0.5;
  double Area = ((P1-P0)^(P2-P0)).norm()*0.5;
  double la = (p1-p2).dot(p1-p2);
  double lb = (p2-p0).dot(p2-p0);
  double lc = (p0-p1).dot(p0-p1);
  double lA = (P1-P2).dot(P1-P2);
  double lB = (P2-P0).dot(P2-P0);
  double lC = (P0-P1).dot(P0-P1);
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
  CVec3d dECd0 = t00*p0 + t01*p1 + t02*p2;
  CVec3d dECd1 = t01*p0 + t11*p1 + t12*p2;
  CVec3d dECd2 = t02*p0 + t12*p1 + t22*p2;
  ////
  dE[0][0]=dECd0.x; dE[0][1]=dECd0.y; dE[0][2]=dECd0.z;
  dE[1][0]=dECd1.x; dE[1][1]=dECd1.y; dE[1][2]=dECd1.z;
  dE[2][0]=dECd2.x; dE[2][1]=dECd2.y; dE[2][2]=dECd2.z;
  
  CMat3d (*op)(const CVec3d&, const CVec3d&) = Mat3_OuterProduct;
  
  double tmp1 = 0.25/area;
  CVec3d dad0 = ((v20.dot(v12))*v01-(v01.dot(v12))*v20)*tmp1;
  CVec3d dad1 = ((v01.dot(v20))*v12-(v12.dot(v20))*v01)*tmp1;
  CVec3d dad2 = ((v12.dot(v01))*v20-(v20.dot(v01))*v12)*tmp1;
  CMat3d ddad0d0 = (CMat3d::Identity(v12.dot(v12)) - op(v12,v12)                   - 4*op(dad0,dad0))*tmp1;
  CMat3d ddad0d1 = (CMat3d::Identity(v20.dot(v12)) - op(v20,v12-v01) - op(v01,v20) - 4*op(dad0,dad1))*tmp1;
  CMat3d ddad0d2 = (CMat3d::Identity(v01.dot(v12)) - op(v01,v12-v20) - op(v20,v01) - 4*op(dad0,dad2))*tmp1;
  CMat3d ddad1d0 = (CMat3d::Identity(v12.dot(v20)) - op(v12,v20-v01) - op(v01,v12) - 4*op(dad1,dad0))*tmp1;
  CMat3d ddad1d1 = (CMat3d::Identity(v20.dot(v20)) - op(v20,v20)                   - 4*op(dad1,dad1))*tmp1;
  CMat3d ddad1d2 = (CMat3d::Identity(v01.dot(v20)) - op(v01,v20-v12) - op(v12,v01) - 4*op(dad1,dad2))*tmp1;
  CMat3d ddad2d0 = (CMat3d::Identity(v12.dot(v01)) - op(v12,v01-v20) - op(v20,v12) - 4*op(dad2,dad0))*tmp1;
  CMat3d ddad2d1 = (CMat3d::Identity(v20.dot(v01)) - op(v20,v01-v12) - op(v12,v20) - 4*op(dad2,dad1))*tmp1;
  CMat3d ddad2d2 = (CMat3d::Identity(v01.dot(v01)) - op(v01,v01)                   - 4*op(dad2,dad2))*tmp1;
  
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
  CMat3d ddEd0d0 = EC*dEA*ddad0d0 + EC*ddEA*op(dad0,dad0) + EA*CMat3d::Identity(t00) + dEA*op(dad0,dECd0)*2.0;
  CMat3d ddEd0d1 = EC*dEA*ddad0d1 + EC*ddEA*op(dad0,dad1) + EA*CMat3d::Identity(t01) + dEA*op(dad0,dECd1) + dEA*op(dad1,dECd0);
  CMat3d ddEd0d2 = EC*dEA*ddad0d2 + EC*ddEA*op(dad0,dad2) + EA*CMat3d::Identity(t02) + dEA*op(dad0,dECd2) + dEA*op(dad2,dECd0);
  CMat3d ddEd1d0 = EC*dEA*ddad1d0 + EC*ddEA*op(dad1,dad0) + EA*CMat3d::Identity(t01) + dEA*op(dad1,dECd0) + dEA*op(dad0,dECd1);
  CMat3d ddEd1d1 = EC*dEA*ddad1d1 + EC*ddEA*op(dad1,dad1) + EA*CMat3d::Identity(t11) + dEA*op(dad1,dECd1)*2.0;
  CMat3d ddEd1d2 = EC*dEA*ddad1d2 + EC*ddEA*op(dad1,dad2) + EA*CMat3d::Identity(t12) + dEA*op(dad1,dECd2) + dEA*op(dad2,dECd1);
  CMat3d ddEd2d0 = EC*dEA*ddad2d0 + EC*ddEA*op(dad2,dad0) + EA*CMat3d::Identity(t02) + dEA*op(dad2,dECd0) + dEA*op(dad0,dECd2);
  CMat3d ddEd2d1 = EC*dEA*ddad2d1 + EC*ddEA*op(dad2,dad1) + EA*CMat3d::Identity(t12) + dEA*op(dad2,dECd1) + dEA*op(dad1,dECd2);
  CMat3d ddEd2d2 = EC*dEA*ddad2d2 + EC*ddEA*op(dad2,dad2) + EA*CMat3d::Identity(t22) + dEA*op(dad2,dECd2)*2.0;
  
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
