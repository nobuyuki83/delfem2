/*
 * Copyright (c) 2019 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */

#include "delfem2/femcloth.h"

#include <cmath>


// --------------------------------------------------------

DFM2_INLINE void MakeConstMatrix3D(
    double C[6][6],
    double lambda,
    double myu,
    const double Gu[3][3])
{
  const double GuGu2[6] = {
    delfem2::femutil::Dot3D(Gu[0],Gu[0]), // 0 xx
    delfem2::femutil::Dot3D(Gu[1],Gu[1]), // 1 yy
    delfem2::femutil::Dot3D(Gu[2],Gu[2]), // 2 zz
    delfem2::femutil::Dot3D(Gu[0],Gu[1]), // 3 xy
    delfem2::femutil::Dot3D(Gu[1],Gu[2]), // 4 yz
    delfem2::femutil::Dot3D(Gu[2],Gu[0])  // 5 zx
  };
  C[0][0] = lambda*GuGu2[0]*GuGu2[0] + 2*myu*(GuGu2[0]*GuGu2[0]); // 00(0):00(0) 00(0):00(0)
  C[0][1] = lambda*GuGu2[0]*GuGu2[1] + 2*myu*(GuGu2[3]*GuGu2[3]); // 00(0):11(1) 01(3):01(3)
  C[0][2] = lambda*GuGu2[0]*GuGu2[2] + 2*myu*(GuGu2[5]*GuGu2[5]); // 00(0):22(2) 02(5):02(5)
  C[0][3] = lambda*GuGu2[0]*GuGu2[3] + 2*myu*(GuGu2[0]*GuGu2[3]); // 00(0):01(3) 00(0):01(3)
  C[0][4] = lambda*GuGu2[0]*GuGu2[4] + 2*myu*(GuGu2[3]*GuGu2[5]); // 00(0):12(4) 01(3):02(5)
  C[0][5] = lambda*GuGu2[0]*GuGu2[5] + 2*myu*(GuGu2[0]*GuGu2[5]); // 00(0):20(5) 00(0):02(5)
  C[1][0] = lambda*GuGu2[1]*GuGu2[0] + 2*myu*(GuGu2[3]*GuGu2[3]); // 11(1):00(0) 01(3):01(3)
  C[1][1] = lambda*GuGu2[1]*GuGu2[1] + 2*myu*(GuGu2[1]*GuGu2[1]); // 11(1):11(1) 11(1):11(1)
  C[1][2] = lambda*GuGu2[1]*GuGu2[2] + 2*myu*(GuGu2[4]*GuGu2[4]); // 11(1):22(2) 12(4):12(4)
  C[1][3] = lambda*GuGu2[1]*GuGu2[3] + 2*myu*(GuGu2[1]*GuGu2[3]); // 11(1):01(3) 11(1):01(3)
  C[1][4] = lambda*GuGu2[1]*GuGu2[4] + 2*myu*(GuGu2[1]*GuGu2[4]); // 11(1):12(4) 11(1):12(4)
  C[1][5] = lambda*GuGu2[1]*GuGu2[5] + 2*myu*(GuGu2[3]*GuGu2[4]); // 11(1):20(5) 12(4):10(3)
  C[2][0] = lambda*GuGu2[2]*GuGu2[0] + 2*myu*(GuGu2[5]*GuGu2[5]); // 22(2):00(0) 02(5):02(5)
  C[2][1] = lambda*GuGu2[2]*GuGu2[1] + 2*myu*(GuGu2[4]*GuGu2[4]); // 22(2):11(1) 12(4):12(4)
  C[2][2] = lambda*GuGu2[2]*GuGu2[2] + 2*myu*(GuGu2[2]*GuGu2[2]); // 22(2):22(2) 22(2):22(2)
  C[2][3] = lambda*GuGu2[2]*GuGu2[3] + 2*myu*(GuGu2[4]*GuGu2[5]); // 22(2):01(3) 12(4):02(5)
  C[2][4] = lambda*GuGu2[2]*GuGu2[4] + 2*myu*(GuGu2[2]*GuGu2[4]); // 22(2):12(4) 22(2):12(4)
  C[2][5] = lambda*GuGu2[2]*GuGu2[5] + 2*myu*(GuGu2[2]*GuGu2[5]); // 22(2):02(5) 22(2):02(5)
  C[3][0] = lambda*GuGu2[3]*GuGu2[0] + 2*myu*(GuGu2[3]*GuGu2[0]); // 01(3):00(0) 00(0):01(3)
  C[3][1] = lambda*GuGu2[3]*GuGu2[1] + 2*myu*(GuGu2[3]*GuGu2[1]); // 01(3):11(1) 11(1):01(3)
  C[3][2] = lambda*GuGu2[3]*GuGu2[2] + 2*myu*(GuGu2[4]*GuGu2[5]); // 01(3):22(2) 12(4):02(5)
  C[3][3] = lambda*GuGu2[3]*GuGu2[3] + 1*myu*(GuGu2[0]*GuGu2[1] + GuGu2[3]*GuGu2[3]); // 01(3):01(3) 00(0):11(1) 01(3):01(3)
  C[3][4] = lambda*GuGu2[3]*GuGu2[4] + 1*myu*(GuGu2[3]*GuGu2[4] + GuGu2[1]*GuGu2[5]); // 01(3):12(4) 01(3):12(4) 11(1):02(5)
  C[3][5] = lambda*GuGu2[3]*GuGu2[5] + 1*myu*(GuGu2[5]*GuGu2[3] + GuGu2[0]*GuGu2[4]); // 01(3):20(5) 02(5):10(3) 00(0):12(4)
  C[4][0] = lambda*GuGu2[4]*GuGu2[0] + 2*myu*(GuGu2[3]*GuGu2[5]); // 12(4):00(0) 01(3):02(5)
  C[4][1] = lambda*GuGu2[4]*GuGu2[1] + 2*myu*(GuGu2[1]*GuGu2[4]); // 12(4):11(1) 11(1):12(4)
  C[4][2] = lambda*GuGu2[4]*GuGu2[2] + 2*myu*(GuGu2[2]*GuGu2[4]); // 12(4):22(2) 22(2):12(4)
  C[4][3] = lambda*GuGu2[4]*GuGu2[3] + 1*myu*(GuGu2[3]*GuGu2[4] + GuGu2[1]*GuGu2[5]); // 12(4):01(3) 10(3):21(4) 11(1):20(5)
  C[4][4] = lambda*GuGu2[4]*GuGu2[4] + 1*myu*(GuGu2[1]*GuGu2[2] + GuGu2[4]*GuGu2[4]); // 12(4):12(4) 11(1):22(2) 12(4):21(4)
  C[4][5] = lambda*GuGu2[4]*GuGu2[5] + 1*myu*(GuGu2[4]*GuGu2[5] + GuGu2[3]*GuGu2[2]); // 12(4):20(5) 12(4):20(5) 10(3):22(2)
  C[5][0] = lambda*GuGu2[5]*GuGu2[0] + 2*myu*(GuGu2[0]*GuGu2[5]); // 02(5):00(0) 00(0):02(5)
  C[5][1] = lambda*GuGu2[5]*GuGu2[1] + 2*myu*(GuGu2[3]*GuGu2[4]); // 02(5):11(1) 10(3):12(4)
  C[5][2] = lambda*GuGu2[5]*GuGu2[2] + 2*myu*(GuGu2[2]*GuGu2[5]); // 02(5):22(2) 22(2):02(5)
  C[5][3] = lambda*GuGu2[5]*GuGu2[3] + 1*myu*(GuGu2[0]*GuGu2[4] + GuGu2[3]*GuGu2[5]); // 02(5):01(3) 00(0):21(4) 01(3):20(5)
  C[5][4] = lambda*GuGu2[5]*GuGu2[4] + 1*myu*(GuGu2[3]*GuGu2[2] + GuGu2[5]*GuGu2[4]); // 02(5):12(4) 01(3):22(2) 02(5):21(4)
  C[5][5] = lambda*GuGu2[5]*GuGu2[5] + 1*myu*(GuGu2[5]*GuGu2[5] + GuGu2[0]*GuGu2[2]); // 02(5):20(5) 02(5):20(5) 00(0):22(2)
}

DFM2_INLINE void MakePositiveDefinite_Sim22(
    const double s2[3],
    double s3[3])
{
  const double b = (s2[0]+s2[1])*0.5;
  const double d = (s2[0]-s2[1])*(s2[0]-s2[1])*0.25 + s2[2]*s2[2];
  const double e = sqrt(d);
  if( b-e > 1.0e-20 ){
    s3[0] = s2[0];
    s3[1] = s2[1];
    s3[2] = s2[2];
    return;
  }
  if( b+e < 0 ){
    s3[0] = 0;
    s3[1] = 0;
    s3[2] = 0;
    return;
  }
  const double l = b+e;
  double t0[2] = { s2[0]-l, s2[2]   };
  double t1[2] = { s2[2],   s2[1]-l };
  //  std::cout << t0[0]*t1[1]-t0[1]*t1[0] << std::endl;
  const double sqlen_t0 = t0[0]*t0[0]+t0[1]*t0[1];
  const double sqlen_t1 = t1[0]*t1[0]+t1[1]*t1[1];
  if( sqlen_t0 > sqlen_t1 ){
    if( sqlen_t0 < 1.0e-20 ){
      s3[0] = 0;
      s3[1] = 0;
      s3[2] = 0;
      return;
    }
    const double invlen_t0 = 1.0/sqrt(sqlen_t0);
    t0[0] *= invlen_t0;
    t0[1] *= invlen_t0;
    s3[0] = l*t0[0]*t0[0];
    s3[1] = l*t0[1]*t0[1];
    s3[2] = l*t0[0]*t0[1];
  }
  else{
    if( sqlen_t1 < 1.0e-20 ){
      s3[0] = 0;
      s3[1] = 0;
      s3[2] = 0;
      return;
    }
    const double invlen_t1 = 1.0/sqrt(sqlen_t1);
    t1[0] *= invlen_t1;
    t1[1] *= invlen_t1;
    s3[0] = l*t1[0]*t1[0];
    s3[1] = l*t1[1]*t1[1];
    s3[2] = l*t1[0]*t1[1];
  }
}

// compute energy and its 1st and 2nd derivative for cloth bending
DFM2_INLINE void delfem2::WdWddW_Bend(
    double& W,  // (out) strain energy
    double dW[4][3], // (out) 1st derivative of energy
    double ddW[4][4][3][3], // (out) 2nd derivative of energy
    //
    const double C[4][3], // (in) undeformed triangle vertex positions
    const double c[4][3], // (in) deformed triangle vertex positions
    double stiff)
{
  const double A0 = femutil::TriArea3D(C[0],C[2],C[3]);
  const double A1 = femutil::TriArea3D(C[1],C[3],C[2]);
  const double L0 = femutil::Distance3D(C[2],C[3]);
  const double H0 = A0*2.0/L0;
  const double H1 = A1*2.0/L0;
  const double e23[3] = { C[3][0]-C[2][0], C[3][1]-C[2][1], C[3][2]-C[2][2] };
  const double e02[3] = { C[2][0]-C[0][0], C[2][1]-C[0][1], C[2][2]-C[0][2] };
  const double e03[3] = { C[3][0]-C[0][0], C[3][1]-C[0][1], C[3][2]-C[0][2] };
  const double e12[3] = { C[2][0]-C[1][0], C[2][1]-C[1][1], C[2][2]-C[1][2] };
  const double e13[3] = { C[3][0]-C[1][0], C[3][1]-C[1][1], C[3][2]-C[1][2] };
  double cot023, cot032;
  {
    const double r2 = -femutil::Dot3D(e02,e23);
    const double r3 = +femutil::Dot3D(e03,e23);
    cot023 = r2/H0;
    cot032 = r3/H0;
  }
  double cot123, cot132;
  {
    const double r2 = -femutil::Dot3D(e12,e23);
    const double r3 = +femutil::Dot3D(e13,e23);
    cot123 = r2/H1;
    cot132 = r3/H1;
  }
  const double tmp0 = stiff/((A0+A1)*L0*L0);
  const double K[4] = { -cot023-cot032, -cot123-cot132, cot032+cot132, cot023+cot123 };
  
  // compute 2nd derivative of energy
  for(int i=0;i<4*4*3*3;i++){ (&ddW[0][0][0][0])[i] = 0; }
  for(int ino=0;ino<4;ino++){
    for(int jno=0;jno<4;jno++){
      const double tmp = K[ino]*K[jno]*tmp0;
      ddW[ino][jno][0][0] = tmp;
      ddW[ino][jno][1][1] = tmp;
      ddW[ino][jno][2][2] = tmp;
    }
  }
  // compute 1st derivative of energy
  W = 0.0;
  for(int ino=0;ino<4;ino++){
    for(int idim=0;idim<3;idim++){
      dW[ino][idim] = 0;
      for(int jno=0;jno<4;jno++){
        for(int jdim=0;jdim<3;jdim++){
          dW[ino][idim] += ddW[ino][jno][idim][jdim]*c[jno][jdim];
        }
      }
      W += dW[ino][idim]*c[ino][idim];
    }
  }
}

DFM2_INLINE void delfem2::WdWddW_CST(
    double& W, // (out) energy
    double dW[3][3], // (out) 1st derivative of energy
    double ddW[3][3][3][3], // (out) 2nd derivative of energy
    //
    const double C[3][3], // (in) undeformed triangle vertex positions
    const double c[3][3], // (in) deformed triangle vertex positions
    const double lambda, // (in) Lame's 1st parameter
    const double myu)     // (in) Lame's 2nd parameter
{
  double Gd[3][3] = { // undeformed edge vector
    { C[1][0]-C[0][0], C[1][1]-C[0][1], C[1][2]-C[0][2] },
    { C[2][0]-C[0][0], C[2][1]-C[0][1], C[2][2]-C[0][2] }, { 0,0,0 } };
  double Area;
  femutil::UnitNormalAreaTri3D(Gd[2], Area, C[0], C[1], C[2]);
  
  double Gu[2][3]; // inverse of Gd
  {
    femutil::Cross3D(Gu[0], Gd[1], Gd[2]);
    const double invtmp1 = 1.0/femutil::Dot3D(Gu[0],Gd[0]);
    Gu[0][0] *= invtmp1;	Gu[0][1] *= invtmp1;	Gu[0][2] *= invtmp1;
    //
    femutil::Cross3D(Gu[1], Gd[2], Gd[0]);
    const double invtmp2 = 1.0/femutil::Dot3D(Gu[1],Gd[1]);
    Gu[1][0] *= invtmp2;	Gu[1][1] *= invtmp2;	Gu[1][2] *= invtmp2;
  }
  
  const double gd[2][3] = { // deformed edge vector
    { c[1][0]-c[0][0], c[1][1]-c[0][1], c[1][2]-c[0][2] },
    { c[2][0]-c[0][0], c[2][1]-c[0][1], c[2][2]-c[0][2] } };
  
  const double E2[3] = {  // green lagrange strain (with engineer's notation)
    0.5*( femutil::Dot3D(gd[0],gd[0]) - femutil::Dot3D(Gd[0],Gd[0]) ),
    0.5*( femutil::Dot3D(gd[1],gd[1]) - femutil::Dot3D(Gd[1],Gd[1]) ),
    1.0*( femutil::Dot3D(gd[0],gd[1]) - femutil::Dot3D(Gd[0],Gd[1]) ) };
  const double GuGu2[3] = { femutil::Dot3D(Gu[0],Gu[0]), femutil::Dot3D(Gu[1],Gu[1]), femutil::Dot3D(Gu[1],Gu[0]) };
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
  W = 0.5*Area*(E2[0]*S2[0] + E2[1]*S2[1] + E2[2]*S2[2]);
  
  // compute 1st derivative
  const double dNdr[3][2] = { {-1.0, -1.0}, {+1.0, +0.0}, {+0.0, +1.0} };
  for(int ino=0;ino<3;ino++){
    for(int idim=0;idim<3;idim++){
      dW[ino][idim] = Area*
      (+S2[0]*gd[0][idim]*dNdr[ino][0]
       +S2[2]*gd[0][idim]*dNdr[ino][1]
       +S2[2]*gd[1][idim]*dNdr[ino][0]
       +S2[1]*gd[1][idim]*dNdr[ino][1]);
    }
  }
  
  double S3[3] = { S2[0], S2[1], S2[2] };
  MakePositiveDefinite_Sim22(S2,S3);
  
  // compute second derivative
  for(int ino=0;ino<3;ino++){
    for(int jno=0;jno<3;jno++){
      for(int idim=0;idim<3;idim++){
        for(int jdim=0;jdim<3;jdim++){
          double dtmp0 = 0;
          dtmp0 += gd[0][idim]*dNdr[ino][0]*Cons2[0][0]*gd[0][jdim]*dNdr[jno][0];
          dtmp0 += gd[0][idim]*dNdr[ino][0]*Cons2[0][1]*gd[1][jdim]*dNdr[jno][1];
          dtmp0 += gd[0][idim]*dNdr[ino][0]*Cons2[0][2]*gd[0][jdim]*dNdr[jno][1];
          dtmp0 += gd[0][idim]*dNdr[ino][0]*Cons2[0][2]*gd[1][jdim]*dNdr[jno][0];
          dtmp0 += gd[1][idim]*dNdr[ino][1]*Cons2[1][0]*gd[0][jdim]*dNdr[jno][0];
          dtmp0 += gd[1][idim]*dNdr[ino][1]*Cons2[1][1]*gd[1][jdim]*dNdr[jno][1];
          dtmp0 += gd[1][idim]*dNdr[ino][1]*Cons2[1][2]*gd[0][jdim]*dNdr[jno][1];
          dtmp0 += gd[1][idim]*dNdr[ino][1]*Cons2[1][2]*gd[1][jdim]*dNdr[jno][0];
          dtmp0 += gd[0][idim]*dNdr[ino][1]*Cons2[2][0]*gd[0][jdim]*dNdr[jno][0];
          dtmp0 += gd[0][idim]*dNdr[ino][1]*Cons2[2][1]*gd[1][jdim]*dNdr[jno][1];
          dtmp0 += gd[0][idim]*dNdr[ino][1]*Cons2[2][2]*gd[0][jdim]*dNdr[jno][1];
          dtmp0 += gd[0][idim]*dNdr[ino][1]*Cons2[2][2]*gd[1][jdim]*dNdr[jno][0];
          dtmp0 += gd[1][idim]*dNdr[ino][0]*Cons2[2][0]*gd[0][jdim]*dNdr[jno][0];
          dtmp0 += gd[1][idim]*dNdr[ino][0]*Cons2[2][1]*gd[1][jdim]*dNdr[jno][1];
          dtmp0 += gd[1][idim]*dNdr[ino][0]*Cons2[2][2]*gd[0][jdim]*dNdr[jno][1];
          dtmp0 += gd[1][idim]*dNdr[ino][0]*Cons2[2][2]*gd[1][jdim]*dNdr[jno][0];
          ddW[ino][jno][idim][jdim] = dtmp0*Area;
        }
      }
      const double dtmp1 = Area*
      (+S3[0]*dNdr[ino][0]*dNdr[jno][0]
       +S3[2]*dNdr[ino][0]*dNdr[jno][1]
       +S3[2]*dNdr[ino][1]*dNdr[jno][0]
       +S3[1]*dNdr[ino][1]*dNdr[jno][1]);
      ddW[ino][jno][0][0] += dtmp1;
      ddW[ino][jno][1][1] += dtmp1;
      ddW[ino][jno][2][2] += dtmp1;
    }
  }
}



// compute energy and its 1st and 2nd derivative for contact against object
DFM2_INLINE void delfem2::WdWddW_Contact(
    double& W,  // (out) energy
    double dW[3], // (out) 1st derivative of energy
    double ddW[3][3], // (out) 2nd derivative of energy
    //
    const double c[3], // (in) deformed triangle vertex positions
    double stiff_contact,
    double contact_clearance,
    const CInput_Contact& input )
{
  double n[3];
  double pd = input.penetrationNormal(n[0],n[1],n[2], c[0],c[1],c[2]);
  pd += contact_clearance;
  if( pd  < 0 ){
    W = 0;
    dW[0] = 0;  dW[1] = 0;  dW[2] = 0;
    ddW[0][0] = 0;  ddW[0][1] = 0;  ddW[0][2] = 0;
    ddW[1][0] = 0;  ddW[1][1] = 0;  ddW[1][2] = 0;
    ddW[2][0] = 0;  ddW[2][1] = 0;  ddW[2][2] = 0;
    return;
  }
  W = 0.5*stiff_contact*pd*pd;
  
  dW[0] = -stiff_contact*pd*n[0];
  dW[1] = -stiff_contact*pd*n[1];
  dW[2] = -stiff_contact*pd*n[2];
  
  ddW[0][0] = stiff_contact*n[0]*n[0];
  ddW[0][1] = stiff_contact*n[0]*n[1];
  ddW[0][2] = stiff_contact*n[0]*n[2];
  ddW[1][0] = stiff_contact*n[1]*n[0];
  ddW[1][1] = stiff_contact*n[1]*n[1];
  ddW[1][2] = stiff_contact*n[1]*n[2];
  ddW[2][0] = stiff_contact*n[2]*n[0];
  ddW[2][1] = stiff_contact*n[2]*n[1];
  ddW[2][2] = stiff_contact*n[2]*n[2];
}


