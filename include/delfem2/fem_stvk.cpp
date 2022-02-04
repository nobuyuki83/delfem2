//
// Created by Nobuyuki Umetani on 2021/11/20.
//

#include "delfem2/fem_stvk.h"

#include "delfem2/vec3.h"
#include "delfem2/geo_vec3.h"

/**
 *
 * @param[out] C
 * @param[out] dCdp
 * @param[in] P undeformed triangle vertex positions
 * @param[in] p deformed triangle vertex positions
 */
DFM2_INLINE void delfem2::CdC_StVK(
  double C[3],
  double dCdp[3][9],
  const double P[3][2],
  const double p[3][3]) {
  const CVec3d Gd0(P[1][0] - P[0][0], P[1][1] - P[0][1], 0.0);
  const CVec3d Gd1(P[2][0] - P[0][0], P[2][1] - P[0][1], 0.0);
  CVec3d Gd2 = Cross(Gd0, Gd1);
  const double Area = Gd2.norm() * 0.5;
  Gd2 /= (Area * 2.0);

  CVec3d Gu0 = Cross(Gd1, Gd2);
  Gu0 /= Gu0.dot(Gd0);
  CVec3d Gu1 = Cross(Gd2, Gd0);
  Gu1 /= Gu1.dot(Gd1);

  const CVec3d gd0(p[1][0] - p[0][0], p[1][1] - p[0][1], p[1][2] - p[0][2]);
  const CVec3d gd1(p[2][0] - p[0][0], p[2][1] - p[0][1], p[2][2] - p[0][2]);

  const double E2[3] = {  // green lagrange strain (with engineer's notation)
    0.5 * (gd0.dot(gd0) - Gd0.dot(Gd0)),
    0.5 * (gd1.dot(gd1) - Gd1.dot(Gd1)),
    1.0 * (gd0.dot(gd1) - Gd0.dot(Gd1))};

  const double GuGu_xx[3] = {Gu0.x * Gu0.x, Gu1.x * Gu1.x, Gu0.x * Gu1.x};
  const double GuGu_yy[3] = {Gu0.y * Gu0.y, Gu1.y * Gu1.y, Gu0.y * Gu1.y};
  const double GuGu_xy[3] = {2.0 * Gu0.x * Gu0.y, 2.0 * Gu1.x * Gu1.y, Gu0.x * Gu1.y + Gu0.y * Gu1.x};
  C[0] = E2[0] * GuGu_xx[0] + E2[1] * GuGu_xx[1] + E2[2] * GuGu_xx[2];
  C[1] = E2[0] * GuGu_yy[0] + E2[1] * GuGu_yy[1] + E2[2] * GuGu_yy[2];
  C[2] = E2[0] * GuGu_xy[0] + E2[1] * GuGu_xy[1] + E2[2] * GuGu_xy[2];

  const CVec3d dC0dp0 = -(GuGu_xx[0] + GuGu_xx[2]) * gd0 - (GuGu_xx[1] + GuGu_xx[2]) * gd1;
  const CVec3d dC0dp1 = GuGu_xx[0] * gd0 + GuGu_xx[2] * gd1;
  const CVec3d dC0dp2 = GuGu_xx[1] * gd1 + GuGu_xx[2] * gd0;
  const CVec3d dC1dp0 = -(GuGu_yy[0] + GuGu_yy[2]) * gd0 - (GuGu_yy[1] + GuGu_yy[2]) * gd1;
  const CVec3d dC1dp1 = GuGu_yy[0] * gd0 + GuGu_yy[2] * gd1;
  const CVec3d dC1dp2 = GuGu_yy[1] * gd1 + GuGu_yy[2] * gd0;
  const CVec3d dC2dp0 = -(GuGu_xy[0] + GuGu_xy[2]) * gd0 - (GuGu_xy[1] + GuGu_xy[2]) * gd1;
  const CVec3d dC2dp1 = GuGu_xy[0] * gd0 + GuGu_xy[2] * gd1;
  const CVec3d dC2dp2 = GuGu_xy[1] * gd1 + GuGu_xy[2] * gd0;

  dC0dp0.CopyTo(dCdp[0] + 0 * 3);
  dC0dp1.CopyTo(dCdp[0] + 1 * 3);
  dC0dp2.CopyTo(dCdp[0] + 2 * 3);
  dC1dp0.CopyTo(dCdp[1] + 0 * 3);
  dC1dp1.CopyTo(dCdp[1] + 1 * 3);
  dC1dp2.CopyTo(dCdp[1] + 2 * 3);
  dC2dp0.CopyTo(dCdp[2] + 0 * 3);
  dC2dp1.CopyTo(dCdp[2] + 1 * 3);
  dC2dp2.CopyTo(dCdp[2] + 2 * 3);
}

DFM2_INLINE void delfem2::CdC_EnergyStVK(
  double &C,
  double dCdp[9],
  const double P[3][2],
  const double p[3][3],
  const double lambda,
  const double myu) {
  const CVec3d Gd0(P[1][0] - P[0][0], P[1][1] - P[0][1], 0.0);
  const CVec3d Gd1(P[2][0] - P[0][0], P[2][1] - P[0][1], 0.0);
  CVec3d Gd2 = Cross(Gd0, Gd1);
  const double Area = Gd2.norm() * 0.5;
  Gd2 /= (Area * 2.0);

  CVec3d Gu0 = Cross(Gd1, Gd2);
  Gu0 /= Gu0.dot(Gd0);
  CVec3d Gu1 = Cross(Gd2, Gd0);
  Gu1 /= Gu1.dot(Gd1);

  const CVec3d gd0(p[1][0] - p[0][0], p[1][1] - p[0][1], p[1][2] - p[0][2]);
  const CVec3d gd1(p[2][0] - p[0][0], p[2][1] - p[0][1], p[2][2] - p[0][2]);

  const double E2[3] = {  // green lagrange strain (with engineer's notation)
    0.5 * (gd0.dot(gd0) - Gd0.dot(Gd0)),
    0.5 * (gd1.dot(gd1) - Gd1.dot(Gd1)),
    1.0 * (gd0.dot(gd1) - Gd0.dot(Gd1))};

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
  const double GuGu2[3] = {Gu0.dot(Gu0), Gu1.dot(Gu1), Gu1.dot(Gu0)};
  const double Cons2[3][3] = { // constitutive tensor
    {lambda * GuGu2[0] * GuGu2[0] + 2 * myu * (GuGu2[0] * GuGu2[0]),
     lambda * GuGu2[0] * GuGu2[1] + 2 * myu * (GuGu2[2] * GuGu2[2]),
     lambda * GuGu2[0] * GuGu2[2] + 2 * myu * (GuGu2[0] * GuGu2[2])},
    {lambda * GuGu2[1] * GuGu2[0] + 2 * myu * (GuGu2[2] * GuGu2[2]),
     lambda * GuGu2[1] * GuGu2[1] + 2 * myu * (GuGu2[1] * GuGu2[1]),
     lambda * GuGu2[1] * GuGu2[2] + 2 * myu * (GuGu2[2] * GuGu2[1])},
    {lambda * GuGu2[2] * GuGu2[0] + 2 * myu * (GuGu2[0] * GuGu2[2]),
     lambda * GuGu2[2] * GuGu2[1] + 2 * myu * (GuGu2[2] * GuGu2[1]),
     lambda * GuGu2[2] * GuGu2[2] + 1 * myu * (GuGu2[0] * GuGu2[1] + GuGu2[2] * GuGu2[2])}};
  const double S2[3] = {  // 2nd Piola-Kirchhoff stress
    Cons2[0][0] * E2[0] + Cons2[0][1] * E2[1] + Cons2[0][2] * E2[2],
    Cons2[1][0] * E2[0] + Cons2[1][1] * E2[1] + Cons2[1][2] * E2[2],
    Cons2[2][0] * E2[0] + Cons2[2][1] * E2[1] + Cons2[2][2] * E2[2]};

  // compute energy
  C = 0.5 * Area * (E2[0] * S2[0] + E2[1] * S2[1] + E2[2] * S2[2]);

  // compute 1st derivative
  const double dNdr[3][2] = {{-1.0, -1.0}, {+1.0, +0.0}, {+0.0, +1.0}};
  for (int ino = 0; ino < 3; ino++) {
    dCdp[ino * 3 + 0] = Area
      * (+S2[0] * gd0[0] * dNdr[ino][0] + S2[2] * gd0[0] * dNdr[ino][1] + S2[2] * gd1[0] * dNdr[ino][0]
        + S2[1] * gd1[0] * dNdr[ino][1]);
    dCdp[ino * 3 + 1] = Area
      * (+S2[0] * gd0[1] * dNdr[ino][0] + S2[2] * gd0[1] * dNdr[ino][1] + S2[2] * gd1[1] * dNdr[ino][0]
        + S2[1] * gd1[1] * dNdr[ino][1]);
    dCdp[ino * 3 + 2] = Area
      * (+S2[0] * gd0[2] * dNdr[ino][0] + S2[2] * gd0[2] * dNdr[ino][1] + S2[2] * gd1[2] * dNdr[ino][0]
        + S2[1] * gd1[2] * dNdr[ino][1]);
  }
}



DFM2_INLINE void MakeConstMatrix3D(
    double C[6][6],
    double lambda,
    double myu,
    const double Gu[3][3])
{
  const double GuGu2[6] = {
    delfem2::Dot3(Gu[0],Gu[0]), // 0 xx
    delfem2::Dot3(Gu[1],Gu[1]), // 1 yy
    delfem2::Dot3(Gu[2],Gu[2]), // 2 zz
    delfem2::Dot3(Gu[0],Gu[1]), // 3 xy
    delfem2::Dot3(Gu[1],Gu[2]), // 4 yz
    delfem2::Dot3(Gu[2],Gu[0])  // 5 zx
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
  UnitNormalAreaTri3(Gd[2], Area, C[0], C[1], C[2]);
  
  double Gu[2][3]; // inverse of Gd
  {
    Cross(Gu[0], Gd[1], Gd[2]);
    const double invtmp1 = 1.0/Dot3(Gu[0],Gd[0]);
    Gu[0][0] *= invtmp1;  Gu[0][1] *= invtmp1;  Gu[0][2] *= invtmp1;
    //
    Cross(Gu[1], Gd[2], Gd[0]);
    const double invtmp2 = 1.0/Dot3(Gu[1],Gd[1]);
    Gu[1][0] *= invtmp2;  Gu[1][1] *= invtmp2;  Gu[1][2] *= invtmp2;
  }
  
  const double gd[2][3] = { // deformed edge vector
    { c[1][0]-c[0][0], c[1][1]-c[0][1], c[1][2]-c[0][2] },
    { c[2][0]-c[0][0], c[2][1]-c[0][1], c[2][2]-c[0][2] } };
  
  const double E2[3] = {  // green lagrange strain (with engineer's notation)
    0.5*( Dot3(gd[0],gd[0]) - Dot3(Gd[0],Gd[0]) ),
    0.5*( Dot3(gd[1],gd[1]) - Dot3(Gd[1],Gd[1]) ),
    1.0*( Dot3(gd[0],gd[1]) - Dot3(Gd[0],Gd[1]) ) };
  const double GuGu2[3] = {
      Dot3(Gu[0],Gu[0]),
      Dot3(Gu[1],Gu[1]),
      Dot3(Gu[1],Gu[0]) };
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

void delfem2::WdWddW_CST_Sensitivity(
  double Kmat[3][3][3][3], 
  double Res[3][3], 
  double dRdC[3][3][3][3],
  double lambda, 
  double myu,
  const double C[3][3], 
  const double c[3][3])
{
  double Gd[3][3] = {
    { C[1][0]-C[0][0], C[1][1]-C[0][1], C[1][2]-C[0][2] },
    { C[2][0]-C[0][0], C[2][1]-C[0][1], C[2][2]-C[0][2] }, { 0,0,0 } };
  double Area;
  UnitNormalAreaTri3(Gd[2], Area, C[0], C[1], C[2]);
  
  double Gu[2][3];
  {
    double invArea = 0.5/Area;
    Cross(Gu[0], Gd[1], Gd[2]);
    Gu[0][0] *= invArea;  Gu[0][1] *= invArea;  Gu[0][2] *= invArea;
    //  
    Cross(Gu[1], Gd[2], Gd[0]);
    Gu[1][0] *= invArea ; Gu[1][1] *= invArea;  Gu[1][2] *= invArea;
  }
  
  const double gd[2][3] = { 
    { c[1][0]-c[0][0], c[1][1]-c[0][1], c[1][2]-c[0][2] },
    { c[2][0]-c[0][0], c[2][1]-c[0][1], c[2][2]-c[0][2] } };
  
  const double E2[3] = {  // engineering strain
    0.5*( Dot3(gd[0],gd[0]) - Dot3(Gd[0],Gd[0]) ),
    0.5*( Dot3(gd[1],gd[1]) - Dot3(Gd[1],Gd[1]) ),
    1.0*( Dot3(gd[0],gd[1]) - Dot3(Gd[0],Gd[1]) ) };
  
  const double E2dC[3][3][3] = { // E2dC[a][i][k] : derivative of E2[k] w.r.t. C[a][i]
    { {+Gd[0][0],+Gd[1][0],+Gd[0][0]+Gd[1][0]},
      {+Gd[0][1],+Gd[1][1],+Gd[0][1]+Gd[1][1]},
      {+Gd[0][2],+Gd[1][2],+Gd[0][2]+Gd[1][2]} },
    { {-Gd[0][0],0,        -Gd[1][0]},
      {-Gd[0][1],0,        -Gd[1][1]},
      {-Gd[0][2],0,        -Gd[1][2]} },
    { {0,        -Gd[1][0],-Gd[0][0]},
      {0,        -Gd[1][1],-Gd[0][1]},
      {0,        -Gd[1][2],-Gd[0][2]} } };
  
  double DAdC[3][3];
  {
    const double D[3][3] = {
      { C[2][0]-C[1][0], C[2][1]-C[1][1], C[2][2]-C[1][2] },
      { C[0][0]-C[2][0], C[0][1]-C[2][1], C[0][2]-C[2][2] },
      { C[1][0]-C[0][0], C[1][1]-C[0][1], C[1][2]-C[0][2] } };
    Cross(DAdC[0], Gd[2], D[0]);
    Cross(DAdC[1], Gd[2], D[1]);
    Cross(DAdC[2], Gd[2], D[2]);
  }
  
  double GuDC[3][3][2][3];  // GuDC[a][i][b][j] : derivative of Gu[b][j] w.r.t. C[a][i]
  {
    const double invSqDArea = 0.25/(Area*Area);
    double G101[3]; Cross(G101,  Gd[1], Gd[2]);
    double G010[3]; Cross(G010,  Gd[2], Gd[0]);
    for(unsigned int idim=0;idim<3;idim++){
      for(unsigned int jdim=0;jdim<3;jdim++){      
        GuDC[0][idim][0][jdim] = invSqDArea*( -2.0*G101[jdim]*DAdC[0][idim] - 2.0*Gd[1][idim]*Gd[0][jdim] + (Gd[0][idim]+Gd[1][idim])*Gd[1][jdim] );
        GuDC[1][idim][0][jdim] = invSqDArea*( -2.0*G101[jdim]*DAdC[1][idim] - Gd[1][idim]*Gd[1][jdim] );
        GuDC[2][idim][0][jdim] = invSqDArea*( -2.0*G101[jdim]*DAdC[2][idim] + 2.0*Gd[1][idim]*Gd[0][jdim] - Gd[0][idim]*Gd[1][jdim] );
        GuDC[0][idim][1][jdim] = invSqDArea*( -2.0*G010[jdim]*DAdC[0][idim] - 2.0*Gd[0][idim]*Gd[1][jdim] + (Gd[0][idim]+Gd[1][idim])*Gd[0][jdim] );
        GuDC[1][idim][1][jdim] = invSqDArea*( -2.0*G010[jdim]*DAdC[1][idim] + 2.0*Gd[0][idim]*Gd[1][jdim] - Gd[1][idim]*Gd[0][jdim] );
        GuDC[2][idim][1][jdim] = invSqDArea*( -2.0*G010[jdim]*DAdC[2][idim] - Gd[0][idim]*Gd[0][jdim] );
      }
      GuDC[0][idim][0][idim] += invSqDArea*( -Dot3(Gd[1],Gd[1]) + Dot3(Gd[0],Gd[1]) );
      GuDC[1][idim][0][idim] += invSqDArea*( +Dot3(Gd[1],Gd[1]) );
      GuDC[2][idim][0][idim] += invSqDArea*( -Dot3(Gd[1],Gd[0]) );
      GuDC[0][idim][1][idim] += invSqDArea*( -Dot3(Gd[0],Gd[0]) + Dot3(Gd[0],Gd[1]) );
      GuDC[1][idim][1][idim] += invSqDArea*( -Dot3(Gd[1],Gd[0]) );
      GuDC[2][idim][1][idim] += invSqDArea*( +Dot3(Gd[0],Gd[0]) );
    }
  }
  
  double GuGu2dC[3][3][3];  // GuDC[a][i][k] : derivative of GuGu[k] w.r.t. C[a][i]
  for(unsigned int jno=0;jno<3;jno++){
  for(unsigned int jdim=0;jdim<3;jdim++){
    GuGu2dC[jno][jdim][0] = 2.0*Dot3(GuDC[jno][jdim][0],Gu[0]);
    GuGu2dC[jno][jdim][1] = 2.0*Dot3(GuDC[jno][jdim][1],Gu[1]);
    GuGu2dC[jno][jdim][2] = Dot3(GuDC[jno][jdim][0],Gu[1]) + Dot3(GuDC[jno][jdim][1],Gu[0]);
  }
  }
  const double GuGu2[3] = {
      Dot3(Gu[0],Gu[0]),
      Dot3(Gu[1],Gu[1]),
      Dot3(Gu[1],Gu[0]) };
  
  const double Cons2[6] = {
    (lambda+2*myu)*GuGu2[0]*GuGu2[0],
    (lambda+2*myu)*GuGu2[1]*GuGu2[1],
    (lambda+  myu)*GuGu2[2]*GuGu2[2] + myu*GuGu2[0]*GuGu2[1],
    (lambda+2*myu)*GuGu2[1]*GuGu2[2],
    (lambda+2*myu)*GuGu2[0]*GuGu2[2],
    lambda*GuGu2[0]*GuGu2[1] + 2*myu*(GuGu2[2]*GuGu2[2]) };   
  
  double Cons2dC[3][3][6];
  {
    for(unsigned int jno=0;jno<3;jno++){
    for(unsigned int jdim=0;jdim<3;jdim++){
      Cons2dC[jno][jdim][0] = 2*(lambda+2*myu)*GuGu2[0]*GuGu2dC[jno][jdim][0];
      Cons2dC[jno][jdim][1] = 2*(lambda+2*myu)*GuGu2[1]*GuGu2dC[jno][jdim][1];
      Cons2dC[jno][jdim][2] = 2*(lambda+  myu)*GuGu2[2]*GuGu2dC[jno][jdim][2] + myu*( GuGu2[1]*GuGu2dC[jno][jdim][0] + GuGu2[0]*GuGu2dC[jno][jdim][1] );
      Cons2dC[jno][jdim][3] =   (lambda+2*myu)*( GuGu2[1]*GuGu2dC[jno][jdim][2] + GuGu2[2]*GuGu2dC[jno][jdim][1] );
      Cons2dC[jno][jdim][4] =   (lambda+2*myu)*( GuGu2[0]*GuGu2dC[jno][jdim][2] + GuGu2[2]*GuGu2dC[jno][jdim][0] );
      Cons2dC[jno][jdim][5] = lambda*( GuGu2[1]*GuGu2dC[jno][jdim][0] + GuGu2[0]*GuGu2dC[jno][jdim][1] ) + 4*myu*GuGu2[2]*GuGu2dC[jno][jdim][2];
    }
    }    
  }

  const double S2[3] = {  // stress
    Cons2[0]*E2[0] + Cons2[5]*E2[1] + Cons2[4]*E2[2], 
    Cons2[5]*E2[0] + Cons2[1]*E2[1] + Cons2[3]*E2[2], 
    Cons2[4]*E2[0] + Cons2[3]*E2[1] + Cons2[2]*E2[2] };  
  
  const double dNdr[3][2] = { {-1.0, -1.0}, {+1.0, +0.0}, {+0.0, +1.0} };    
  for(unsigned int ino=0;ino<3;ino++){
    for(unsigned int idim=0;idim<3;idim++){
      Res[ino][idim] = Area*
      (+S2[0]*gd[0][idim]*dNdr[ino][0]
       +S2[2]*gd[0][idim]*dNdr[ino][1]
       +S2[2]*gd[1][idim]*dNdr[ino][0]
       +S2[1]*gd[1][idim]*dNdr[ino][1]);
    }
  }  
  
  double S2dC[3][3][3];
  for(unsigned int jno=0;jno<3;jno++){
  for(unsigned int jdim=0;jdim<3;jdim++){
    S2dC[jno][jdim][0]
    = Cons2[0]*E2dC[jno][jdim][0] + Cons2dC[jno][jdim][0]*E2[0]
    + Cons2[5]*E2dC[jno][jdim][1] + Cons2dC[jno][jdim][5]*E2[1]
    + Cons2[4]*E2dC[jno][jdim][2] + Cons2dC[jno][jdim][4]*E2[2];
    S2dC[jno][jdim][1]
    = Cons2[5]*E2dC[jno][jdim][0] + Cons2dC[jno][jdim][5]*E2[0]
    + Cons2[1]*E2dC[jno][jdim][1] + Cons2dC[jno][jdim][1]*E2[1]
    + Cons2[3]*E2dC[jno][jdim][2] + Cons2dC[jno][jdim][3]*E2[2];    
    S2dC[jno][jdim][2]
    = Cons2[4]*E2dC[jno][jdim][0] + Cons2dC[jno][jdim][4]*E2[0]
    + Cons2[3]*E2dC[jno][jdim][1] + Cons2dC[jno][jdim][3]*E2[1]
    + Cons2[2]*E2dC[jno][jdim][2] + Cons2dC[jno][jdim][2]*E2[2];        
  }    
  }

  for(unsigned int ino=0;ino<3;ino++){
  for(unsigned int idim=0;idim<3;idim++){
    for(unsigned int jno=0;jno<3;jno++){
    for(unsigned int jdim=0;jdim<3;jdim++){
      dRdC[ino][jno][idim][jdim] = 
      +(Area*S2dC[jno][jdim][0]+0.5*DAdC[jno][jdim]*S2[0])*gd[0][idim]*dNdr[ino][0]
      +(Area*S2dC[jno][jdim][2]+0.5*DAdC[jno][jdim]*S2[2])*gd[0][idim]*dNdr[ino][1]
      +(Area*S2dC[jno][jdim][2]+0.5*DAdC[jno][jdim]*S2[2])*gd[1][idim]*dNdr[ino][0]
      +(Area*S2dC[jno][jdim][1]+0.5*DAdC[jno][jdim]*S2[1])*gd[1][idim]*dNdr[ino][1];
    }
    }
  }
  }    
  
  // ---------------------
  // Calc stiffness matrix 
  
  double S3[3] = { S2[0], S2[1], S2[2] };
  MakePositiveDefinite_Sim22(S2,S3);
  for(unsigned int ino=0;ino<3;ino++){    
    for(unsigned int jno=0;jno<ino+1;jno++){
      for(unsigned int idim=0;idim<3;idim++){
        for(unsigned int jdim=0;jdim<3;jdim++){ 
          double dtmp0 = 0;         
          dtmp0 += gd[0][idim]*(dNdr[ino][0]*Cons2[0]+dNdr[ino][1]*Cons2[4])*gd[0][jdim]*dNdr[jno][0];
          dtmp0 += gd[0][idim]*(dNdr[ino][0]*Cons2[5]+dNdr[ino][1]*Cons2[3])*gd[1][jdim]*dNdr[jno][1];
          dtmp0 += gd[0][idim]*(dNdr[ino][0]*Cons2[4]+dNdr[ino][1]*Cons2[2])*gd[0][jdim]*dNdr[jno][1];
          dtmp0 += gd[0][idim]*(dNdr[ino][0]*Cons2[4]+dNdr[ino][1]*Cons2[2])*gd[1][jdim]*dNdr[jno][0];
          dtmp0 += gd[1][idim]*(dNdr[ino][0]*Cons2[4]+dNdr[ino][1]*Cons2[5])*gd[0][jdim]*dNdr[jno][0];
          dtmp0 += gd[1][idim]*(dNdr[ino][0]*Cons2[3]+dNdr[ino][1]*Cons2[1])*gd[1][jdim]*dNdr[jno][1];
          dtmp0 += gd[1][idim]*(dNdr[ino][0]*Cons2[2]+dNdr[ino][1]*Cons2[3])*gd[0][jdim]*dNdr[jno][1];
          dtmp0 += gd[1][idim]*(dNdr[ino][0]*Cons2[2]+dNdr[ino][1]*Cons2[3])*gd[1][jdim]*dNdr[jno][0];
          Kmat[ino][jno][idim][jdim] = dtmp0*Area;
        }
      }
      const double dtmp1 = Area*
      (+S3[0]*dNdr[ino][0]*dNdr[jno][0]
       +S3[2]*dNdr[ino][0]*dNdr[jno][1]
       +S3[2]*dNdr[ino][1]*dNdr[jno][0]
       +S3[1]*dNdr[ino][1]*dNdr[jno][1]);
      Kmat[ino][jno][0][0] += dtmp1;
      Kmat[ino][jno][1][1] += dtmp1;
      Kmat[ino][jno][2][2] += dtmp1;
    }
  }
  
  for(unsigned int ino=0;ino<3;ino++){    
    for(unsigned int jno=ino+1;jno<3;jno++){
      for(unsigned int idim=0;idim<3;idim++){
        for(unsigned int jdim=0;jdim<3;jdim++){
          Kmat[ino][jno][idim][jdim] = Kmat[jno][ino][jdim][idim];
        }
      }
    }
  }                        
}
