//
// Created by Nobuyuki Umetani on 2021/11/20.
//

#include "delfem2/fem_stvk.h"

#include "delfem2/vec3.h"

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
  Gu0 /= Dot(Gu0, Gd0);
  CVec3d Gu1 = Cross(Gd2, Gd0);
  Gu1 /= Dot(Gu1, Gd1);

  const CVec3d gd0(p[1][0] - p[0][0], p[1][1] - p[0][1], p[1][2] - p[0][2]);
  const CVec3d gd1(p[2][0] - p[0][0], p[2][1] - p[0][1], p[2][2] - p[0][2]);

  const double E2[3] = {  // green lagrange strain (with engineer's notation)
    0.5 * (Dot(gd0, gd0) - Dot(Gd0, Gd0)),
    0.5 * (Dot(gd1, gd1) - Dot(Gd1, Gd1)),
    1.0 * (Dot(gd0, gd1) - Dot(Gd0, Gd1))};

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
  Gu0 /= Dot(Gu0, Gd0);
  CVec3d Gu1 = Cross(Gd2, Gd0);
  Gu1 /= Dot(Gu1, Gd1);

  const CVec3d gd0(p[1][0] - p[0][0], p[1][1] - p[0][1], p[1][2] - p[0][2]);
  const CVec3d gd1(p[2][0] - p[0][0], p[2][1] - p[0][1], p[2][2] - p[0][2]);

  const double E2[3] = {  // green lagrange strain (with engineer's notation)
    0.5 * (Dot(gd0, gd0) - Dot(Gd0, Gd0)),
    0.5 * (Dot(gd1, gd1) - Dot(Gd1, Gd1)),
    1.0 * (Dot(gd0, gd1) - Dot(Gd0, Gd1))};

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
