//
// Created by Nobuyuki Umetani on 2021-11-20.
//

#include "delfem2/fem_quadratic_bending.h"

#include "delfem2/vec3.h"

DFM2_INLINE void delfem2::CdC_QuadBend(
    double C[3],
    double dCdp[3][4][3],
    const double P[4][3],
    const double p[4][3]) {
  const double A0 = Area_Tri3(P[0], P[2], P[3]);
  const double A1 = Area_Tri3(P[1], P[3], P[2]);
  const double L = Distance3(P[2], P[3]);
  const double H0 = 2.0 * A0 / L;
  const double H1 = 2.0 * A1 / L;
  const CVec3d e23(P[3][0] - P[2][0], P[3][1] - P[2][1], P[3][2] - P[2][2]);
  const CVec3d e02(P[2][0] - P[0][0], P[2][1] - P[0][1], P[2][2] - P[0][2]);
  const CVec3d e03(P[3][0] - P[0][0], P[3][1] - P[0][1], P[3][2] - P[0][2]);
  const CVec3d e12(P[2][0] - P[1][0], P[2][1] - P[1][1], P[2][2] - P[1][2]);
  const CVec3d e13(P[3][0] - P[1][0], P[3][1] - P[1][1], P[3][2] - P[1][2]);
  const double cot023 = -(e02.dot(e23)) / H0;
  const double cot032 = +(e03.dot(e23)) / H0;
  const double cot123 = -(e12.dot(e23)) / H1;
  const double cot132 = +(e13.dot(e23)) / H1;
  const double tmp0 = sqrt(3.0) / (sqrt(A0 + A1) * L);
  const double K[4] = {
      (-cot023 - cot032) * tmp0,
      (-cot123 - cot132) * tmp0,
      (+cot032 + cot132) * tmp0,
      (+cot023 + cot123) * tmp0};
  C[0] = K[0] * p[0][0] + K[1] * p[1][0] + K[2] * p[2][0] + K[3] * p[3][0];
  C[1] = K[0] * p[0][1] + K[1] * p[1][1] + K[2] * p[2][1] + K[3] * p[3][1];
  C[2] = K[0] * p[0][2] + K[1] * p[1][2] + K[2] * p[2][2] + K[3] * p[3][2];
  std::fill_n((&dCdp[0][0][0]), 36, 0.0);
  for (int idim = 0; idim < 3; ++idim) {
    dCdp[idim][0][idim] = K[0];
    dCdp[idim][1][idim] = K[1];
    dCdp[idim][2][idim] = K[2];
    dCdp[idim][3][idim] = K[3];
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
  const double A0 = Area_Tri3(C[0],C[2],C[3]);
  const double A1 = Area_Tri3(C[1],C[3],C[2]);
  const double L0 = Distance3(C[2],C[3]);
  const double H0 = A0*2.0/L0;
  const double H1 = A1*2.0/L0;
  const double e23[3] = { C[3][0]-C[2][0], C[3][1]-C[2][1], C[3][2]-C[2][2] };
  const double e02[3] = { C[2][0]-C[0][0], C[2][1]-C[0][1], C[2][2]-C[0][2] };
  const double e03[3] = { C[3][0]-C[0][0], C[3][1]-C[0][1], C[3][2]-C[0][2] };
  const double e12[3] = { C[2][0]-C[1][0], C[2][1]-C[1][1], C[2][2]-C[1][2] };
  const double e13[3] = { C[3][0]-C[1][0], C[3][1]-C[1][1], C[3][2]-C[1][2] };
  double cot023, cot032;
  {
    const double r2 = -Dot3(e02,e23);
    const double r3 = +Dot3(e03,e23);
    cot023 = r2/H0;
    cot032 = r3/H0;
  }
  double cot123, cot132;
  {
    const double r2 = -Dot3(e12,e23);
    const double r3 = +Dot3(e13,e23);
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