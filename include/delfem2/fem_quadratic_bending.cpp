//
// Created by Nobuyuki Umetani on 2021-11-20.
//

#include "delfem2/fem_quadratic_bending.h"

#include "delfem2/vec3.h"
#include "delfem2/geo_vec3.h"

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
  double &W,  // (out) strain energy
  double dW[4][3], // (out) 1st derivative of energy
  double ddW[4][4][3][3], // (out) 2nd derivative of energy
  //
  const double C[4][3], // (in) undeformed triangle vertex positions
  const double c[4][3], // (in) deformed triangle vertex positions
  double stiff) {
  const double A0 = Area_Tri3(C[0], C[2], C[3]);
  const double A1 = Area_Tri3(C[1], C[3], C[2]);
  const double L0 = Distance3(C[2], C[3]);
  const double H0 = A0 * 2.0 / L0;
  const double H1 = A1 * 2.0 / L0;
  const double e23[3] = {C[3][0] - C[2][0], C[3][1] - C[2][1], C[3][2] - C[2][2]};
  const double e02[3] = {C[2][0] - C[0][0], C[2][1] - C[0][1], C[2][2] - C[0][2]};
  const double e03[3] = {C[3][0] - C[0][0], C[3][1] - C[0][1], C[3][2] - C[0][2]};
  const double e12[3] = {C[2][0] - C[1][0], C[2][1] - C[1][1], C[2][2] - C[1][2]};
  const double e13[3] = {C[3][0] - C[1][0], C[3][1] - C[1][1], C[3][2] - C[1][2]};
  double cot023, cot032;
  {
    const double r2 = -Dot3(e02, e23);
    const double r3 = +Dot3(e03, e23);
    cot023 = r2 / H0;
    cot032 = r3 / H0;
  }
  double cot123, cot132;
  {
    const double r2 = -Dot3(e12, e23);
    const double r3 = +Dot3(e13, e23);
    cot123 = r2 / H1;
    cot132 = r3 / H1;
  }
  const double tmp0 = stiff / ((A0 + A1) * L0 * L0);
  const double K[4] = {-cot023 - cot032, -cot123 - cot132, cot032 + cot132, cot023 + cot123};

  // compute 2nd derivative of energy
  for (int i = 0; i < 4 * 4 * 3 * 3; i++) { (&ddW[0][0][0][0])[i] = 0; }
  for (int ino = 0; ino < 4; ino++) {
    for (int jno = 0; jno < 4; jno++) {
      const double tmp = K[ino] * K[jno] * tmp0;
      ddW[ino][jno][0][0] = tmp;
      ddW[ino][jno][1][1] = tmp;
      ddW[ino][jno][2][2] = tmp;
    }
  }
  // compute 1st derivative of energy
  W = 0.0;
  for (int ino = 0; ino < 4; ino++) {
    for (int idim = 0; idim < 3; idim++) {
      dW[ino][idim] = 0;
      for (int jno = 0; jno < 4; jno++) {
        for (int jdim = 0; jdim < 3; jdim++) {
          dW[ino][idim] += ddW[ino][jno][idim][jdim] * c[jno][jdim];
        }
      }
      W += dW[ino][idim] * c[ino][idim];
    }
  }
}

namespace delfem2::fem_quadratic_bending {
void DerDoubleAreaTri3D(
  double dAdC[3][3],
  const double c0[3],
  const double c1[3],
  const double c2[3],
  const double un[3])    // unit normal
{
  const double c[3][3] = {
    {c2[0] - c1[0], c2[1] - c1[1], c2[2] - c1[2]},
    {c0[0] - c2[0], c0[1] - c2[1], c0[2] - c2[2]},
    {c1[0] - c0[0], c1[1] - c0[1], c1[2] - c0[2]}};
  delfem2::Cross3(dAdC[0], un, c[0]);
  delfem2::Cross3(dAdC[1], un, c[1]);
  delfem2::Cross3(dAdC[2], un, c[2]);
}

}  // namespace delfem2::fem_quadratic_bending



void delfem2::WdWddW_QuadraticBending_Sensitivity(
  double ddW[4][4][3][3],
  double dW[4][3],
  double dRdC[4][4][3][3],
  const double C[4][3],
  const double c[4][3],
  double stiff) {
  namespace dfm2 = delfem2;
  namespace lcl = delfem2::fem_quadratic_bending;
  double A0, UN0[3];
  dfm2::UnitNormalAreaTri3(UN0, A0, C[0], C[2], C[3]);
  double A1, UN1[3];
  dfm2::UnitNormalAreaTri3(UN1, A1, C[1], C[3], C[2]);
  const double L0 = dfm2::Distance3(C[2], C[3]);
  double coeff = stiff / ((A0 + A1) * L0 * L0);
  double dcoeffdC[4][3];
  {
    double dA0dC023[3][3];
    lcl::DerDoubleAreaTri3D(dA0dC023, C[0], C[2], C[3], UN0);  // 023
    double dA1dC132[3][3];
    lcl::DerDoubleAreaTri3D(dA1dC132, C[1], C[3], C[2], UN1);  // 132
    const double t0 = stiff / ((A0 + A1) * (A0 + A1) * L0 * L0 * L0 * L0);
    const double u23[3] = {(C[3][0] - C[2][0]) / L0, (C[3][1] - C[2][1]) / L0, (C[3][2] - C[2][2]) / L0};
    for (unsigned int i = 0; i < 3; i++) {
      dcoeffdC[0][i] = -dA0dC023[0][i] * L0 * L0 * 0.5 * t0;
      dcoeffdC[1][i] = -dA1dC132[0][i] * L0 * L0 * 0.5 * t0;
      dcoeffdC[2][i] = -((dA0dC023[1][i] + dA1dC132[2][i]) * L0 * L0 * 0.5 - (A0 + A1) * 2.0 * L0 * u23[i]) * t0;
      dcoeffdC[3][i] = -((dA0dC023[2][i] + dA1dC132[1][i]) * L0 * L0 * 0.5 + (A0 + A1) * 2.0 * L0 * u23[i]) * t0;
    }
  }
  const double H0 = A0 * 2.0 / L0;
  const double H1 = A1 * 2.0 / L0;
  const double e23[3] = {C[3][0] - C[2][0], C[3][1] - C[2][1], C[3][2] - C[2][2]};
  const double e02[3] = {C[2][0] - C[0][0], C[2][1] - C[0][1], C[2][2] - C[0][2]};
  const double e03[3] = {C[3][0] - C[0][0], C[3][1] - C[0][1], C[3][2] - C[0][2]};
  const double e12[3] = {C[2][0] - C[1][0], C[2][1] - C[1][1], C[2][2] - C[1][2]};
  const double e13[3] = {C[3][0] - C[1][0], C[3][1] - C[1][1], C[3][2] - C[1][2]};
  double dKdC[4][3][4];
  double cot023, cot032, cot123, cot132;
  {
    const double d02 = -dfm2::Dot3(e02, e23);
    const double d03 = +dfm2::Dot3(e03, e23);
    cot023 = d02 / H0;
    cot032 = d03 / H0;
    const double r02 = d03 / (d02 + d03);
    const double invH0 = 1.0 / H0;
    const double uvH0[3] = {
      (r02 * C[2][0] + (1 - r02) * C[3][0] - C[0][0]) * invH0,
      (r02 * C[2][1] + (1 - r02) * C[3][1] - C[0][1]) * invH0,
      (r02 * C[2][2] + (1 - r02) * C[3][2] - C[0][2]) * invH0};
//    std::cout << uvH0[0]*uvH0[0] + uvH0[1]*uvH0[1] + uvH0[2]*uvH0[2] << std::endl;
    const double dcotdC023[3][3] = {
      {+e23[0] * invH0 + d02 * uvH0[0] * invH0 * invH0,
       +e23[1] * invH0 + d02 * uvH0[1] * invH0 * invH0,
       +e23[2] * invH0 + d02 * uvH0[2] * invH0 * invH0},
      {(e02[0] - e23[0]) * invH0 - r02 * d02 * uvH0[0] * invH0 * invH0,
       (e02[1] - e23[1]) * invH0 - r02 * d02 * uvH0[1] * invH0 * invH0,
       (e02[2] - e23[2]) * invH0 - r02 * d02 * uvH0[2] * invH0 * invH0},
      {-e02[0] * invH0 - (1 - r02) * d02 * uvH0[0] * invH0 * invH0,
       -e02[1] * invH0 - (1 - r02) * d02 * uvH0[1] * invH0 * invH0,
       -e02[2] * invH0 - (1 - r02) * d02 * uvH0[2] * invH0 * invH0}
    };
    const double dcotdC032[3][3] = {
      {-e23[0] * invH0 + d03 * uvH0[0] * invH0 * invH0,
       -e23[1] * invH0 + d03 * uvH0[1] * invH0 * invH0,
       -e23[2] * invH0 + d03 * uvH0[2] * invH0 * invH0},
      {(e23[0] + e03[0]) * invH0 - (1 - r02) * d03 * uvH0[0] * invH0 * invH0,
       (e23[1] + e03[1]) * invH0 - (1 - r02) * d03 * uvH0[1] * invH0 * invH0,
       (e23[2] + e03[2]) * invH0 - (1 - r02) * d03 * uvH0[2] * invH0 * invH0},
      {-e03[0] * invH0 - r02 * d03 * uvH0[0] * invH0 * invH0,
       -e03[1] * invH0 - r02 * d03 * uvH0[1] * invH0 * invH0,
       -e03[2] * invH0 - r02 * d03 * uvH0[2] * invH0 * invH0}
    };
    const double d12 = -dfm2::Dot3(e12, e23);
    const double d13 = +dfm2::Dot3(e13, e23);
    cot123 = d12 / H1;
    cot132 = d13 / H1;
    const double r12 = d13 / (d12 + d13);
    const double invH1 = 1.0 / H1;
    const double uvH1[3] = {
      (r12 * C[2][0] + (1 - r12) * C[3][0] - C[1][0]) * invH1,
      (r12 * C[2][1] + (1 - r12) * C[3][1] - C[1][1]) * invH1,
      (r12 * C[2][2] + (1 - r12) * C[3][2] - C[1][2]) * invH1};
//    std::cout << uvH0[0]*uvH0[0] + uvH0[1]*uvH0[1] + uvH0[2]*uvH0[2] << std::endl;
    const double dcotdC123[3][3] = {
      {+e23[0] * invH1 + d12 * uvH1[0] * invH1 * invH1,
       +e23[1] * invH1 + d12 * uvH1[1] * invH1 * invH1,
       +e23[2] * invH1 + d12 * uvH1[2] * invH1 * invH1},
      {(e12[0] - e23[0]) * invH1 - r12 * d12 * uvH1[0] * invH1 * invH1,
       (e12[1] - e23[1]) * invH1 - r12 * d12 * uvH1[1] * invH1 * invH1,
       (e12[2] - e23[2]) * invH1 - r12 * d12 * uvH1[2] * invH1 * invH1},
      {-e12[0] * invH1 - (1 - r12) * d12 * uvH1[0] * invH1 * invH1,
       -e12[1] * invH1 - (1 - r12) * d12 * uvH1[1] * invH1 * invH1,
       -e12[2] * invH1 - (1 - r12) * d12 * uvH1[2] * invH1 * invH1}
    };
    const double dcotdC132[3][3] = {
      {-e23[0] * invH1 + d13 * uvH1[0] * invH1 * invH1,
       -e23[1] * invH1 + d13 * uvH1[1] * invH1 * invH1,
       -e23[2] * invH1 + d13 * uvH1[2] * invH1 * invH1},
      {(e23[0] + e13[0]) * invH1 - (1 - r12) * d13 * uvH1[0] * invH1 * invH1,
       (e23[1] + e13[1]) * invH1 - (1 - r12) * d13 * uvH1[1] * invH1 * invH1,
       (e23[2] + e13[2]) * invH1 - (1 - r12) * d13 * uvH1[2] * invH1 * invH1},
      {-e13[0] * invH1 - r12 * d13 * uvH1[0] * invH1 * invH1,
       -e13[1] * invH1 - r12 * d13 * uvH1[1] * invH1 * invH1,
       -e13[2] * invH1 - r12 * d13 * uvH1[2] * invH1 * invH1}
    };
    for (unsigned int i = 0; i < 3; i++) {
      dKdC[0][i][0] = -dcotdC023[0][i] - dcotdC032[0][i];
      dKdC[1][i][0] = 0;
      dKdC[2][i][0] = -dcotdC023[1][i] - dcotdC032[2][i];
      dKdC[3][i][0] = -dcotdC023[2][i] - dcotdC032[1][i];
      //
      dKdC[0][i][1] = 0;
      dKdC[1][i][1] = -dcotdC123[0][i] - dcotdC132[0][i];
      dKdC[2][i][1] = -dcotdC123[1][i] - dcotdC132[2][i];
      dKdC[3][i][1] = -dcotdC123[2][i] - dcotdC132[1][i];
      //
      dKdC[0][i][2] = +dcotdC032[0][i];
      dKdC[1][i][2] = +dcotdC132[0][i];
      dKdC[2][i][2] = +dcotdC032[2][i] + dcotdC132[2][i];
      dKdC[3][i][2] = +dcotdC032[1][i] + dcotdC132[1][i];
      //
      dKdC[0][i][3] = +dcotdC023[0][i];
      dKdC[1][i][3] = +dcotdC123[0][i];
      dKdC[2][i][3] = +dcotdC023[1][i] + dcotdC123[1][i];
      dKdC[3][i][3] = +dcotdC023[2][i] + dcotdC123[2][i];
    }
  }
  const double K[4] = {-cot023 - cot032, -cot123 - cot132, cot032 + cot132, cot023 + cot123};
  std::fill_n(&ddW[0][0][0][0], 144, 0.0);
  for (unsigned int ino = 0; ino < 4; ino++) {
    for (unsigned int jno = 0; jno < 4; jno++) {
      const double tmp = K[ino] * K[jno] * coeff;
      ddW[ino][jno][0][0] = tmp;
      ddW[ino][jno][1][1] = tmp;
      ddW[ino][jno][2][2] = tmp;
    }
  }
  for (unsigned int ino = 0; ino < 4; ino++) {
    for (unsigned int idim = 0; idim < 3; idim++) {
      dW[ino][idim] = 0;
      for (unsigned int jno = 0; jno < 4; jno++) {
        for (unsigned int jdim = 0; jdim < 3; jdim++) {
          dW[ino][idim] += ddW[ino][jno][idim][jdim] * c[jno][jdim];
        }
      }
    }
  }

  const double Kc[3] = {
    K[0] * c[0][0] + K[1] * c[1][0] + K[2] * c[2][0] + K[3] * c[3][0],
    K[0] * c[0][1] + K[1] * c[1][1] + K[2] * c[2][1] + K[3] * c[3][1],
    K[0] * c[0][2] + K[1] * c[1][2] + K[2] * c[2][2] + K[3] * c[3][2]};
  for (unsigned int ino = 0; ino < 4; ino++) {
    for (unsigned int jno = 0; jno < 4; jno++) {
      for (unsigned int idim = 0; idim < 3; idim++) {
        for (unsigned int jdim = 0; jdim < 3; jdim++) {
          dRdC[ino][jno][idim][jdim] =
            dcoeffdC[jno][jdim] * K[ino] * Kc[idim]
              + coeff * dKdC[jno][jdim][ino] * Kc[idim]
              + coeff * K[ino] * (
                dKdC[jno][jdim][0] * c[0][idim] +
                  dKdC[jno][jdim][1] * c[1][idim] +
                  dKdC[jno][jdim][2] * c[2][idim] +
                  dKdC[jno][jdim][3] * c[3][idim]);

        }
      }
    }
  }
}
