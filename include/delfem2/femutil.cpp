/*
 * Copyright (c) 2019-2021 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */

#include "delfem2/femutil.h"

#include <complex>
#include <cassert>

// area of a triangle
DFM2_INLINE double delfem2::femutil::TriArea2D(
    const double p0[],
    const double p1[],
    const double p2[]) {
  return 0.5 * ((p1[0] - p0[0]) * (p2[1] - p0[1]) - (p2[0] - p0[0]) * (p1[1] - p0[1]));
}

DFM2_INLINE double delfem2::femutil::TriArea3D(
    const double v1[3],
    const double v2[3],
    const double v3[3]) {
  double x, y, z;
  x = (v2[1] - v1[1]) * (v3[2] - v1[2]) - (v3[1] - v1[1]) * (v2[2] - v1[2]);
  y = (v2[2] - v1[2]) * (v3[0] - v1[0]) - (v3[2] - v1[2]) * (v2[0] - v1[0]);
  z = (v2[0] - v1[0]) * (v3[1] - v1[1]) - (v3[0] - v1[0]) * (v2[1] - v1[1]);
  return 0.5 * sqrt(x * x + y * y + z * z);
}

DFM2_INLINE double delfem2::femutil::Distance3D(
    const double p0[3],
    const double p1[3]) {
  return sqrt(
      (p1[0] - p0[0]) * (p1[0] - p0[0]) + (p1[1] - p0[1]) * (p1[1] - p0[1]) + (p1[2] - p0[2]) * (p1[2] - p0[2]));
}

DFM2_INLINE void delfem2::femutil::UnitNormalAreaTri3D(
    double n[3],
    double &a,
    const double v1[3],
    const double v2[3],
    const double v3[3]) {
  n[0] = (v2[1] - v1[1]) * (v3[2] - v1[2]) - (v3[1] - v1[1]) * (v2[2] - v1[2]);
  n[1] = (v2[2] - v1[2]) * (v3[0] - v1[0]) - (v3[2] - v1[2]) * (v2[0] - v1[0]);
  n[2] = (v2[0] - v1[0]) * (v3[1] - v1[1]) - (v3[0] - v1[0]) * (v2[1] - v1[1]);
  a = sqrt(n[0] * n[0] + n[1] * n[1] + n[2] * n[2]) * 0.5;
  const double invlen = 0.5 / a;
  n[0] *= invlen;
  n[1] *= invlen;
  n[2] *= invlen;
}

// area coordinate inside a triangle
DFM2_INLINE void delfem2::femutil::TriAreaCoord(
    double vc_p[],
    const double p0[],
    const double p1[],
    const double p2[],
    const double pb[]) {
  vc_p[0] = TriArea2D(pb, p1, p2);
  vc_p[1] = TriArea2D(p0, pb, p2);
  vc_p[2] = TriArea2D(p0, p1, pb);

  const double area = TriArea2D(p0, p1, p2);
  const double inv_area = 1.0 / area;

  vc_p[0] = vc_p[0] * inv_area;
  vc_p[1] = vc_p[1] * inv_area;
  vc_p[2] = vc_p[2] * inv_area;

  assert(fabs(vc_p[0] + vc_p[1] + vc_p[2] - 1.0) < 1.0e-15);
}

DFM2_INLINE double delfem2::femutil::Dot3D(const double a[], const double b[]) {
  return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}

DFM2_INLINE void delfem2::femutil::Cross3D(double r[3], const double v1[3], const double v2[3]) {
  r[0] = v1[1] * v2[2] - v2[1] * v1[2];
  r[1] = v1[2] * v2[0] - v2[2] * v1[0];
  r[2] = v1[0] * v2[1] - v2[0] * v1[1];
}

DFM2_INLINE void delfem2::femutil::MatVec3(
    double y[3],
    const double m[9], const double x[3]) {
  y[0] = m[0] * x[0] + m[1] * x[1] + m[2] * x[2];
  y[1] = m[3] * x[0] + m[4] * x[1] + m[5] * x[2];
  y[2] = m[6] * x[0] + m[7] * x[1] + m[8] * x[2];
}

DFM2_INLINE void delfem2::femutil::MatMat3(
    double *C,
    const double *A,
    const double *B) {
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      C[i * 3 + j] = A[i * 3 + 0] * B[0 * 3 + j] + A[i * 3 + 1] * B[1 * 3 + j] + A[i * 3 + 2] * B[2 * 3 + j];
    }
  }
}

DFM2_INLINE void delfem2::femutil::MatMatTrans3(
    double *C,
    const double *A,
    const double *B) {
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      C[i * 3 + j] = A[i * 3 + 0] * B[j * 3 + 0] + A[i * 3 + 1] * B[j * 3 + 1] + A[i * 3 + 2] * B[j * 3 + 2];
    }
  }
}

namespace delfem2 {
namespace femutil {

void HexVox(
    double &detjac,
    double dndx[][3],
    [[maybe_unused]] double an[],
    const double coords[][3],
    double dndr[8][3]) {
  double dxdr[3][3] = {
      {0.0, 0.0, 0.0},
      {0.0, 0.0, 0.0},
      {0.0, 0.0, 0.0},
  };

  for (int inode = 0; inode < 8; inode++) {
    dxdr[0][0] += coords[inode][0] * dndr[inode][0];
    dxdr[0][1] += coords[inode][0] * dndr[inode][1];
    dxdr[0][2] += coords[inode][0] * dndr[inode][2];
    dxdr[1][0] += coords[inode][1] * dndr[inode][0];
    dxdr[1][1] += coords[inode][1] * dndr[inode][1];
    dxdr[1][2] += coords[inode][1] * dndr[inode][2];
    dxdr[2][0] += coords[inode][2] * dndr[inode][0];
    dxdr[2][1] += coords[inode][2] * dndr[inode][1];
    dxdr[2][2] += coords[inode][2] * dndr[inode][2];
  }

  detjac = +dxdr[0][0] * dxdr[1][1] * dxdr[2][2]
      + dxdr[1][0] * dxdr[2][1] * dxdr[0][2]
      + dxdr[2][0] * dxdr[0][1] * dxdr[1][2]
      - dxdr[0][0] * dxdr[2][1] * dxdr[1][2]
      - dxdr[1][0] * dxdr[0][1] * dxdr[2][2]
      - dxdr[2][0] * dxdr[1][1] * dxdr[0][2];

  const double inv_jac = 1.0 / detjac;

  double drdx[3][3];
  drdx[0][0] = inv_jac * (dxdr[1][1] * dxdr[2][2] - dxdr[1][2] * dxdr[2][1]);
  drdx[0][1] = inv_jac * (dxdr[0][2] * dxdr[2][1] - dxdr[0][1] * dxdr[2][2]);
  drdx[0][2] = inv_jac * (dxdr[0][1] * dxdr[1][2] - dxdr[0][2] * dxdr[1][1]);
  drdx[1][0] = inv_jac * (dxdr[1][2] * dxdr[2][0] - dxdr[1][0] * dxdr[2][2]);
  drdx[1][1] = inv_jac * (dxdr[0][0] * dxdr[2][2] - dxdr[0][2] * dxdr[2][0]);
  drdx[1][2] = inv_jac * (dxdr[0][2] * dxdr[1][0] - dxdr[0][0] * dxdr[1][2]);
  drdx[2][0] = inv_jac * (dxdr[1][0] * dxdr[2][1] - dxdr[1][1] * dxdr[2][0]);
  drdx[2][1] = inv_jac * (dxdr[0][1] * dxdr[2][0] - dxdr[0][0] * dxdr[2][1]);
  drdx[2][2] = inv_jac * (dxdr[0][0] * dxdr[1][1] - dxdr[0][1] * dxdr[1][0]);

  for (int inode = 0; inode < 8; inode++) {
    dndx[inode][0] = dndr[inode][0] * drdx[0][0] + dndr[inode][1] * drdx[1][0] + dndr[inode][2] * drdx[2][0];
    dndx[inode][1] = dndr[inode][0] * drdx[0][1] + dndr[inode][1] * drdx[1][1] + dndr[inode][2] * drdx[2][1];
    dndx[inode][2] = dndr[inode][0] * drdx[0][2] + dndr[inode][1] * drdx[1][2] + dndr[inode][2] * drdx[2][2];
  }
}

}  // namespace femutil
}  // namespace delfem2

// =======================================================================
// below: tri

// derivative of a shape function of a triangle and constant compornent 
DFM2_INLINE void delfem2::TriDlDx(
    double dldx[][2],
    double const_term[],
    const double p0[],
    const double p1[],
    const double p2[]) {
  const double area = ::delfem2::femutil::TriArea2D(p0, p1, p2);
  const double tmp1 = 0.5 / area;

  const_term[0] = tmp1 * (p1[0] * p2[1] - p2[0] * p1[1]);
  const_term[1] = tmp1 * (p2[0] * p0[1] - p0[0] * p2[1]);
  const_term[2] = tmp1 * (p0[0] * p1[1] - p1[0] * p0[1]);

  dldx[0][0] = tmp1 * (p1[1] - p2[1]);
  dldx[1][0] = tmp1 * (p2[1] - p0[1]);
  dldx[2][0] = tmp1 * (p0[1] - p1[1]);

  dldx[0][1] = tmp1 * (p2[0] - p1[0]);
  dldx[1][1] = tmp1 * (p0[0] - p2[0]);
  dldx[2][1] = tmp1 * (p1[0] - p0[0]);
}

// ------------------------------
// below: tet

DFM2_INLINE double delfem2::femutil::TetVolume3D(
    const double v1[3],
    const double v2[3],
    const double v3[3],
    const double v4[3]) {
  return
      ((v2[0] - v1[0]) * ((v3[1] - v1[1]) * (v4[2] - v1[2]) - (v4[1] - v1[1]) * (v3[2] - v1[2]))
          - (v2[1] - v1[1]) * ((v3[0] - v1[0]) * (v4[2] - v1[2]) - (v4[0] - v1[0]) * (v3[2] - v1[2]))
          + (v2[2] - v1[2]) * ((v3[0] - v1[0]) * (v4[1] - v1[1]) - (v4[0] - v1[0]) * (v3[1] - v1[1]))
      ) * 0.16666666666666666666666666666667;
}

// caluculate Derivative of Area Coord
DFM2_INLINE void delfem2::TetDlDx(
    double dldx[][3],
    double a[],
    const double p0[],
    const double p1[],
    const double p2[],
    const double p3[]) {
  const double vol = femutil::TetVolume3D(p0, p1, p2, p3);
  const double dtmp1 = 1.0 / (vol * 6.0);

  a[0] = +dtmp1 * (p1[0] * (p2[1] * p3[2] - p3[1] * p2[2]) - p1[1] * (p2[0] * p3[2] - p3[0] * p2[2])
      + p1[2] * (p2[0] * p3[1] - p3[0] * p2[1]));
  a[1] = -dtmp1 * (p2[0] * (p3[1] * p0[2] - p0[1] * p3[2]) - p2[1] * (p3[0] * p0[2] - p0[0] * p3[2])
      + p2[2] * (p3[0] * p0[1] - p0[0] * p3[1]));
  a[2] = +dtmp1 * (p3[0] * (p0[1] * p1[2] - p1[1] * p0[2]) - p3[1] * (p0[0] * p1[2] - p1[0] * p0[2])
      + p3[2] * (p0[0] * p1[1] - p1[0] * p0[1]));
  a[3] = -dtmp1 * (p0[0] * (p1[1] * p2[2] - p2[1] * p1[2]) - p0[1] * (p1[0] * p2[2] - p2[0] * p1[2])
      + p0[2] * (p1[0] * p2[1] - p2[0] * p1[1]));

  dldx[0][0] = -dtmp1 * ((p2[1] - p1[1]) * (p3[2] - p1[2]) - (p3[1] - p1[1]) * (p2[2] - p1[2]));
  dldx[0][1] = +dtmp1 * ((p2[0] - p1[0]) * (p3[2] - p1[2]) - (p3[0] - p1[0]) * (p2[2] - p1[2]));
  dldx[0][2] = -dtmp1 * ((p2[0] - p1[0]) * (p3[1] - p1[1]) - (p3[0] - p1[0]) * (p2[1] - p1[1]));

  dldx[1][0] = +dtmp1 * ((p3[1] - p2[1]) * (p0[2] - p2[2]) - (p0[1] - p2[1]) * (p3[2] - p2[2]));
  dldx[1][1] = -dtmp1 * ((p3[0] - p2[0]) * (p0[2] - p2[2]) - (p0[0] - p2[0]) * (p3[2] - p2[2]));
  dldx[1][2] = +dtmp1 * ((p3[0] - p2[0]) * (p0[1] - p2[1]) - (p0[0] - p2[0]) * (p3[1] - p2[1]));

  dldx[2][0] = -dtmp1 * ((p0[1] - p3[1]) * (p1[2] - p3[2]) - (p1[1] - p3[1]) * (p0[2] - p3[2]));
  dldx[2][1] = +dtmp1 * ((p0[0] - p3[0]) * (p1[2] - p3[2]) - (p1[0] - p3[0]) * (p0[2] - p3[2]));
  dldx[2][2] = -dtmp1 * ((p0[0] - p3[0]) * (p1[1] - p3[1]) - (p1[0] - p3[0]) * (p0[1] - p3[1]));

  dldx[3][0] = +dtmp1 * ((p1[1] - p0[1]) * (p2[2] - p0[2]) - (p2[1] - p0[1]) * (p1[2] - p0[2]));
  dldx[3][1] = -dtmp1 * ((p1[0] - p0[0]) * (p2[2] - p0[2]) - (p2[0] - p0[0]) * (p1[2] - p0[2]));
  dldx[3][2] = +dtmp1 * ((p1[0] - p0[0]) * (p2[1] - p0[1]) - (p2[0] - p0[0]) * (p1[1] - p0[1]));

  //  std::cout << dldx[0][0]+dldx[1][0]+dldx[2][0]+dldx[3][0] << std::endl;
  //  std::cout << dldx[0][1]+dldx[1][1]+dldx[2][1]+dldx[3][1] << std::endl;
  //  std::cout << dldx[0][2]+dldx[1][2]+dldx[2][2]+dldx[3][2] << std::endl;

  //  std::cout << a[0]+dldx[0][0]*p0[0]+dldx[0][1]*p0[1]+dldx[0][2]*p0[2] << std::endl;
  //  std::cout << a[1]+dldx[1][0]*p1[0]+dldx[1][1]*p1[1]+dldx[1][2]*p1[2] << std::endl;
  //  std::cout << a[2]+dldx[2][0]*p2[0]+dldx[2][1]*p2[1]+dldx[2][2]*p2[2] << std::endl;
  //  std::cout << a[3]+dldx[3][0]*p3[0]+dldx[3][1]*p3[1]+dldx[3][2]*p3[2] << std::endl;
}

DFM2_INLINE void delfem2::SetEmatConsistentMassTet(
    double eM[4][4][3][3],
    double w0) {
  const double dtmp1 = w0 * 0.05;
  for (int ino = 0; ino < 4; ino++) {
    for (int jno = 0; jno < 4; jno++) {
      eM[ino][jno][0][0] = dtmp1;
      eM[ino][jno][1][1] = dtmp1;
      eM[ino][jno][2][2] = dtmp1;
      eM[ino][jno][0][1] = 0.0;
      eM[ino][jno][0][2] = 0.0;
      eM[ino][jno][1][0] = 0.0;
      eM[ino][jno][1][2] = 0.0;
      eM[ino][jno][2][0] = 0.0;
      eM[ino][jno][2][1] = 0.0;
    }
    eM[ino][ino][0][0] += dtmp1;
    eM[ino][ino][1][1] += dtmp1;
    eM[ino][ino][2][2] += dtmp1;
  }
}

DFM2_INLINE void delfem2::SetEMatLaplaceTet(
    double C[4][4],
    double w0,
    const double dldx[4][3]) {
  for (int ino = 0; ino < 4; ino++) {
    for (int jno = 0; jno < 4; jno++) {
      C[ino][jno] = w0 * (dldx[jno][0] * dldx[ino][0] + dldx[jno][1] * dldx[ino][1] + dldx[jno][2] * dldx[ino][2]);
    }
  }
}

DFM2_INLINE void delfem2::SetEMatLaplaceTet(
    double C[4][4][3][3],
    double w0,
    const double dldx[4][3]) {
  for (int ino = 0; ino < 4; ino++) {
    for (int jno = 0; jno < 4; jno++) {
      const double
          dtmp1 = w0 * (dldx[jno][0] * dldx[ino][0] + dldx[jno][1] * dldx[ino][1] + dldx[jno][2] * dldx[ino][2]);
      C[ino][jno][0][0] = dtmp1;
      C[ino][jno][0][1] = 0.0;
      C[ino][jno][0][2] = 0.0;
      C[ino][jno][1][0] = 0.0;
      C[ino][jno][1][1] = dtmp1;
      C[ino][jno][1][2] = 0.0;
      C[ino][jno][2][0] = 0.0;
      C[ino][jno][2][1] = 0.0;
      C[ino][jno][2][2] = dtmp1;
    }
  }
}

// above: tet
// -----------
// below: vox

DFM2_INLINE void delfem2::ShapeFunc_Vox8(
    const double &r0,
    const double &r1,
    const double &r2,
    const double coords[][3],
    double &detjac,
    double dndx[][3],
    double an[]) {
  an[0] = 0.125 * (1.0 - r0) * (1.0 - r1) * (1.0 - r2);
  an[1] = 0.125 * (1.0 + r0) * (1.0 - r1) * (1.0 - r2);
  an[2] = 0.125 * (1.0 - r0) * (1.0 + r1) * (1.0 - r2);
  an[3] = 0.125 * (1.0 + r0) * (1.0 + r1) * (1.0 - r2);
  an[4] = 0.125 * (1.0 - r0) * (1.0 - r1) * (1.0 + r2);
  an[5] = 0.125 * (1.0 + r0) * (1.0 - r1) * (1.0 + r2);
  an[6] = 0.125 * (1.0 - r0) * (1.0 + r1) * (1.0 + r2);
  an[7] = 0.125 * (1.0 + r0) * (1.0 + r1) * (1.0 + r2);

  double dndr[8][3];
  dndr[0][0] = -0.125 * (1.0 - r1) * (1.0 - r2);
  dndr[1][0] = -dndr[0][0];
  dndr[2][0] = -0.125 * (1.0 + r1) * (1.0 - r2);
  dndr[3][0] = -dndr[2][0];
  dndr[4][0] = -0.125 * (1.0 - r1) * (1.0 + r2);
  dndr[5][0] = -dndr[4][0];
  dndr[6][0] = -0.125 * (1.0 + r1) * (1.0 + r2);
  dndr[7][0] = -dndr[6][0];

  dndr[0][1] = -0.125 * (1.0 - r0) * (1.0 - r2);
  dndr[1][1] = -0.125 * (1.0 + r0) * (1.0 - r2);
  dndr[2][1] = -dndr[0][1];
  dndr[3][1] = -dndr[1][1];
  dndr[4][1] = -0.125 * (1.0 - r0) * (1.0 + r2);
  dndr[5][1] = -0.125 * (1.0 + r0) * (1.0 + r2);
  dndr[6][1] = -dndr[4][1];
  dndr[7][1] = -dndr[5][1];

  dndr[0][2] = -0.125 * (1.0 - r0) * (1.0 - r1);
  dndr[1][2] = -0.125 * (1.0 + r0) * (1.0 - r1);
  dndr[2][2] = -0.125 * (1.0 - r0) * (1.0 + r1);
  dndr[3][2] = -0.125 * (1.0 + r0) * (1.0 + r1);
  dndr[4][2] = -dndr[0][2];
  dndr[5][2] = -dndr[1][2];
  dndr[6][2] = -dndr[2][2];
  dndr[7][2] = -dndr[3][2];

  femutil::HexVox(
      detjac, dndx, an,
      coords, dndr);
}

// vox
// ---------------
// hex

DFM2_INLINE void delfem2::ShapeFunc_Hex8(
    const double &r0,
    const double &r1,
    const double &r2,
    const double coords[][3],
    double &detjac,
    double dndx[][3],
    double an[]) {
  an[0] = 0.125 * (1.0 - r0) * (1.0 - r1) * (1.0 - r2);
  an[1] = 0.125 * (1.0 + r0) * (1.0 - r1) * (1.0 - r2);
  an[2] = 0.125 * (1.0 + r0) * (1.0 + r1) * (1.0 - r2);
  an[3] = 0.125 * (1.0 - r0) * (1.0 + r1) * (1.0 - r2);
  an[4] = 0.125 * (1.0 - r0) * (1.0 - r1) * (1.0 + r2);
  an[5] = 0.125 * (1.0 + r0) * (1.0 - r1) * (1.0 + r2);
  an[6] = 0.125 * (1.0 + r0) * (1.0 + r1) * (1.0 + r2);
  an[7] = 0.125 * (1.0 - r0) * (1.0 + r1) * (1.0 + r2);

  double dndr[8][3];
  dndr[0][0] = -0.125 * (1.0 - r1) * (1.0 - r2);
  dndr[1][0] = +0.125 * (1.0 - r1) * (1.0 - r2);
  dndr[2][0] = +0.125 * (1.0 + r1) * (1.0 - r2);
  dndr[3][0] = -0.125 * (1.0 + r1) * (1.0 - r2);
  dndr[4][0] = -0.125 * (1.0 - r1) * (1.0 + r2);
  dndr[5][0] = +0.125 * (1.0 - r1) * (1.0 + r2);
  dndr[6][0] = +0.125 * (1.0 + r1) * (1.0 + r2);
  dndr[7][0] = -0.125 * (1.0 + r1) * (1.0 + r2);

  dndr[0][1] = -0.125 * (1.0 - r0) * (1.0 - r2);
  dndr[1][1] = -0.125 * (1.0 + r0) * (1.0 - r2);
  dndr[2][1] = +0.125 * (1.0 + r0) * (1.0 - r2);
  dndr[3][1] = +0.125 * (1.0 - r0) * (1.0 - r2);
  dndr[4][1] = -0.125 * (1.0 - r0) * (1.0 + r2);
  dndr[5][1] = -0.125 * (1.0 + r0) * (1.0 + r2);
  dndr[6][1] = +0.125 * (1.0 + r0) * (1.0 + r2);
  dndr[7][1] = +0.125 * (1.0 - r0) * (1.0 + r2);

  dndr[0][2] = -0.125 * (1.0 - r0) * (1.0 - r1);
  dndr[1][2] = -0.125 * (1.0 + r0) * (1.0 - r1);
  dndr[2][2] = -0.125 * (1.0 + r0) * (1.0 + r1);
  dndr[3][2] = -0.125 * (1.0 - r0) * (1.0 + r1);
  dndr[4][2] = +0.125 * (1.0 - r0) * (1.0 - r1);
  dndr[5][2] = +0.125 * (1.0 + r0) * (1.0 - r1);
  dndr[6][2] = +0.125 * (1.0 + r0) * (1.0 + r1);
  dndr[7][2] = +0.125 * (1.0 - r0) * (1.0 + r1);

  femutil::HexVox(
      detjac, dndx, an,
      coords, dndr);
}

DFM2_INLINE double delfem2::DiffShapeFuncAtQuadraturePoint_Hex(
    double dndx[8][3],
    int iGauss,
    int ir1,
    int ir2,
    int ir3,
    const double aP0[8][3]) {
  const double r1 = LineGauss<double>[iGauss][ir1][0];
  const double r2 = LineGauss<double>[iGauss][ir2][0];
  const double r3 = LineGauss<double>[iGauss][ir3][0];
  double an[8], detjac;
  ShapeFunc_Hex8(r1, r2, r3, aP0, detjac, dndx, an);
  return detjac *
      LineGauss<double>[iGauss][ir1][1] *
      LineGauss<double>[iGauss][ir2][1] *
      LineGauss<double>[iGauss][ir3][1];
}