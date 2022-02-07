/*
 * Copyright (c) 2019-2021 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */

#include "delfem2/pbd_geo3.h"

#include "delfem2/mat2.h"
#include "delfem2/svd3.h"
#include "delfem2/geo3_v23m34q.h"
#include "delfem2/geo_vec3.h"
#include "delfem2/geo_mat3.h"

// =======================================

DFM2_INLINE void delfem2::PBD_Post(
    std::vector<double> &aXYZ,
    std::vector<double> &aUVW,
    double dt,
    const std::vector<double> &aXYZt,
    const std::vector<int> &aBCFlag) {
  const size_t np = aXYZ.size() / 3;
  assert(aBCFlag.size() == np);
  assert(aXYZ.size() == np * 3);
  assert(aUVW.size() == np * 3);
  for (unsigned int ip = 0; ip < np; ++ip) {
    if (aBCFlag[ip] != 0) continue;
    aUVW[ip * 3 + 0] = (aXYZt[ip * 3 + 0] - aXYZ[ip * 3 + 0]) / dt;
    aUVW[ip * 3 + 1] = (aXYZt[ip * 3 + 1] - aXYZ[ip * 3 + 1]) / dt;
    aUVW[ip * 3 + 2] = (aXYZt[ip * 3 + 2] - aXYZ[ip * 3 + 2]) / dt;
    aXYZ[ip * 3 + 0] = aXYZt[ip * 3 + 0];
    aXYZ[ip * 3 + 1] = aXYZt[ip * 3 + 1];
    aXYZ[ip * 3 + 2] = aXYZt[ip * 3 + 2];
  }
}

DFM2_INLINE void delfem2::PBD_Pre3D(
    std::vector<double> &aXYZt,
    double dt,
    const double gravity[3],
    const std::vector<double> &aXYZ,
    const std::vector<double> &aUVW,
    const std::vector<int> &aBCFlag) {
  const size_t np = aXYZ.size() / 3;
  assert(aBCFlag.size() == np);
  assert(aUVW.size() == np * 3);
  aXYZt.resize(np * 3);
  for (unsigned int ip = 0; ip < np; ++ip) {
    if (aBCFlag[ip] != 0) {
      aXYZt[ip * 3 + 0] = aXYZ[ip * 3 + 0];
      aXYZt[ip * 3 + 1] = aXYZ[ip * 3 + 1];
      aXYZt[ip * 3 + 2] = aXYZ[ip * 3 + 2];
      continue;
    }
    aXYZt[ip * 3 + 0] = aXYZ[ip * 3 + 0] + dt * aUVW[ip * 3 + 0] + dt * dt * gravity[0];
    aXYZt[ip * 3 + 1] = aXYZ[ip * 3 + 1] + dt * aUVW[ip * 3 + 1] + dt * dt * gravity[1];
    aXYZt[ip * 3 + 2] = aXYZ[ip * 3 + 2] + dt * aUVW[ip * 3 + 2] + dt * dt * gravity[2];
  }
}

DFM2_INLINE void delfem2::PBD_Update_Const3(
    double *aXYZt,
    const int np,
    const int ndim,
    const double *m,
    const double *C,
    const double *dCdp,
    const unsigned int *aIP,
    double ratio) {
  std::vector<double> mi(np);
  for (int ip = 0; ip < np; ++ip) { mi[ip] = 1.0 / m[ip]; }
  //
  const int nc = 3;
  std::vector<double> MinvC(nc * np * ndim);
  for (int ic = 0; ic < nc; ++ic) {
    for (int ine = 0; ine < np; ++ine) {
      MinvC[ic * np * ndim + ine * 3 + 0] = dCdp[ic * np * ndim + ine * ndim + 0] * mi[ine];
      MinvC[ic * np * ndim + ine * 3 + 1] = dCdp[ic * np * ndim + ine * ndim + 1] * mi[ine];
      MinvC[ic * np * ndim + ine * 3 + 2] = dCdp[ic * np * ndim + ine * ndim + 2] * mi[ine];
    }
  }
  double A[nc * nc] = {0, 0, 0, 0, 0, 0, 0, 0, 0};
  for (int i = 0; i < np * ndim; ++i) {
    A[0 * 3 + 0] += MinvC[0 * np * ndim + i] * dCdp[0 * np * ndim + i];
    A[0 * 3 + 1] += MinvC[0 * np * ndim + i] * dCdp[1 * np * ndim + i];
    A[0 * 3 + 2] += MinvC[0 * np * ndim + i] * dCdp[2 * np * ndim + i];
    A[1 * 3 + 0] += MinvC[1 * np * ndim + i] * dCdp[0 * np * ndim + i];
    A[1 * 3 + 1] += MinvC[1 * np * ndim + i] * dCdp[1 * np * ndim + i];
    A[1 * 3 + 2] += MinvC[1 * np * ndim + i] * dCdp[2 * np * ndim + i];
    A[2 * 3 + 0] += MinvC[2 * np * ndim + i] * dCdp[0 * np * ndim + i];
    A[2 * 3 + 1] += MinvC[2 * np * ndim + i] * dCdp[1 * np * ndim + i];
    A[2 * 3 + 2] += MinvC[2 * np * ndim + i] * dCdp[2 * np * ndim + i];
  }
  double Ainv[nc * nc];
  Inverse_Mat3(Ainv, A);
  double lmd[nc];
  MatVec3(lmd, Ainv, C);
  for (int ine = 0; ine < np; ++ine) {
    const unsigned int ip0 = aIP[ine];
    for (int ic = 0; ic < nc; ++ic) {
      for (int idim = 0; idim < ndim; ++idim) {
        aXYZt[ip0 * 3 + idim] -= ratio * MinvC[ic * np * ndim + ine * ndim + idim] * lmd[ic];
      }
    }

  }
}

DFM2_INLINE void delfem2::PBD_ConstProj_Rigid3D(
    double *aXYZt,
    double stiffness,
    const int *clstr_ind,
    int nclstr_ind,
    const int *clstr,
    [[maybe_unused]] int nclstr0,
    [[maybe_unused]] const double *aXYZ0,
    [[maybe_unused]] int nXYZ0) {
  const int nclstr = nclstr_ind - 1;
  for (int iclstr = 0; iclstr < nclstr; ++iclstr) {
    CVec3d pc(0, 0, 0), qc(0, 0, 0);
    for (int iip = clstr_ind[iclstr]; iip < clstr_ind[iclstr + 1]; iip++) {
      const int ip = clstr[iip];
      assert(ip < nclstr0);
      qc += CVec3d(aXYZ0[ip * 3 + 0], aXYZ0[ip * 3 + 1], aXYZ0[ip * 3 + 2]);
      pc += CVec3d(aXYZt[ip * 3 + 0], aXYZt[ip * 3 + 1], aXYZt[ip * 3 + 2]);
    }
    qc /= (clstr_ind[iclstr + 1] - clstr_ind[iclstr]);
    pc /= (clstr_ind[iclstr + 1] - clstr_ind[iclstr]);

    double A[9] = {0, 0, 0, 0, 0, 0, 0, 0, 0};
    for (int iip = clstr_ind[iclstr]; iip < clstr_ind[iclstr + 1]; iip++) {
      const int ip = clstr[iip];
      const CVec3d dq = CVec3d(aXYZ0[ip * 3 + 0], aXYZ0[ip * 3 + 1], aXYZ0[ip * 3 + 2]) - qc; // undeform
      const CVec3d dp = CVec3d(aXYZt[ip * 3 + 0], aXYZt[ip * 3 + 1], aXYZ0[ip * 3 + 2]) - pc; // deform
      A[0 * 3 + 0] += dp[0] * dq[0];
      A[0 * 3 + 1] += dp[0] * dq[1];
      A[0 * 3 + 2] += dp[0] * dq[2];
      A[1 * 3 + 0] += dp[1] * dq[0];
      A[1 * 3 + 1] += dp[1] * dq[1];
      A[1 * 3 + 2] += dp[1] * dq[2];
      A[2 * 3 + 0] += dp[2] * dq[0];
      A[2 * 3 + 1] += dp[2] * dq[1];
      A[2 * 3 + 2] += dp[2] * dq[2];
    }
    double R[9];
    GetRotPolarDecomp(R, A, 20);
    //    std::cout << R[0] << " " << R[1] << " " << R[2] << " " << R[3] << std::endl;

    for (int iip = clstr_ind[iclstr]; iip < clstr_ind[iclstr + 1]; iip++) {
      const int ip = clstr[iip];
      CVec3d dq = CVec3d(aXYZ0[ip * 3 + 0], aXYZ0[ip * 3 + 1], aXYZ0[ip * 3 + 2]) - qc;
      CVec3d pg = pc + Mat3Vec(R, dq); // goal position
      CVec3d pg2 = stiffness * pg + (1 - stiffness) * CVec3d(aXYZt[ip * 3 + 0], aXYZt[ip * 3 + 1], aXYZt[ip * 3 + 2]);
      aXYZt[ip * 3 + 0] = pg2.p[0];
      aXYZt[ip * 3 + 1] = pg2.p[1];
      aXYZt[ip * 3 + 2] = pg2.p[2];
    }
  }
}

DFM2_INLINE void delfem2::PBD_ConstProj_Rigid2D(
    double *aXYt,
    double stiffness,
    const unsigned int *clstr_ind,
    size_t nclstr_ind,
    const unsigned int *clstr,
    [[maybe_unused]] size_t nclstr0,
    const double *aXY0,
    [[maybe_unused]] size_t nXY0) {
  const size_t nclstr = nclstr_ind - 1;
  for (unsigned int iclstr = 0; iclstr < nclstr; ++iclstr) {
    CVec2d pc(0, 0), qc(0, 0);
    for (unsigned int iip = clstr_ind[iclstr]; iip < clstr_ind[iclstr + 1]; iip++) {
      assert(iip < nclstr0);
      const unsigned int ip = clstr[iip];
      assert(ip < nXY0);
      qc += CVec2d(aXY0[ip * 2 + 0], aXY0[ip * 2 + 1]);
      pc += CVec2d(aXYt[ip * 2 + 0], aXYt[ip * 2 + 1]);
    }
    qc /= (clstr_ind[iclstr + 1] - clstr_ind[iclstr]);
    pc /= (clstr_ind[iclstr + 1] - clstr_ind[iclstr]);

    double A[4] = {0, 0, 0, 0};
    for (unsigned int iip = clstr_ind[iclstr]; iip < clstr_ind[iclstr + 1]; iip++) {
      const unsigned int ip = clstr[iip];
      const CVec2d dq = CVec2d(aXY0[ip * 2 + 0], aXY0[ip * 2 + 1]) - qc; // undeform
      const CVec2d dp = CVec2d(aXYt[ip * 2 + 0], aXYt[ip * 2 + 1]) - pc; // deform
      A[0 * 2 + 0] += dp[0] * dq[0];
      A[0 * 2 + 1] += dp[0] * dq[1];
      A[1 * 2 + 0] += dp[1] * dq[0];
      A[1 * 2 + 1] += dp[1] * dq[1];
    }
    double R[4];
    RotationalComponentOfMatrix2(R, A);
    //    std::cout << R[0] << " " << R[1] << " " << R[2] << " " << R[3] << std::endl;

    for (unsigned int iip = clstr_ind[iclstr]; iip < clstr_ind[iclstr + 1]; iip++) {
      const unsigned int ip = clstr[iip];
      CVec2d dq = CVec2d(aXY0[ip * 2 + 0], aXY0[ip * 2 + 1]) - qc;
      CVec2d pg = pc + Mat2Vec(R, dq); // goal position
      CVec2d pg2 = stiffness * pg + (1 - stiffness) * CVec2d(aXYt[ip * 2 + 0], aXYt[ip * 2 + 1]);
      aXYt[ip * 2 + 0] = pg2.x;
      aXYt[ip * 2 + 1] = pg2.y;
    }
  }
}


// -------------------------

DFM2_INLINE void delfem2::PBD_Seam(
  double *aXYZt,
  [[maybe_unused]] size_t nXYZ,
  const unsigned int *aLine,
  size_t nline) {
  for (unsigned int il = 0; il < nline; ++il) {
    const unsigned int ip0 = aLine[il * 2 + 0];
    const unsigned int ip1 = aLine[il * 2 + 1];
    const double p[2][3] = {
      {aXYZt[ip0 * 3 + 0], aXYZt[ip0 * 3 + 1], aXYZt[ip0 * 3 + 2]},
      {aXYZt[ip1 * 3 + 0], aXYZt[ip1 * 3 + 1], aXYZt[ip1 * 3 + 2]}};
    double d0 = Distance3(p[0], p[1]);
    double dLen = 0.01;
    if (d0 > dLen) {
      double n01[3] = {p[1][0] - p[0][0], p[1][1] - p[0][1], p[1][2] - p[0][2]};
      const double l01 = Length3(n01);
      const double invl01 = 1.0 / l01;
      n01[0] *= invl01;
      n01[1] *= invl01;
      n01[2] *= invl01;
      aXYZt[ip0 * 3 + 0] += n01[0] * dLen * 0.5;
      aXYZt[ip0 * 3 + 1] += n01[1] * dLen * 0.5;
      aXYZt[ip0 * 3 + 2] += n01[2] * dLen * 0.5;
      aXYZt[ip1 * 3 + 0] -= n01[0] * dLen * 0.5;
      aXYZt[ip1 * 3 + 1] -= n01[1] * dLen * 0.5;
      aXYZt[ip1 * 3 + 2] -= n01[2] * dLen * 0.5;
    } else {
      aXYZt[ip0 * 3 + 0] = (p[0][0] + p[1][0]) * 0.5;
      aXYZt[ip0 * 3 + 1] = (p[0][1] + p[1][1]) * 0.5;
      aXYZt[ip0 * 3 + 2] = (p[0][2] + p[1][2]) * 0.5;
      aXYZt[ip1 * 3 + 0] = (p[0][0] + p[1][0]) * 0.5;
      aXYZt[ip1 * 3 + 1] = (p[0][1] + p[1][1]) * 0.5;
      aXYZt[ip1 * 3 + 2] = (p[0][2] + p[1][2]) * 0.5;
    }
  }
}

