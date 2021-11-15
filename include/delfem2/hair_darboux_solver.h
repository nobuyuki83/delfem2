/*
 * Copyright (c) 2020 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */

#ifndef DFM2_HAIR_DARBOUX_SOLVER_H
#define DFM2_HAIR_DARBOUX_SOLVER_H

#include "delfem2/dfm2_inline.h"
#include "delfem2/vec3.h"
#include "delfem2/mat3.h"
#include "delfem2/lsmats.h"
#include "delfem2/hair_darboux_util.h"
#include "delfem2/fem_distance3.h"
#include "delfem2/fem_rod3_darboux.h"

namespace delfem2 {

/**
 * Merge linear system for hair
 * @param vec_r
 * @param[in,out] mats sparse matrix
 * @param[in] stiff_stretch
 * @param[in] stiff_bendtwist
 * @param[in] aIP_HairRoot
 * @param[in] aP
 * @param[in] aS
 * @param[in] aP0 vertex position of rest shape
 * @param[in] aS0
 * @return
 */
template <class BLOCK_LINEAR_SOLVER>
double MergeLinSys_Hair(
    BLOCK_LINEAR_SOLVER &sparse,
    const double stiff_stretch,
    const double stiff_bendtwist[3],
    const std::vector<unsigned int> &aIP_HairRoot,
    const std::vector<delfem2::CVec3d> &aP,
    const std::vector<delfem2::CVec3d> &aS,
    const std::vector<delfem2::CVec3d> &aP0,
    const std::vector<delfem2::CVec3d> &aS0) {
  using namespace delfem2;
  double W = 0;
  for (unsigned int ihair = 0; ihair < aIP_HairRoot.size() - 1; ++ihair) {
    unsigned int ips = aIP_HairRoot[ihair];
    unsigned int ns = aIP_HairRoot[ihair + 1] - aIP_HairRoot[ihair] - 1;
    for (unsigned int is = 0; is < ns; ++is) {
      const unsigned int ip0 = ips + is + 0;
      const unsigned int ip1 = ips + is + 1;
      const unsigned int aINoel[2] = {ip0, ip1};
      const double L0 = (aP0[ip0] - aP0[ip1]).norm();
      const CVec3d aPE[2] = {aP[ip0], aP[ip1]};
      // --------------
      CVec3d dW_dP[2];
      CMat3d ddW_ddP[2][2];
      W += WdWddW_SquareLengthLineseg3D(
          dW_dP, ddW_ddP,
          stiff_stretch, aPE, L0);
      {
        double eM[2][2][4][4];
        for (unsigned int i = 0; i < 2 * 2 * 4 * 4; ++i) { (&eM[0][0][0][0])[i] = 0.0; }
        for (int in = 0; in < 2; ++in) {
          for (int jn = 0; jn < 2; ++jn) {
            ddW_ddP[in][jn].CopyToMat4(&eM[in][jn][0][0]);
          }
        }
        sparse.template Merge<2, 2, 4, 4>(aINoel, aINoel, eM);
      }
      for (int in = 0; in < 2; in++) {
        const unsigned int ip = aINoel[in];
        sparse.vec_r[ip * 4 + 0] -= dW_dP[in].x;
        sparse.vec_r[ip * 4 + 1] -= dW_dP[in].y;
        sparse.vec_r[ip * 4 + 2] -= dW_dP[in].z;
      }
    }
  }
  // --------------------------
  for (unsigned int ihair = 0; ihair < aIP_HairRoot.size() - 1; ++ihair) {
    const unsigned int ips = aIP_HairRoot[ihair];
    const unsigned int nr = aIP_HairRoot[ihair + 1] - aIP_HairRoot[ihair] - 2;
    for (unsigned int ir = 0; ir < nr; ++ir) {
      const unsigned int ip0 = ips + ir + 0;
      const unsigned int ip1 = ips + ir + 1;
      const unsigned int ip2 = ips + ir + 2;
      const unsigned int aINoel[3] = {ip0, ip1, ip2};
      const CVec3d aPE[3] = {aP[ip0], aP[ip1], aP[ip2]};
      const CVec3d aSE[2] = {aS[ip0], aS[ip1]};
      CVec3d Darboux0;
      { // Daboux vector of initial configuration
        const CVec3d aPE0[3] = {aP0[ip0], aP0[ip1], aP0[ip2]};
        const CVec3d aSE0[2] = {aS0[ip0], aS0[ip1]};
        Darboux0 = Darboux_Rod(aPE0, aSE0);
      }
      // ------
      CVec3d dW_dP[3];
      double dW_dt[2];
      CMat3d ddW_ddP[3][3];
      CVec3d ddW_dtdP[2][3];
      double ddW_ddt[2][2];
      W += WdWddW_Rod(
          dW_dP, dW_dt, ddW_ddP, ddW_dtdP, ddW_ddt,
          stiff_bendtwist,
          aPE, aSE, Darboux0, false);
      {
        double eM[3][3][4][4];
        for (unsigned int i = 0; i < 3 * 3 * 4 * 4; ++i) { (&eM[0][0][0][0])[i] = 0.0; }
        for (int in = 0; in < 3; ++in) {
          for (int jn = 0; jn < 3; ++jn) {
            ddW_ddP[in][jn].CopyToMat4(&eM[in][jn][0][0]);
          }
        }
        for (int in = 0; in < 3; ++in) {
          for (int jn = 0; jn < 2; ++jn) {
            eM[in][jn][0][3] = eM[jn][in][3][0] = ddW_dtdP[jn][in].x;
            eM[in][jn][1][3] = eM[jn][in][3][1] = ddW_dtdP[jn][in].y;
            eM[in][jn][2][3] = eM[jn][in][3][2] = ddW_dtdP[jn][in].z;
          }
        }
        for (int in = 0; in < 2; ++in) {
          for (int jn = 0; jn < 2; ++jn) {
            eM[in][jn][3][3] = ddW_ddt[in][jn];
          }
        }
        sparse.template Merge<3, 3, 4, 4>(aINoel, aINoel, eM);
      }
      {
        for (int ino = 0; ino < 3; ino++) {
          const unsigned int ip = aINoel[ino];
          sparse.vec_r[ip * 4 + 0] -= dW_dP[ino].x;
          sparse.vec_r[ip * 4 + 1] -= dW_dP[ino].y;
          sparse.vec_r[ip * 4 + 2] -= dW_dP[ino].z;
        }
        for (int in = 0; in < 2; in++) {
          const unsigned int in0 = aINoel[in];
          sparse.vec_r[in0 * 4 + 3] -= dW_dt[in];
        }
      }
    }
  }
  return W;
}

/**
 * @brief static minimization of the deformation energy
 * @param aP (in&out) position of the vertices of the rods
 * @param aS (in&out) director vectors
 * @param mats (in&out) sparse matrix
 * @param mdtt (in) mass divided by square of timestep (mass/dt/dt)
 * @param aP0 (in) initial positions of the vertices of the rods
 * @param aS0 (in) initial darboux vectors
 * @param aBCFlag (in) boundary condition flag. Non zero value means fixed value
 * @param aIP_HairRoot (in) indeces of the root points
 */
template <class BLOCK_LINEAR_SOLVER>
DFM2_INLINE void Solve_RodHair(
    std::vector<CVec3d> &aP,
    std::vector<CVec3d> &aS,
    BLOCK_LINEAR_SOLVER &sparse,
    const double stiff_stretch,
    const double stiff_bendtwist[3],
    double mdtt,
    const std::vector<CVec3d> &aP0,
    const std::vector<CVec3d> &aS0,
    const std::vector<unsigned int> &aIP_HairRoot) {
  assert( sparse.ndim() == 4 && sparse.nblk() == aP.size() );
  assert(aP0.size() == aP.size());
  assert(aS0.size() == aP.size());
  assert(aP.size() == aP.size());
  assert(aS.size() == aP.size());
  sparse.BeginMerge();
  double W = MergeLinSys_Hair(
      sparse,
      stiff_stretch, stiff_bendtwist,
      aIP_HairRoot, aP, aS, aP0, aS0);
  for (unsigned int ip = 0; ip < aP.size(); ++ip) {
    sparse.AddValueToDiagonal(ip,0,mdtt);
    sparse.AddValueToDiagonal(ip,1,mdtt);
    sparse.AddValueToDiagonal(ip,2,mdtt);
  }
  std::cout << "energy:" << W << std::endl;
  //    std::cout << "sym: " << CheckSymmetry(mats) << std::endl;
  sparse.Solve();
  UpdateSolutionHair(
      aP, aS,
      sparse.vec_x, aIP_HairRoot, sparse.dof_bcflag);
}

} // namespace delfem2

#endif  /* DFM2_HAIR_DARBOUX_SOLVER_H */
