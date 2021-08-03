/*
 * Copyright (c) 2020 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */

#ifndef DFM2_FEMROD_H
#define DFM2_FEMROD_H

#include "delfem2/dfm2_inline.h"
#include "delfem2/vec3.h"
#include "delfem2/mat3.h"
#include "delfem2/lsmats.h"

namespace delfem2 {

DFM2_INLINE void WdWddW_Rod2(
    double &W,
    double dW[3][2],
    double ddW[3][3][2][2],
    const double ap[3][2],
    const double aL[2],
    double stiff_stretch01,
    double stiff_stretch12,
    double stiff1_bend);


DFM2_INLINE double WdWddW_Rod3BendStraight(
    delfem2::CVec3d dW_dP[3],
    delfem2::CMat3d ddW_ddP[3][3],
    const delfem2::CVec3d vec_pos[3]);

/**
 * @brief energy W and its derivative dW and second derivative ddW
 * where W = a^T R(dn) b(theta)
 */
DFM2_INLINE void RodFrameTrans(
    CVec3d frm[3],
    const CVec3d &S0,
    const CVec3d &V01,
    const CVec3d &du,
    double dtheta);

/**
 * @param[out] dF_dv how i-th frame axis moves w.r.t the movement of the edge
 * @param[out] dF_dt how i-th frame axis moves w.r.t the rotation along 2nd axis
 * @param[in] l01 length of the edge
 * @param[in] Frm frame axis vectors
 */
DFM2_INLINE void DiffFrameRod(
    CMat3d dF_dv[3],
    CVec3d dF_dt[3],
    //
    double l01,
    const CVec3d Frm[3]);

/**
 * @brief second derivative (ddW) of W, where W = Q^T Frm[iaxis]
 */
DFM2_INLINE void DifDifFrameRod(
    CMat3d &ddW_ddv,
    CVec3d &ddW_dvdt,
    double &ddW_dtt,
    //
    unsigned int iaxis,
    double l01,
    const CVec3d &Q,
    const CVec3d Frm[3]);

DFM2_INLINE double WdWddW_DotFrame(
    CVec3d dV_dP[3],
    double dV_dt[2],
    CMat3d ddV_ddP[3][3],
    CVec3d ddV_dtdP[2][3],
    double ddV_ddt[2][2],
    //
    const CVec3d P[3],
    const CVec3d S[2],
    const double off[3]);

DFM2_INLINE CVec3d Darboux_Rod(
    const CVec3d P[3],
    const CVec3d S[2]);

/**
 * @param[in] P points of rod
 * @param[in] S director vectors on the edges
 * @param[in] daboux0 Daboux vector in the rest shape
 * @param[in] is_eaxct whether the hessian is exact or not
 */
DFM2_INLINE double WdWddW_Rod(
    CVec3d dW_dP[3],
    double dW_dt[2],
    CMat3d ddW_ddP[3][3],
    CVec3d ddW_dtdP[2][3],
    double ddW_ddt[2][2],
    const double stiff_bendtwist[3],
    const CVec3d P[3],
    const CVec3d S[2],
    const CVec3d &darboux0,
    bool is_exact);

/**
 * compute energy and its graident and hessian for the 3D spring
 * @param[out] dW_dP
 * @param[out] ddW_ddP
 * @param[in] stiff stiffness
 * @param[in] P current positions
 * @param[in] L0 rest length
 * @return energy
 */
DFM2_INLINE double WdWddW_SquareLengthLineseg3D(
    CVec3d dW_dP[2],
    CMat3d ddW_ddP[2][2],
    //
    const double stiff,
    const CVec3d P[2],
    double L0);

// ------------------------

DFM2_INLINE void Solve_DispRotSeparate(
    std::vector<CVec3d> &aP,
    std::vector<CVec3d> &aS,
    CMatrixSparse<double> &mats,
    const double stiff_stretch,
    const double stiff_bendtwist[3],
    const std::vector<CVec3d> &aP0,
    const std::vector<CVec3d> &aDarboux0,
    const std::vector<unsigned int> &aElemSeg,
    const std::vector<unsigned int> &aElemRod,
    const std::vector<int> &aBCFlag);


// --------------
// below: rod hair

DFM2_INLINE void ParallelTransport_RodHair(
    std::vector<CVec3d> &aP0,
    std::vector<CVec3d> &aS0,
    const std::vector<unsigned int> &aIP_HairRoot);

DFM2_INLINE void MakeBCFlag_RodHair(
    std::vector<int> &aBCFlag,
    const std::vector<unsigned int> &aIP_HairRoot);

DFM2_INLINE void MakeSparseMatrix_RodHair(
    CMatrixSparse<double> &mats,
    const std::vector<unsigned int> &aIP_HairRoot,
    unsigned int ndof_par_node);

DFM2_INLINE void MakeDirectorOrthogonal_RodHair(
    std::vector<CVec3d> &aS,
    const std::vector<CVec3d> &aP);

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
DFM2_INLINE double MergeLinSys_Hair(
    std::vector<double> &vec_r,
    CMatrixSparse<double> &mats,
    const double stiff_stretch,
    const double stiff_bendtwist[3],
    const std::vector<unsigned int> &aIP_HairRoot,
    const std::vector<CVec3d> &aP,
    const std::vector<CVec3d> &aS,
    const std::vector<CVec3d> &aP0,
    const std::vector<CVec3d> &aS0);

DFM2_INLINE void UpdateSolutionHair(
    std::vector<CVec3d> &aP,
    std::vector<CVec3d> &aS,
    const std::vector<double> &vec_x,
    const std::vector<unsigned int> &aIP_HairRoot,
    const std::vector<int> &aBCFlag);

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
DFM2_INLINE void Solve_RodHair(
    std::vector<CVec3d> &aP,
    std::vector<CVec3d> &aS,
    CMatrixSparse<double> &mats,
    const double stiff_stretch,
    const double stiff_bendtwist[3],
    double mdtt,
    const std::vector<CVec3d> &aP0,
    const std::vector<CVec3d> &aS0,
    const std::vector<int> &aBCFlag,
    const std::vector<unsigned int> &aIP_HairRoot);

class CContactHair {
 public:
  unsigned int ip0, ip1;
  double s;
  unsigned int iq0, iq1;
  double t;
  CVec3d norm;
 public:
  CVec3d Direction(const std::vector<CVec3d> &aP) const {
    return (1 - s) * aP[ip0] + s * aP[ip1] - (1 - t) * aP[iq0] - t * aP[iq1];
  }
};

DFM2_INLINE void Solve_RodHairContact(
    std::vector<CVec3d> &aP,
    std::vector<CVec3d> &aS,
    CMatrixSparse<double> &mats,
    const double stiff_stretch,
    const double stiff_bendtwist[3],
    double mdtt,
    const std::vector<CVec3d> &aPt0,
    const std::vector<CVec3d> &aP0,
    const std::vector<CVec3d> &aS0,
    const std::vector<int> &aBCFlag,
    const std::vector<unsigned int> &aIP_HairRoot,
    const double clearance,
    const double stiff_contact,
    const std::vector<CContactHair> &aContact);

} // namespace delfem2

#ifdef DFM2_HEADER_ONLY
#include "delfem2/femrod.cpp"
#endif

#endif
