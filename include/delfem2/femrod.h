//
//  femrod.h
//  examples_glfwold_hdronly
//
//  Created by Nobuyuki Umetani on 2020-05-01.
//

#ifndef DFM2_FEMROD_H
#define DFM2_FEMROD_H

#include "delfem2/dfm2_inline.h"
#include "delfem2/vec3.h"
#include "delfem2/mat3.h"
#include "delfem2/mats.h"

namespace delfem2 {

/**
 * @brief energy W and its derivative dW and second derivative ddW
 * where W = a^T R(dn) b(theta)
 */
DFM2_INLINE void RodFrameTrans
(CVec3d frm[3],
 const CVec3d& S0,
 const CVec3d& V01,
 const CVec3d& du,
 double dtheta);

/**
 * @param dF_dv (out) how i-th frame axis moves w.r.t the movement of the edge
 * @param dF_dt (out) how i-th frame axis moves w.r.t the rotation along 2nd axis
 * @param l01 (in) length of the edge
 * @param Frm (in) frame axis vectors
 */
DFM2_INLINE void DiffFrameRod
(CMat3d dF_dv[3],
 CVec3d dF_dt[3],
 //
 double l01,
 const CVec3d Frm[3]);

/**
 * @brief second derivative (ddW) of W, where W = Q^T Frm[iaxis]
 */
DFM2_INLINE void DifDifFrameRod
(CMat3d& ddW_ddv,
 CVec3d& ddW_dvdt,
 double& ddW_dtt,
 //
 unsigned int iaxis,
 double l01,
 const CVec3d& Q,
 const CVec3d Frm[3]);

DFM2_INLINE double WdWddW_DotFrame
(CVec3d dV_dP[3],
 double dV_dt[2],
 CMat3d ddV_ddP[3][3],
 CVec3d ddV_dtdP[2][3],
 double ddV_ddt[2][2],
 //
 const CVec3d P[3],
 const CVec3d S[2],
 const double off[3]);

DFM2_INLINE CVec3d Darboux_Rod
(const CVec3d P[3],
 const CVec3d S[2]);

/**
 * @param P[3] (in) point of rod
 * @param S[2] (in) director vectors on the edges
 * @param off[3] (in) Daboux vector in the material frame
 * @param is_eaxct
 */
DFM2_INLINE double WdWddW_Rod
(CVec3d dW_dP[3],
 double dW_dt[2],
 CMat3d ddW_ddP[3][3],
 CVec3d ddW_dtdP[2][3],
 double ddW_ddt[2][2],
 //
 const CVec3d P[3],
 const CVec3d S[2],
 const CVec3d& darboux0,
 bool is_exact);

DFM2_INLINE double WdWddW_SquareLengthLineseg3D
(CVec3d dW_dP[2],
 CMat3d ddW_ddP[2][2],
 //
 const CVec3d P[2],
 double L0);

// ------------------------

DFM2_INLINE void Solve_DispRotSeparate
 (std::vector<CVec3d>& aP,
  std::vector<CVec3d>& aS,
  CMatrixSparse<double>& mats,
  const std::vector<CVec3d>& aP0,
  const std::vector<CVec3d>& aDarboux0,
  const std::vector<unsigned int>& aElemSeg,
  const std::vector<unsigned int>& aElemRod,
  const std::vector<int>& aBCFlag);

DFM2_INLINE void Solve_DispRotCombined
 (std::vector<CVec3d>& aP,
  std::vector<CVec3d>& aS,
  CMatrixSparse<double>& mats,
  const std::vector<CVec3d>& aP0,
  const std::vector<CVec3d>& aDarboux0,
  const std::vector<int>& aBCFlag,
  const std::vector<unsigned int>& aIP_HairRoot);

} // namespace delfem2

#ifdef DFM2_HEADER_ONLY
#  include "delfem2/femrod.cpp"
#endif

#endif /* femrod_h */
