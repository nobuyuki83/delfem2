/*
 * Copyright (c) 2020 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */

#ifndef DFM2_FEM_ROD3_DARBOUX_H
#define DFM2_FEM_ROD3_DARBOUX_H

#include "delfem2/dfm2_inline.h"
#include "delfem2/vec3.h"
#include "delfem2/mat3.h"

namespace delfem2 {

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

} // namespace delfem2

#ifndef DFM2_STATIC_LIBRARY
#  include "delfem2/fem_rod3_darboux.cpp"
#endif

#endif  /* DFM2_FEM_ROD3_ENERGY_DARBOUX_H */
