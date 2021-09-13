/*
 * Copyright (c) 2020 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */

#ifndef DFM2_FEM_DISTANCE3_H
#define DFM2_FEM_DISTANCE3_H

#include "delfem2/dfm2_inline.h"
#include "delfem2/vec3.h"
#include "delfem2/mat3.h"

namespace delfem2 {
/**
 * compute energy and its graident and hessian for the 3D spring
 * @param[out] dw_dp
 * @param[out] ddw_ddp
 * @param[in] stiffness stiffness
 * @param[in] vtx_xyz_ini current positions
 * @param[in] edge_length_ini rest length
 * @return energy
 */
DFM2_INLINE double WdWddW_SquareLengthLineseg3D(
    CVec3d dw_dp[2],
    CMat3d ddw_ddp[2][2],
    //
    double stiffness,
    const CVec3d vtx_xyz_ini[2],
    double edge_length_ini);

template <typename T>
DFM2_INLINE void CdC_SquareLengthLineseg3D(
    T& c,
    T dc_dpos[2][3],
    const T pos_xyz[2][3],
    T length_ini);

} // namespace delfem2

#ifndef DFM2_STATIC_LIBRARY
#  include "delfem2/fem_distance3.cpp"
#endif

#endif  /* DFM2_FEM_ROD3_ENERGY_STRAIGHT_H */
