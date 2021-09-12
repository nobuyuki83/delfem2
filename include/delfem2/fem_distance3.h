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
    double stiff,
    const CVec3d P[2],
    double L0);

template <typename T>
DFM2_INLINE T WdW_SquareLengthLineseg3D(
    T dW_dP[2][3],
    T stiff,
    const T p[2][3],
    T L0);

} // namespace delfem2

#ifndef DFM2_STATIC_LIBRARY
#  include "delfem2/fem_distance3.cpp"
#endif

#endif  /* DFM2_FEM_ROD3_ENERGY_STRAIGHT_H */
