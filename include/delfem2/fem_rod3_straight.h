/*
 * Copyright (c) 2020 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */

#ifndef DFM2_FEM_ROD3_STRAIGHT_H
#define DFM2_FEM_ROD3_STRAIGHT_H

#include "delfem2/dfm2_inline.h"
#include "delfem2/vec3.h"
#include "delfem2/mat3.h"

namespace delfem2 {

template <typename T>
DFM2_INLINE void CdC_Rod3BendStraight(
    T C[3],
    T dC_dP[3][3][3],
    const T vec_pos[3][3]);

DFM2_INLINE double WdWddW_Rod3BendStraight(
    delfem2::CVec3d dW_dP[3],
    delfem2::CMat3d ddW_ddP[3][3],
    const delfem2::CVec3d vec_pos[3]);

} // namespace delfem2

#ifndef DFM2_STATIC_LIBRARY
#  include "delfem2/fem_rod3_straight.cpp"
#endif

#endif  /* DFM2_FEM_ROD3_ENERGY_STRAIGHT_H */
