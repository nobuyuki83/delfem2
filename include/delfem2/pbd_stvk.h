/*
 * Copyright (c) 2019 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */


#ifndef DFM2_PBD_STVK_H
#define DFM2_PBD_STVK_H

#include "delfem2/dfm2_inline.h"
#include "delfem2/fem_stvk.h"
#include "delfem2/pbd_geo3.h"

// ------------------------------

namespace delfem2 {

/**
 *
 * @param aXYZt
 * @param nXYZ
 * @param[in] aETri triangle index
 * @param[in] aVec2 2D positions
 */
template<typename Array2ui, typename Array2d>
void PBD_TriStrain(
    double *aXYZt,
    size_t nXYZ,
    const Array2ui &tri_vtx,
    const Array2d &vtx_xy) {
  for (const auto &etri: tri_vtx) {
    unsigned int i0 = etri[0];
    unsigned int i1 = etri[1];
    unsigned int i2 = etri[2];
    const double P[3][2] = { // rest shape coordinate
        {vtx_xy[i0][0], vtx_xy[i0][1]},
        {vtx_xy[i1][0], vtx_xy[i1][1]},
        {vtx_xy[i2][0], vtx_xy[i2][1]}};
    const double p[3][3] = {
        {aXYZt[i0 * 3 + 0], aXYZt[i0 * 3 + 1], aXYZt[i0 * 3 + 2]},
        {aXYZt[i1 * 3 + 0], aXYZt[i1 * 3 + 1], aXYZt[i1 * 3 + 2]},
        {aXYZt[i2 * 3 + 0], aXYZt[i2 * 3 + 1], aXYZt[i2 * 3 + 2]}};
    double C[3], dCdp[3][9];
    CdC_StVK(C, dCdp, P, p);
    const double mass[3] = {1, 1, 1};
    const unsigned int aIP[3] = {i0, i1, i2};
    PBD_Update_Const3(
        aXYZt,
        3, 3, mass, C, &dCdp[0][0], aIP, 1.0);
  }
}


} // namespace delfem2


#endif /* DFM2_PBD_STVK_H */
