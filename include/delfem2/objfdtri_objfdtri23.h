/*
 * Copyright (c) 2019 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */


#ifndef DFM2_OBJFDTRI_V23DTRI_H
#define DFM2_OBJFDTRI_V23DTRI_H

#include "delfem2/dfm2_inline.h"
#include "delfem2/dtri2_v2dtri.h"

// ------------------------------

namespace delfem2 {

/**
 *
 * @param aXYZt
 * @param nXYZ
 * @param[in] aETri triangle index
 * @param[in] aVec2 2D positions
 */
void PBD_TriStrain(
    double* aXYZt,
    unsigned int nXYZ,
    const std::vector<delfem2::CDynTri>& aETri,
    const std::vector<CVec2d>& aVec2);

/**
 *
 * @param aXYZt
 * @param nXYZ
 * @param aETri
 * @param aVec2
 * @param ratio
 */
void PBD_Bend(
    double* aXYZt,
    unsigned int nXYZ,
    const std::vector<delfem2::CDynTri>& aETri,
    const std::vector<CVec2d>& aVec2,
    double ratio);

} // namespace delfem2

#ifdef DFM2_HEADER_ONLY
#  include "delfem2/objfdtri_objfdtri23.cpp"
#endif

#endif /* pbd_v23_h */
