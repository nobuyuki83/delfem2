/*
 * Copyright (c) 2019 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */


#ifndef OBJFUNC_V23DTRI_h
#define OBJFUNC_V23DTRI_h

#include <vector>

#include "delfem2/dtri.h"
#include "delfem2/dtri_v2.h"

void PBD_TriStrain(double* aXYZt,
                   unsigned int nXYZ,
                   const std::vector<delfem2::CDynTri>& aETri,
                   const std::vector<CVector2>& aVec2);

void PBD_Bend(double* aXYZt,
              unsigned int nXYZ,
              const std::vector<delfem2::CDynTri>& aETri,
              const std::vector<CVector2>& aVec2);

#endif /* pbd_v23_h */
