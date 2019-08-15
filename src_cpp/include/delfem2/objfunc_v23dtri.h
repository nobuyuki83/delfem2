/*
 * Copyright (c) 2019 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */


#ifndef OBJFUNC_V23DTRI_h
#define OBJFUNC_V23DTRI_h

#include <vector>

#include "delfem2/dyntri.h"
#include "delfem2/dyntri_v2.h"

void PBD_TriStrain(std::vector<double>& aXYZt,
                   const std::vector<ETri>& aETri,
                   const std::vector<CVector2>& aVec2);

void PBD_Bend(std::vector<double>& aXYZt,
              const std::vector<ETri>& aETri,
              const std::vector<CVector2>& aVec2);

#endif /* pbd_v23_h */
