/*
 * Copyright (c) 2019 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */


#ifndef DFM2_DEFARAPENERGY_GEO3_H
#define DFM2_DEFARAPENERGY_GEO3_H

#include <vector>

#include "delfem2/dfm2_inline.h"
#include "delfem2/mat3.h"
#include "delfem2/vec3.h"

namespace delfem2 {
  
DFM2_INLINE double W_ArapEnergy(
    const std::vector<double>& aXYZ0,
    const std::vector<double>& aXYZ1,
    const std::vector<double>& aQuat1,
    const std::vector<unsigned int>& psup_ind,
    const std::vector<unsigned int>& psup);

DFM2_INLINE void dW_ArapEnergy(
    std::vector<double>& aRes,
    const std::vector<double>& aXYZ0,
    const std::vector<double>& aXYZ1,
    const std::vector<double>& aQuat1,
    const std::vector<unsigned int>& psup_ind,
    const std::vector<unsigned int>& psup);

DFM2_INLINE void ddW_ArapEnergy(
    std::vector<double>& eM,
    const std::vector<unsigned int>& aIP,
    const std::vector<double>& aXYZ0,
    const std::vector<double>& aQuat1);

}

#ifndef DFM2_STATIC_LIBRARY
#  include "delfem2/defarapenergy_geo3.cpp"
#endif

#endif /* pbd_v23_h */
