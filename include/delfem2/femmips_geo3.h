/*
 * Copyright (c) 2019 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */


#ifndef DFM2_FEMMIPS_GEO3_H
#define DFM2_FEMMIPS_GEO3_H

#include <vector>

#include "delfem2/dfm2_inline.h"
#include "delfem2/mat3.h"
#include "delfem2/vec3.h"

namespace delfem2 {

DFM2_INLINE void WdWddW_MIPS(
    double& E, double dE[3][3], double ddE[3][3][3][3],
    const double c[3][3],
    const double C[3][3]);  

}

#ifndef DFM2_STATIC_LIBRARY
#  include "delfem2/femmips_geo3.cpp"
#endif

#endif /* pbd_v23_h */
