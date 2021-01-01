 /*
 * Copyright (c) 2019 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */


#ifndef DFM2_GEOCONVHULL3_H
#define DFM2_GEOCONVHULL3_H

#include "delfem2/vec3.h"
#include "delfem2/dfm2_inline.h"
#include <cassert>
#include <cmath>
#include <iostream>
#include <vector>

#define NEARLY_ZERO 1.e-16

namespace delfem2 {

template <typename T>
void ConvexHull(
    std::vector<int>& aTri,
    const std::vector<CVec3<T> >& aXYZ);

} // end namespace delfem2

#ifdef DFM2_HEADER_ONLY
#  include "delfem2/geoconvhull3_v3.cpp"
#endif


#endif // VEC3_H
