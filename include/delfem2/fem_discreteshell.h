/*
 * Copyright (c) 2019 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */


#ifndef DFM2_FEM_DISCRETESHELL_H
#define DFM2_FEM_DISCRETESHELL_H

#include <vector>

#include "delfem2/dfm2_inline.h"
#include "delfem2/mat3.h"
#include "delfem2/vec3.h"

namespace delfem2 {

template<typename T>
DFM2_INLINE void CdC_DiscreteShell(
    double &C,
    CVec3<T> dC[4],
    // -----
    const CVec3<T> &p0,
    const CVec3<T> &p1,
    const CVec3<T> &p2,
    const CVec3<T> &p3);
}

#ifndef DFM2_STATIC_LIBRARY
#  include "delfem2/fem_discreteshell.cpp"
#endif

#endif /* DFM2_FEM_DISCRETESHELL_H */
