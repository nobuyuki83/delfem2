/*
 * Copyright (c) 2020 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */

#ifndef DFM2_FEM_ROD2_H
#define DFM2_FEM_ROD2_H

#include "delfem2/dfm2_inline.h"

namespace delfem2 {

DFM2_INLINE void WdWddW_Rod2(
    double &W,
    double dW[3][2],
    double ddW[3][3][2][2],
    const double ap[3][2],
    const double aL[2],
    double stiff_stretch01,
    double stiff_stretch12,
    double stiff1_bend);

} // namespace delfem2

#ifndef DFM2_STATIC_LIBRARY
#  include "delfem2/fem_rod2.cpp"
#endif

#endif  /* DFM2_FEM_ROD2_H */
