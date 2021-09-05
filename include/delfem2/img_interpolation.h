/*
 * Copyright (c) 2019 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */

#ifndef DFM2_IMG_INTERPOLATION_H
#define DFM2_IMG_INTERPOLATION_H

#include <vector>
#include <string>

#include "delfem2/dfm2_inline.h"

namespace delfem2 {


DFM2_INLINE void ImageInterpolation_Bilinear(
    std::vector<double>& aColor,
    int width,
    int height,
    const unsigned char* img,
    const double* aXY,
    unsigned int nXY);

}  // delfem2

#ifndef DFM2_STATIC_LIBRARY
#  include "delfem2/img_interpolation.cpp"
#endif

#endif
