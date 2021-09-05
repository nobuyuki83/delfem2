/*
 * Copyright (c) 2019 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */

#ifndef DFM2_IMG_IOPPM_H
#define DFM2_IMG_IOPPM_H

#include <vector>
#include <string>

#include "delfem2/dfm2_inline.h"

namespace delfem2 {

DFM2_INLINE bool LoadImage_PPMBinary(
    const std::string &filename,
    std::vector<unsigned char> &image,
    int &width, int &height);

/**
 *
 * @param[out] width
 * @param[out] height
 * @param[out] image binary image data. top left corner is the origin. width first
 * @param[in] fname file path
 * @return 0: success, 1: cannot open 2: format different
 */
DFM2_INLINE int LoadImage_PPMAscii(
    unsigned int &width, unsigned int &height,
    std::vector<unsigned char> &image,
    const std::string &fname);

}  // namespace delfem2

#ifndef DFM2_STATIC_LIBRARY
#  include "delfem2/img_ioppm.cpp"
#endif

#endif
