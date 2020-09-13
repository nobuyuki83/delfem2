/*
 * Copyright (c) 2019 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */

#ifndef DFM2_IMGIO_H
#define DFM2_IMGIO_H

#include "delfem2/dfm2_inline.h"
#include <vector>
#include <string>

namespace delfem2 {

DFM2_INLINE bool LoadImage_PPMBinary(
    const std::string &filename,
    std::vector<unsigned char> &image,
    int &width, int &height);

/**
 *
 * @param width
 * @param height
 * @param image binary image data. top left corner is the origin. width first
 * @param fname file path
 * @return 0: success, 1: cannot open 2: format different
 */
DFM2_INLINE int LoadImage_PPMAscii(
    unsigned int &width, unsigned int &height,
    std::vector<unsigned char> &image,
    const std::string &fname);


class SFile_TGA {
public:
  unsigned char imageTypeCode;
  short int imageWidth;
  short int imageHeight;
  unsigned char bitCount;
  unsigned char *imageData;
};

DFM2_INLINE bool LoadTGAFile(
    const std::string &filename,
    SFile_TGA *tgaFile);

DFM2_INLINE void ImageInterpolation_Bilinear(
    std::vector<double>& aColor,
    int width,
    int height,
    const unsigned char* img,
    const double* aXY,
    unsigned int nXY);

}

#ifdef DFM2_HEADER_ONLY
#  include "delfem2/imgio.cpp"
#endif

#endif
