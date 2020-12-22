/*
 * Copyright (c) 2019 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */

/**
 * @file 4x4 matrix class (CMat4) and functions
 */

#ifndef DFM2_GRID2_H
#define DFM2_GRID2_H

#include "delfem2/dfm2_inline.h"
#include <vector>
#include <cassert>
#include <cmath>
#include <iostream>
#include <limits> // using NaN Check


namespace delfem2 {

void LaplacianGrid(
    std::vector<float> &aZLap,
    const std::vector<float> &aZ,
    unsigned int ny,
    unsigned int nx) {
  for (unsigned int iy1 = 0; iy1 < ny; ++iy1) {
    for (unsigned int ix1 = 0; ix1 < nx; ++ix1) {
      unsigned int iy0 = (iy1 == 0) ? 0 : iy1 - 1;
      unsigned int iy2 = (iy1 == ny - 1) ? ny - 1 : iy1 + 1;
      unsigned int ix0 = (ix1 == 0) ? 0 : ix1 - 1;
      unsigned int ix2 = (ix1 == nx - 1) ? nx - 1 : ix1 + 1;
      float v11 = aZ[iy1 * nx + ix1];
      float v01 = aZ[iy1 * nx + ix0];
      float v10 = aZ[iy0 * nx + ix1];
      float v21 = aZ[iy1 * nx + ix2];
      float v12 = aZ[iy2 * nx + ix1];
      aZLap[iy1 * nx + ix1] = 4 * v11 - v01 - v10 - v21 - v12;
    }
  }
}

void LaplacianSmoothGrid(
    std::vector<float> &aZLap,
    double alpha,
    const std::vector<float> &aZ,
    unsigned int ny,
    unsigned int nx) {
  for (unsigned int iy1 = 0; iy1 < ny; ++iy1) {
    for (unsigned int ix1 = 0; ix1 < nx; ++ix1) {
      unsigned int iy0 = (iy1 == 0) ? 0 : iy1 - 1;
      unsigned int iy2 = (iy1 == ny - 1) ? ny - 1 : iy1 + 1;
      unsigned int ix0 = (ix1 == 0) ? 0 : ix1 - 1;
      unsigned int ix2 = (ix1 == nx - 1) ? nx - 1 : ix1 + 1;
      float v11 = aZ[iy1 * nx + ix1];
      float v01 = aZ[iy1 * nx + ix0];
      float v10 = aZ[iy0 * nx + ix1];
      float v21 = aZ[iy1 * nx + ix2];
      float v12 = aZ[iy2 * nx + ix1];
      float avg = 0.25f * (v01 + v10 + v21 + v12);
      aZLap[iy1 * nx + ix1] = alpha * avg + (1 - alpha) * v11;
    }
  }
}

}


#endif