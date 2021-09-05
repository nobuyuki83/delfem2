/*
 * Copyright (c) 2019 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */

#include <cmath>
#include <cassert>

#include "delfem2/img_interpolation.h"

void delfem2::ImageInterpolation_Bilinear(
    double color[3],
    int width,
    int height,
    const unsigned char *img,
    double x,
    double y) {
  int ix0 = static_cast<int>(floor(x));
  int iy0 = static_cast<int>(floor(y));
  int ix1 = ix0 + 1;
  if (ix1 == width) { ix1 = width - 1; }
  int iy1 = iy0 + 1;
  if (iy1 == height) { iy1 = height - 1; }
  assert(ix0 >= 0 && ix0 < width);
  assert(iy0 >= 0 && iy0 < height);
  assert(ix1 >= 0 && ix1 < width);
  assert(ix1 >= 0 && ix1 < width);
  double rx = x - ix0;
  double ry = y - iy0;
  double w00 = 1.0 / 255.0 * (1 - rx) * (1 - ry);
  double w01 = 1.0 / 255.0 * (1 - rx) * ry;
  double w10 = 1.0 / 255.0 * rx * (1 - ry);
  double w11 = 1.0 / 255.0 * rx * ry;
  for (int i = 0; i < 3; ++i) {
    color[i] = +
        w00 * img[(ix0 + iy0 * width) * 3 + i] +
        w01 * img[(ix0 + iy1 * width) * 3 + i] +
        w10 * img[(ix1 + iy0 * width) * 3 + i] +
        w11 * img[(ix1 + iy1 * width) * 3 + i];
  }
}