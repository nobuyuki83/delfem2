/*
 * Copyright (c) 2020 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */

#ifndef DFM2_SAMPLER_H
#define DFM2_SAMPLER_H

#include <vector>
#include <climits>
#include <cassert>

namespace delfem2 {

/**
 * @brief one step of importance sampling using the Lloyd method. Brute force approach.
 * @param aXY
 * @param ndiv grid resolution
 * @param aD importance defined on the grid
 * @param min_aabb
 * @param max_aabb
 */
void Step_Lloyd2(
    std::vector<double> &aXY,
    //
    const unsigned int ndiv,
    const std::vector<double> &aD,
    const double min_aabb[2],
    const double max_aabb[2]) {
  assert(aD.size() == ndiv * ndiv);
  const unsigned int np = static_cast<unsigned int>(aXY.size() / 2);
  std::vector<unsigned int> aV; // Volonoi, index of point, distance for every pixels
  aV.resize(ndiv * ndiv);
  for (unsigned int ih = 0; ih < ndiv; ++ih) {
    for (unsigned int iw = 0; iw < ndiv; ++iw) {
      const double rx0 = (iw + 0.5) / ndiv;
      const double ry0 = 1.0 - (ih + 0.5) / ndiv;
      const double x0 = min_aabb[0] + (max_aabb[0] - min_aabb[0]) * rx0;
      const double y0 = min_aabb[1] + (max_aabb[1] - min_aabb[1]) * ry0;
      double min_dist = 4.0;
      unsigned int min_ip = UINT_MAX;
      for (unsigned int ip = 0; ip < np; ++ip) {
        const double x1 = aXY[ip * 2 + 0];
        const double y1 = aXY[ip * 2 + 1];
        const double d01 = (x0 - x1) * (x0 - x1) + (y0 - y1) * (y0 - y1);
        if (d01 > min_dist) { continue; }
        min_dist = d01;
        min_ip = ip;
      }
      aV[ih * ndiv + iw] = min_ip;
    }
  }

  std::vector<double> awpw(np * 3, 1.0e-10);
  for (unsigned int ih = 0; ih < ndiv; ++ih) {
    for (unsigned int iw = 0; iw < ndiv; ++iw) {
      const double rx0 = (iw + 0.5) / ndiv;
      const double ry0 = 1.0 - (ih + 0.5) / ndiv;
      const double x0 = min_aabb[0] + (max_aabb[0] - min_aabb[0]) * rx0;
      const double y0 = min_aabb[1] + (max_aabb[1] - min_aabb[1]) * ry0;
      double w0 = aD[ih * ndiv + iw];
      const unsigned int ip0 = aV[ih * ndiv + iw];
      awpw[ip0 * 3 + 0] += w0 * x0;
      awpw[ip0 * 3 + 1] += w0 * y0;
      awpw[ip0 * 3 + 2] += w0;
    }
  }

  for (unsigned int ip = 0; ip < np; ++ip) {
    aXY[ip * 2 + 0] = awpw[ip * 3 + 0] / awpw[ip * 3 + 2];
    aXY[ip * 2 + 1] = awpw[ip * 3 + 1] / awpw[ip * 3 + 2];
  }
}

}

#endif // SAMPLER