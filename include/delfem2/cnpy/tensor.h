/*
 * Copyright (c) 2019 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */

/**
 * @file 4x4 matrix class (CMat4) and functions
 */

#ifndef DFM2_CNPY_TENSOR_H
#define DFM2_CNPY_TENSOR_H

#include "delfem2/dfm2_inline.h"
#include "cnpy/cnpy.h"
#include <vector>
#include <cassert>
#include <cmath>
#include <iostream>
#include <limits> // using NaN Check

// -----------------------------

namespace delfem2 {

template <typename T>
void StableSoftmax(
    double* pout,
    unsigned int nc,
    const T* pin)
{
  assert(nc>0);
  double pmax = pin[0];
  for(unsigned int ic=1;ic<nc;++ic) {
    double p0 = pin[ic];
    if( p0 > pmax ){ pmax = p0; }
  }
  double sum_exp = 0.0;
  for (unsigned int ic = 0; ic < nc; ++ic) {
    sum_exp += exp(pin[ic]-pmax);
  }
  for (unsigned int ic = 0; ic < nc; ++ic) {
    double v0 = exp(pin[ic]-pmax);
    double v1 = v0 / sum_exp;
    pout[ic] = v1;
  }
}


template<typename VAL>
void LoadNumpyArray_BCWH(
    unsigned int &nh,
    unsigned int &nw,
    unsigned int &nc,
    std::vector<VAL> &aFlagImg,
    //
    const std::string &str_npz) {
  cnpy::npz_t npz0 = cnpy::npz_load(str_npz);
  assert(npz0.size() == 1);
  cnpy::NpyArray np0 = npz0.begin()->second;
  assert(np0.shape[0] == 1);
  nc = np0.shape[1];
  nw = np0.shape[2];
  nh = np0.shape[3];
  assert(np0.num_vals == nw * nh * nc);
  assert(np0.num_bytes() / np0.num_vals == sizeof(float));
  const VAL *p = np0.data<VAL>();
  aFlagImg.resize(nh * nw * nc);
  for (unsigned int ih = 0; ih < nh; ++ih) {
    for (unsigned int iw = 0; iw < nw; ++iw) {
      for (unsigned int ic = 0; ic < nc; ++ic) {
        aFlagImg[(ih * nw + iw) * nc + ic] = p[ic * nh * nw + (nh - ih - 1) * nw + iw];
      }
    }
  }
}


} // delfem2


#endif