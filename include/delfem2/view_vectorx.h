/*
 * Copyright (c) 2019 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */

/**
 * @file vector class for making array have same interface as the Eigen::VectorX
 *
 * DOONE(Dec. 25th 2020): splitting this file into "vecx.h" and "itersol.h" in the future
 */

#ifndef DFM2_VIEW_VECTORX_H
#define DFM2_VIEW_VECTORX_H

#include <vector>
#include <cassert>
#include <complex>
#include <iostream>

#include "delfem2/dfm2_inline.h"

namespace delfem2 {

template<typename REAL>
class ViewAsVectorX {
 public:
  ViewAsVectorX(REAL *p_, std::size_t n_) : p(p_), n(n_) {}
  ViewAsVectorX(std::vector<REAL> &v) : p(v.data()), n(v.size()) {}
  //
  REAL dot(const ViewAsVectorX &rhs) const {
    assert(n == rhs.n);
    REAL d = 0;
    for (unsigned int i = 0; i < n; ++i) {
      d += p[i] * rhs.p[i];
    }
    return d;
  }
  /**
   * deep copy (copying values not the pointer)
   * @param rhs
   */
  void operator=(const ViewAsVectorX &rhs) {
    assert(n == rhs.n);
    for (unsigned int i = 0; i < n; ++i) { p[i] = rhs.p[i]; }
  }
  //
  void setZero() {
    for (unsigned int i = 0; i < n; ++i) { p[i] = 0; }
  }
 public:
  REAL *const p;
  const std::size_t n;
};
using ViewAsVectorXd = ViewAsVectorX<double>;
using ViewAsVectorXf = ViewAsVectorX<float>;
using ViewAsVectorXcd = ViewAsVectorX<std::complex<double> >;
using ViewAsVectorXcf = ViewAsVectorX<std::complex<float> >;

template<typename REAL>
void AddScaledVec(
    ViewAsVectorX<REAL> &y,
    double alpha,
    const ViewAsVectorX<REAL> &x) {
  assert(y.n == x.n);
  const std::size_t n = x.n;
  for (unsigned int i = 0; i < n; i++) {
    y.p[i] += alpha * x.p[i];
  }
}

template<typename REAL>
void ScaleAndAddVec(
    ViewAsVectorX<REAL> &y,
    REAL beta,
    const ViewAsVectorX<REAL> &x) {
  assert(y.n == x.n);
  const std::size_t n = x.n;
  for (unsigned int i = 0; i < n; i++) {
    y.p[i] = beta * y.p[i] + x.p[i];
  }
}

template<typename REAL, class MAT>
void AddMatVec(
    ViewAsVectorX<REAL> &lhs,
    REAL scale_lhs,
    REAL scale_rhs,
    const MAT &mat,
    const ViewAsVectorX<REAL> &rhs) {
  assert(lhs.n == rhs.n);
  mat.MatVec(lhs.p,
             scale_rhs, rhs.p, scale_lhs);
}

template<typename REAL, class PREC>
void SolvePrecond(
    ViewAsVectorX<REAL> &Pr_vec,
    const PREC &ilu) {
  ilu.SolvePrecond(Pr_vec.p);
}

template<typename REAL>
REAL Dot(
    const ViewAsVectorX<REAL> &x,
    const ViewAsVectorX<REAL> &y) {
  return x.dot(y);
}

} // delfem2


#endif // DFM2_LSVECX_H
