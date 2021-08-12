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

#ifndef DFM2_LSVECX_H
#define DFM2_LSVECX_H

#include <vector>
#include <cassert>
#include <complex>
#include <iostream>

#include "delfem2/dfm2_inline.h"

namespace delfem2 {

template<typename REAL>
class CVecX {
 public:
  CVecX(REAL *p_, std::size_t n_) : p(p_), n(n_) {}
  CVecX(std::vector<REAL> &v) : p(v.data()), n(v.size()) {}
  //
  REAL dot(const CVecX &rhs) const {
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
  void operator=(const CVecX &rhs) {
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
using CVecXd = CVecX<double>;
using CVecXf = CVecX<float>;
using CVecXcd = CVecX<std::complex<double> >;
using CVecXcf = CVecX<std::complex<float> >;

template<typename REAL>
void AddScaledVec(
    CVecX<REAL> &y,
    double alpha,
    const CVecX<REAL> &x) {
  assert(y.n == x.n);
  const std::size_t n = x.n;
  for (unsigned int i = 0; i < n; i++) {
    y.p[i] += alpha * x.p[i];
  }
}

template<typename REAL>
void ScaleAndAddVec(
    CVecX<REAL> &y,
    REAL beta,
    const CVecX<REAL> &x) {
  assert(y.n == x.n);
  const std::size_t n = x.n;
  for (unsigned int i = 0; i < n; i++) {
    y.p[i] = beta * y.p[i] + x.p[i];
  }
}

template<typename REAL, class MAT>
void AddMatVec(
    CVecX<REAL> &lhs,
    REAL scale_lhs,
    REAL scale_rhs,
    const MAT &mat,
    const CVecX<REAL> &rhs) {
  assert(lhs.n == rhs.n);
  mat.MatVec(lhs.p,
             scale_rhs, rhs.p, scale_lhs);
}

template<typename REAL, class PREC>
void SolvePrecond(
    CVecX<REAL> &Pr_vec,
    const PREC &ilu) {
  ilu.SolvePrecond(Pr_vec.p);
}

template<typename REAL>
REAL Dot(
    const CVecX<REAL> &x,
    const CVecX<REAL> &y) {
  return x.dot(y);
}

} // delfem2


#endif // LSVECX
