/*
 * Copyright (c) 2020 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */

#ifndef DFEM2_MATN_H
#define DFEM2_MATN_H

/**
 * @file template class/function for vector/matrix where size is determined at compiling
 * @detail this header's extension is .hpp because it is purely template header
 */

namespace delfem2 {

template<typename REAL, unsigned int n>
bool Inverse_Matrix(
    REAL *a) {
  for (unsigned i = 0; i < n; i++) {
    if (fabs(a[i * n + i]) < 1.0e-30) {
      return false;
    }
    {
      const REAL tmp1 = 1.0 / a[i * n + i];
      a[i * n + i] = 1;
      for (unsigned int k = 0; k < n; k++) {
        a[i * n + k] *= tmp1;
      }
    }
    for (unsigned int j = 0; j < n; j++) {
      if (j == i) { continue; }
      const REAL tmp2 = a[j * n + i];
      a[j * n + i] = 0;
      for (unsigned int k = 0; k < n; k++) {
        a[j * n + k] -= tmp2 * a[i * n + k];
      }
    }
  }
}

template<>
bool Inverse_Matrix<double, 3>(double *a) {
  const double det =
      +a[0] * a[4] * a[8] + a[3] * a[7] * a[2] + a[6] * a[1] * a[5]
          - a[0] * a[7] * a[5] - a[6] * a[4] * a[2] - a[3] * a[1] * a[8];
  const double inv_det = 1.0 / det;
  const double t[9] = {a[0], a[1], a[2], a[3], a[4], a[5], a[6], a[7], a[8]};
  a[0] = inv_det * (t[4] * t[8] - t[5] * t[7]);
  a[1] = inv_det * (t[2] * t[7] - t[1] * t[8]);
  a[2] = inv_det * (t[1] * t[5] - t[2] * t[4]);
  a[3] = inv_det * (t[5] * t[6] - t[3] * t[8]);
  a[4] = inv_det * (t[0] * t[8] - t[2] * t[6]);
  a[5] = inv_det * (t[2] * t[3] - t[0] * t[5]);
  a[6] = inv_det * (t[3] * t[7] - t[4] * t[6]);
  a[7] = inv_det * (t[1] * t[6] - t[0] * t[7]);
  a[8] = inv_det * (t[0] * t[4] - t[1] * t[3]);
  return true;
}

/**
 * @param a
 * @return
 * @details partial function template specification is not allowed. so I define the case for "float"
 */
template<>
bool Inverse_Matrix<float, 3>(float *a) {
  const float det =
      +a[0] * a[4] * a[8] + a[3] * a[7] * a[2] + a[6] * a[1] * a[5]
          - a[0] * a[7] * a[5] - a[6] * a[4] * a[2] - a[3] * a[1] * a[8];
  const float inv_det = 1.f / det;
  const float t[9] = {a[0], a[1], a[2], a[3], a[4], a[5], a[6], a[7], a[8]};
  a[0] = inv_det * (t[4] * t[8] - t[5] * t[7]);
  a[1] = inv_det * (t[2] * t[7] - t[1] * t[8]);
  a[2] = inv_det * (t[1] * t[5] - t[2] * t[4]);
  a[3] = inv_det * (t[5] * t[6] - t[3] * t[8]);
  a[4] = inv_det * (t[0] * t[8] - t[2] * t[6]);
  a[5] = inv_det * (t[2] * t[3] - t[0] * t[5]);
  a[6] = inv_det * (t[3] * t[7] - t[4] * t[6]);
  a[7] = inv_det * (t[1] * t[6] - t[0] * t[7]);
  a[8] = inv_det * (t[0] * t[4] - t[1] * t[3]);
  return true;
}

template<>
bool Inverse_Matrix<double, 2>(double *B) {
  const double det = B[0] * B[3] - B[1] * B[2];
  if (fabs(det) < 1.0e-10) return false;
  const double invdet = 1.0 / det;
  const double t[4] = {B[0], B[1], B[2], B[3]};
  B[0] = +invdet * t[3];
  B[1] = -invdet * t[1];
  B[2] = -invdet * t[2];
  B[3] = +invdet * t[0];
  return true;
}

template<>
bool Inverse_Matrix<double, 1>(double *B) {
  B[0] = 1.0 / B[0];
  return true;
}

// ==========================

template<typename REAL, unsigned int N>
void MatMat(
    REAL *C,
    const REAL *A,
    const REAL *B) {
  for (unsigned int i = 0; i < N; ++i) {
    for (unsigned int j = 0; j < N; ++j) {
      C[i * N + j] = 0;
      for (unsigned int k = 0; k < N; ++k) {
        C[i * N + j] += A[i * N + k] * B[k * N + j];
      }
    }
  }
}

template<typename REAL, unsigned int N>
void Sub_MatMat(
    REAL *pVal_ij,
    const REAL *pVal_ik,
    const REAL *pVal_kj) {
  for (unsigned int i = 0; i < N; i++) {
    for (unsigned int j = 0; j < N; ++j) {
      for (unsigned int k = 0; k < N; ++k) {
        pVal_ij[i * N + j] -= pVal_ik[i * N + k] * pVal_kj[k * N + j];
      }
    }
  }
}

// -----------------------------
// below: matrix and vector

template<typename REAL, unsigned int nrow, unsigned int ncol>
void MatVec(
    REAL *y,
    const REAL *A,
    const REAL *x) {
  for (unsigned int i = 0; i < nrow; ++i) {
    y[i] = 0;
    for (unsigned int j = 0; j < ncol; ++j) {
      y[i] += A[i * ncol + j] * x[j];
    }
  }
}

template<typename REAL, unsigned int nrow, unsigned int ncol>
void Sub_MatVec(
    double *pTmpVec,
    const double *pVal_ij,
    const double *valj) {
  for (unsigned int i = 0; i < nrow; ++i) {
    for (unsigned int j = 0; j < ncol; ++j) {
      pTmpVec[i] -= pVal_ij[i * ncol + j] * valj[j];
    }
  }
}

}

#endif /* DFEM2_MATN_H */
