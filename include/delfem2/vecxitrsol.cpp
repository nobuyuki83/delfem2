/*
 * Copyright (c) 2019 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */

#include "delfem2/vecxitrsol.h"

#include <cassert>
#include <cmath>
#include <vector>
#include <complex>
#include <climits>

// ----------------------------------------------------------------------

template<typename T>
DFM2_INLINE T
delfem2::Dot(
    const std::vector<T> &r_vec,
    const std::vector<T> &u_vec) {
  assert(r_vec.size() == u_vec.size());
  const std::size_t n = r_vec.size();
  T r = 0;
  for (unsigned int i = 0; i < n; i++) { r += r_vec[i] * u_vec[i]; }
  return r;
}
#ifdef DFM2_STATIC_LIBRARY
template float delfem2::Dot(
    const std::vector<float> &r_vec, const std::vector<float> &u_vec);
template double delfem2::Dot(
    const std::vector<double> &r_vec, const std::vector<double> &u_vec);
#endif

namespace delfem2 {

template<>
DFM2_INLINE std::complex<double>
Dot(
    const std::vector<std::complex<double>> &va,
    const std::vector<std::complex<double>> &vb) {
  const std::size_t n = va.size();
  assert(vb.size() == n);
  double sr = 0.0, si = 0.0;
  for (unsigned int i = 0; i < n; ++i) {
    const std::complex<double> &a = va[i];
    const std::complex<double> &b = vb[i];
    sr += a.real() * b.real() + a.imag() * b.imag();
    si += a.imag() * b.real() - a.real() * b.imag();
  }
  return {sr, si};
}

}

// ----------------------------------------------
// dotx

template<typename T>
DFM2_INLINE T delfem2::DotX(
    const T *va,
    const T *vb,
    size_t n) {
  T r = 0.0;
  for (unsigned int i = 0; i < n; i++) { r += va[i] * vb[i]; }
  return r;
}
#ifdef DFM2_STATIC_LIBRARY
template float delfem2::DotX(const float *va, const float *vb, size_t n);
template double delfem2::DotX(const double *va, const double *vb, size_t n);
#endif

namespace delfem2 {

template<>
DFM2_INLINE std::complex<double> DotX(
    const std::complex<double> *va,
    const std::complex<double> *vb,
    size_t n) {
  double sr = 0.0;
  double si = 0.0;
  for (unsigned int i = 0; i < n; i++) {
    const std::complex<double> &a = va[i];
    const std::complex<double> &b = vb[i];
    sr += a.real() * b.real() + a.imag() * b.imag();
    si += a.imag() * b.real() - a.real() * b.imag();
  }
  return {sr, si};
}

}


// ----------------------------------------------
// distance

template<typename T>
T delfem2::Distance(
    const std::vector<T> &va,
    const std::vector<T> &vb) {
  const size_t n = va.size();
  T r = 0.0;
  for (unsigned i = 0; i < n; i++) { r += (va[i] - vb[i]) * (va[i] - vb[i]); }
  return r;
}
#ifdef DFM2_STATIC_LIBRARY
template float delfem2::Distance(
    const std::vector<float> &va, const std::vector<float> &vb);
template double delfem2::Distance(
    const std::vector<double> &va, const std::vector<double> &vb);
#endif



// -------------------------------------------------

namespace delfem2 {

// {y} = {y} + a * {x}
template<typename VAL>
void AXPY(
    VAL a,
    const std::vector<VAL> &x,
    std::vector<VAL> &y) {
  const std::size_t n = x.size();
  assert(y.size() == n);
  for (unsigned int i = 0; i < n; i++) { y[i] += a * x[i]; }
}
#ifdef DFM2_STATIC_LIBRARY
template void AXPY(
    float a,
    const std::vector<float> &x,
    std::vector<float> &y);
template void AXPY(
    double a,
    const std::vector<double> &x,
    std::vector<double> &y);
template void AXPY(
    std::complex<double> a,
    const std::vector<std::complex<double>> &x,
    std::vector<std::complex<double>> &y);
#endif

}

// -----------------------------------------------------------

namespace delfem2 {

// {y} = {y} + a * {x}
template<typename VAL>
void AXPY(
    VAL a,
    const VAL *x,
    VAL *y,
    unsigned int n) {
  for (unsigned int i = 0; i < n; i++) { y[i] += a * x[i]; }
}
#ifdef DFM2_STATIC_LIBRARY
template void AXPY(float a, const float *x, float *y, unsigned int n);
template void AXPY(double a, const double *x, double *y, unsigned int n);
template void AXPY(
    std::complex<double> a,
    const std::complex<double> *x,
    std::complex<double> *y,
    unsigned int n);
#endif

}


// --------

template<typename VAL>
void delfem2::ScaleX(
    VAL *p0,
    VAL s,
    size_t n) {
  for (unsigned int i = 0; i < n; ++i) { p0[i] *= s; }
}
#ifdef DFM2_STATIC_LIBRARY
template void delfem2::ScaleX(float *p0, float s, size_t n);
template void delfem2::ScaleX(double *p0, double s, size_t n);
#endif



// -----------------------------------------------------------



DFM2_INLINE void delfem2::NormalizeX(
    double *p0,
    size_t n) {
  const double ss = delfem2::DotX(p0, p0, n);
  ScaleX(p0, 1.0 / sqrt(ss), n);
}

DFM2_INLINE void delfem2::OrthogonalizeToUnitVectorX(
    double *p1,
    const double *p0,
    size_t n) {
  double d = delfem2::DotX(p0, p1, n);
  for (unsigned int i = 0; i < n; ++i) { p1[i] -= d * p0[i]; }
}

DFM2_INLINE std::complex<double>
delfem2::MultSumX(
    const std::complex<double> *va,
    const std::complex<double> *vb,
    unsigned int n) {
  std::complex<double> s(0, 0);
  for (unsigned int i = 0; i < n; ++i) {
    const std::complex<double> &a = va[i];
    const std::complex<double> &b = vb[i];
    s += a * b;
  }
  return s;
}

// ---------------------------------------------------------------

namespace delfem2 {

template<>
DFM2_INLINE void XPlusAY(
    std::vector<double> &X,
    const size_t nDoF,
    const std::vector<int> &aBCFlag,
    double alpha,
    const std::vector<double> &Y) {
  for (unsigned int i = 0; i < nDoF; ++i) {
    if (aBCFlag[i] != 0) continue;
    X[i] += alpha * Y[i];
  }
}

template<>
DFM2_INLINE void
XPlusAY(
    std::vector<std::complex<double> > &X,
    const size_t nDoF,
    const std::vector<int> &aBCFlag,
    std::complex<double> alpha,
    const std::vector<std::complex<double> > &Y) {
  for (unsigned int i = 0; i < nDoF; ++i) {
    if (aBCFlag[i] != 0) continue;
    X[i] += alpha * Y[i];
  }
}

} // end namespace delfem2

// -------------------------------------

namespace delfem2 {

template<>
DFM2_INLINE void setRHS_Zero(
    std::vector<double> &vec_b,
    const std::vector<int> &aBCFlag,
    int iflag_nonzero) {
  const std::size_t ndof = vec_b.size();
  for (unsigned int i = 0; i < ndof; ++i) {
    if (aBCFlag[i] == iflag_nonzero) continue;
    vec_b[i] = 0;
  }
}

template<>
DFM2_INLINE void setRHS_Zero(
    std::vector<std::complex<double>> &vec_b,
    const std::vector<int> &aBCFlag,
    int iflag_nonzero) {
  const int ndof = (int) vec_b.size();
  for (int i = 0; i < ndof; ++i) {
    if (aBCFlag[i] == iflag_nonzero) continue;
    vec_b[i] = 0;
  }
}

}

DFM2_INLINE void
delfem2::XPlusAYBZ(
    std::vector<double> &X,
    const size_t nDoF,
    const std::vector<int> &aBCFlag,
    double alpha,
    const std::vector<double> &Y,
    double beta,
    const std::vector<double> &Z) {
  for (unsigned int i = 0; i < nDoF; ++i) {
    if (aBCFlag[i] != 0) continue;
    X[i] += alpha * Y[i] + beta * Z[i];
  }
}

DFM2_INLINE void
delfem2::XPlusAYBZCW(
    std::vector<double> &X,
    const size_t nDoF,
    const std::vector<int> &aBCFlag,
    double alpha,
    const std::vector<double> &Y,
    double beta,
    const std::vector<double> &Z,
    double gamma,
    const std::vector<double> &W) {
  for (unsigned int i = 0; i < nDoF; ++i) {
    if (aBCFlag[i] != 0) continue;
    X[i] += alpha * Y[i] + beta * Z[i] + gamma * W[i];
  }
}

DFM2_INLINE void
delfem2::MatVec(
    double *y,
    const double *A,
    unsigned int ncol,
    unsigned int nrow,
    const double *x) {
  for (unsigned int i = 0; i < ncol; ++i) {
    y[i] = 0;
    for (unsigned int j = 0; j < nrow; ++j) {
      y[i] += A[i * nrow + j] * x[j];
    }
  }
}

DFM2_INLINE void
delfem2::MatTVec(
    double *y,
    const double *A,
    unsigned int ncol,
    unsigned int nrow,
    const double *x) {
  for (unsigned int j = 0; j < nrow; ++j) { y[j] = 0; }
  for (unsigned int i = 0; i < ncol; ++i) {
    for (unsigned int j = 0; j < nrow; ++j) {
      y[j] += A[i * nrow + j] * x[i];
    }
  }
}

DFM2_INLINE void
MatMatX(
    double *M,  // [ni, nj]
    unsigned int ni,
    unsigned int nj,
    const double *A,  // [ni, nk]
    unsigned int nk,
    const double *B)  // [nk, nj]
{
  for (unsigned int i = 0; i < ni; ++i) {
    for (unsigned int j = 0; j < nj; ++j) {
      M[i * nj + j] = 0.0;
      for (unsigned int k = 0; k < nk; ++k) {
        M[i * nj + j] += A[i * nk + k] * B[k * nj + j];
      }
    }
  }
}

DFM2_INLINE void MatMatTX(
    double *M,  // [ni, nj]
    unsigned int ni,
    unsigned int nj,
    const double *A,  // [ni, nk]
    unsigned int nk,
    const double *B)  // [nj, nk]
{
  for (unsigned int i = 0; i < ni; ++i) {
    for (unsigned int j = 0; j < nj; ++j) {
      M[i * nj + j] = 0.0;
      for (unsigned int k = 0; k < nk; ++k) {
        M[i * nj + j] += A[i * nk + k] * B[j * nk + k];
      }
    }
  }
}
