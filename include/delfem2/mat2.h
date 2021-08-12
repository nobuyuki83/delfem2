/*
 * Copyright (c) 2019 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */

#ifndef DFM2_MAT2_H
#define DFM2_MAT2_H

#include <cstdlib>
#include <cassert>
#include <iostream>
#include <cmath>
#include <vector>

#include "delfem2/dfm2_inline.h"

// -----------------------------------------------------

namespace delfem2 {

template<typename T>
void MatVec2(
    T w[2],
    const T A[4],
    const T v[2]);

template<typename T>
void MatMat2(
    T AB[4],
    const T A[4],
    const T B[4]);

DFM2_INLINE bool InverseMat2(
    double invB[4],
    const double B[4]);

DFM2_INLINE void gramian2(
    double AtA[3],
    const double A[4]);

DFM2_INLINE void VLVt2(
    double A[4],
    double l0,
    double l1,
    const double V[4]);

DFM2_INLINE void RotationalComponentOfMatrix2(
    double R[4],
    const double M[4]);

// ----------------------

template<typename T>
class CMat2; // this pre-definition is needed for following functions

template<typename T>
CMat2<T> operator*(
    T d,
    const CMat2<T> &rhs);

/**
 * @brief 2 dimensional vector class
 */
template<typename T>
class CMat2 {
 public:

  CMat2() : p{0, 0, 0, 0} {}
  CMat2(const CMat2 &rhs) {
    this->p[0] = rhs.p[0];
    this->p[1] = rhs.p[1];
    this->p[2] = rhs.p[2];
    this->p[3] = rhs.p[3];
  }
  CMat2(T xx, T xy, T yx, T yy) {
    this->p[0] = xx;
    this->p[1] = xy;
    this->p[2] = yx;
    this->p[3] = yy;
  }
  CMat2(T val_dia) : p{val_dia, 0, 0, val_dia} {}

  // above: constructor / destructor
  // -------------------------------
  // below: operator
  CMat2 operator-() const {
    return CMat2(-p[0], -p[1], -p[2], -p[3]);
  }

  inline CMat2 &operator+=(const CMat2 &rhs) {
    p[0] += rhs.p[0];
    p[1] += rhs.p[1];
    p[2] += rhs.p[2];
    p[3] += rhs.p[3];
    return *this;
  }

  inline CMat2 &operator-=(const CMat2 &rhs) {
    p[0] -= rhs.p[0];
    p[1] -= rhs.p[1];
    p[2] -= rhs.p[2];
    p[3] -= rhs.p[3];
    return *this;
  }

  inline CMat2 &operator*=(T scale) {
    p[0] *= scale;
    p[1] *= scale;
    p[2] *= scale;
    p[3] *= scale;
    return *this;
  }

  inline CMat2 &operator/=(T d) {
    if (fabs(d) < 1.0e-6) {
      assert(0);
      return *this;
    }
    p[0] /= d;
    p[1] /= d;
    p[2] /= d;
    p[3] /= d;
    return *this;
  }

  inline CMat2 operator+(const CMat2 &rhs) const {
    CMat2 v = *this;
    v += rhs;
    return v;
  }

  inline CMat2 operator-(const CMat2 &rhs) const {
    CMat2 v = *this;
    v -= rhs;
    return v;
  }

  inline T operator[](int i) const {
    assert(i < 4);
    return p[i];
  }

  inline T &operator[](int i) {
    assert(i < 4);
    return p[i];
  }

  inline T operator()(int i, int j) const {
    assert(i < 2 && j < 2);
    return p[i * 2 + j];
  }

  inline T &operator()(int i, int j) {
    assert(i < 2 && j < 2);
    return p[i * 2 + j];
  }

  // above: operator
  // ---------------
  // below: function

  //! @brief set all value to zero
  //! @details naming inspired by Eigen
  inline void setZero() {
    p[0] = T(0);
    p[1] = T(0);
    p[2] = T(0);
    p[3] = T(0);
  }

  //! @brief transpose matrix
  //! @details naming inspired by Eigen
  CMat2 transpose() const {
    return CMat2(p[0], p[2], p[1], p[3]);
  }

  void transposeInPlace() {
    const T t1 = p[1];
    p[1] = p[2];
    p[2] = t1;
  }

  //! @details naming inspired by Eigen
  T determinant() const {
    return p[0] * p[3] - p[1] * p[2];
  }

  template<typename S>
  CMat2<S> cast() const {
    return CMat2<S>(
        static_cast<S>(p[0]),
        static_cast<S>(p[1]),
        static_cast<S>(p[2]),
        static_cast<S>(p[3]));
  };

  void multiply_vec2(T vo[2], const T vi[2]) const {
    MatVec2(
        vo,
        p, vi);
  }

  static CMat2<T> outer_product(const T v0[2], const T v1[2]) {
    return CMat2<T>(v0[0] * v1[0], v0[0] * v1[1], v0[1] * v1[0], v0[1] * v1[1]);
  }

 public:
  T p[4];
};

using CMat2d = CMat2<double>;
using CMat2f = CMat2<float>;
using CMat2i = CMat2<int>;

// --------------------

template<typename T>
void polar_decomposition(
    CMat2<T> &R,
    CMat2<T> &S,
    const CMat2<T> &m);

template<typename T>
void svd(CMat2<T> &U,
         CMat2<T> &sig,
         CMat2<T> &V,
         const CMat2<T> &m);

} // delfem2

#ifndef DFM2_STATIC_LIBRARY
#  include "delfem2/mat2.cpp"
#endif

#endif // VEC_2


