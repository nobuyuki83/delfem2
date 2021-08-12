/*
 * Copyright (c) 2019-2021 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */

/**
 * @file stand alone implementation of the quaternion function and class
 * @details the  order of the parameters on memory is (w,x,y,z)
 */

#ifndef DFM2_QUAT_H
#define DFM2_QUAT_H

#if defined(_MSC_VER)
#pragma warning( push )
#pragma warning( disable : 4201 )
#endif

#include <random>

#include "delfem2/dfm2_inline.h"

namespace delfem2 {

template<typename T>
DFM2_INLINE void Normalize_Quat(
    T q[4]);

/**
 *
 * @tparam T
 * @param[out] qinv (x,y,z,w)
 * @param[in] q quaternion (x,y,z,w)
 */
template<typename T>
DFM2_INLINE void Inverse_Quat(
    T qinv[4],
    T q[4]) {
  const T sqlen = q[0] * q[0] + q[1] * q[1] + q[2] * q[2] + q[3] * q[3];
  const T sqleninv = 1 / sqlen;
  qinv[0] = -q[0] * sqleninv;
  qinv[1] = -q[1] * sqleninv;
  qinv[2] = -q[2] * sqleninv;
  qinv[3] = +q[3] * sqleninv;
}

/**
 * @brief Set Identity in the quaternion
 * @param[out] q (x,y,z,w)
 */
template<typename T>
DFM2_INLINE void Quat_Identity(
    T q[4]);

/**
 * rotation arounc axis (order x->y->z)
 * @tparam REAL
 * @param[out] q (x,y,z,w)
 * @param[in] x radian
 * @param[in] y radian
 * @param[in] z radian
 */
template<typename REAL>
DFM2_INLINE void Quat_Bryant(
    REAL q[4],
    REAL x, REAL y, REAL z);

/**
 * Quaternion for cartesian rotation angle
 * (3D axis with magnitude of rotation angle)
 * @tparam REAL float and double
 * @param[out] q (x,y,z,w)
 * @param[in] a 3D vector (x,y,z)
 */
template<typename REAL>
DFM2_INLINE void Quat_CartesianAngle(
    REAL q[4],
    const REAL a[3]);

/**
 * copy quaternion
 * @tparam REAL float or double
 */
template<typename REAL>
DFM2_INLINE void Copy_Quat(
    REAL r[],
    const REAL p[]);

/**
 * multiply two quaternion
 * @tparam REAL float or double
 * @param r (out)
 * @param p (in) lhs quaternion as 4D array (x,y,z,w)
 * @param q (in) rhs quaternion as 4D array (x,y,z,w)
 * @details quaternions don't commute (qp!=pq)
 */
template<typename REAL>
DFM2_INLINE void QuatQuat(
    REAL r[],
    const REAL p[],
    const REAL q[]);

/**
 * inner product of two quaternion
 * @tparam REAL
 * @param[in] p (x,y,z,w)
 * @param[in] q (x,y,z,w)
 * @return
 */
template<typename REAL>
DFM2_INLINE REAL Dot_Quat(
    const REAL p[],
    const REAL q[]);

template<typename T>
DFM2_INLINE T Length_Quat(
    const T q[]);

// -----------------------------
// below: quaternion and vector

/**
 * transform a 3D vector with quaternion vo  = q*vi*adj(q)
 * @tparam REAL float or double
 */
template<typename REAL>
DFM2_INLINE void QuatVec(
    REAL vo[],
    const REAL q[],
    const REAL vi[]);

/**
 * trasnsform a 3D vector with conjusgate of a quaternion
 * @param[out] vo
 * @param[in] q quaternion
 * @param[in] vi vector
 */
template<typename T>
DFM2_INLINE void QuatConjVec(
    T vo[],
    const T q[],
    const T vi[]);

// -------------------------------------------------------

template<typename T>
class CQuat;

template<typename T>
CQuat<T> operator+(const CQuat<T> &, const CQuat<T> &);

template<typename T>
CQuat<T> operator-(const CQuat<T> &, const CQuat<T> &);

template<typename T>
CQuat<T> operator*(const CQuat<T> &, T);    //!< multiply scalar

template<typename T>
CQuat<T> operator/(const CQuat<T> &, T);    //!< divide by scalar

template<typename T>
CQuat<T> operator*(const CQuat<T> &, const CQuat<T> &);

template<typename T>
std::ostream &operator<<(std::ostream &output, const CQuat<T> &q);

template<typename T>
CQuat<T> SphericalLinearInterp(const CQuat<T> &, const CQuat<T> &, T);

template<typename T>
CQuat<T> operator*(T, const CQuat<T> &);  //!< multiply scalar


/**
 * @class class of Quaternion
 **/
template<typename T>
class CQuat {
 public:
  CQuat() : p{0, 0, 0, 1} {}

  explicit CQuat(const T rhs[4]) : p{rhs[0], rhs[1], rhs[2], rhs[3]} {};

  /**
   * the order of the argument is (r,x,y,z)
   * but internal storage order is (x,y,z,r)
   */
  CQuat(T w, T x, T y, T z) : p{x, y, z, w} {};

  ~CQuat() = default;
  // -----------
  static CQuat Random(T a) {
    CQuat<T> q;
    q.p[0] = 2 * a * rand() / (RAND_MAX + 1.0) - a;
    q.p[1] = 2 * a * rand() / (RAND_MAX + 1.0) - a;
    q.p[2] = 2 * a * rand() / (RAND_MAX + 1.0) - a;
    q.p[3] = 1.0;
    Normalize_Quat(q.p);
    return q;
  }

  /**
   * raw pointer copy (order of x,y,z,w)
   */
  template<typename T1>
  void CopyTo(T1 *q1) const {
    q1[0] = static_cast<T1>(p[0]);
    q1[1] = static_cast<T1>(p[1]);
    q1[2] = static_cast<T1>(p[2]);
    q1[3] = static_cast<T1>(p[3]);
  }
  template<typename T1>
  CQuat<T1> cast() const {
    // initialization is in the order of (w,x,y,z)
    return CQuat<T1>(static_cast<T1>(w),
                     static_cast<T1>(x),
                     static_cast<T1>(y),
                     static_cast<T1>(z));
  }
  CQuat<T> conjugate() const {
    // initialization is in the order of (w,x,y,z)
    return CQuat<T>(+w, -x, -y, -z);
  }

  T norm() const {
    return std::sqrt(x * x + y * y + z * z + w * w);
  }

  void normalize();

  void SetSmallerRotation();

  static CQuat<T> Identity() {
    return CQuat<T>(1, 0, 0, 0); // initialization is (w,x,y,z)
  }
 public:
  union {
    T p[4];
    struct {
      T x, y, z, w;
    };
  };
};
using CQuatd = CQuat<double>;
using CQuatf = CQuat<float>;

}

#if defined(_MSC_VER)
#pragma warning( pop )
#endif

#ifndef DFM2_STATIC_LIBRARY
#include "delfem2/quat.cpp"
#endif

#endif // !defined(DFM2_QUAT_H)
