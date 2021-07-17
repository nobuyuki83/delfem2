/*
 * Copyright (c) 2019 Nobuyuki Umetani
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

#include <random>
#include "delfem2/dfm2_inline.h"

namespace delfem2 {

template <typename T>
DFM2_INLINE void Normalize_Quat(
    T q[4]);

template <typename T>
DFM2_INLINE void Inverse_Quat(
    T qinv[4],
    T q[4])
{
  const T sqlen = q[0]*q[0]+q[1]*q[1]+q[2]*q[2]+q[3]*q[3];
  const T sqleninv = 1/sqlen;
  qinv[0] = +q[0]*sqleninv;
  qinv[1] = -q[1]*sqleninv;
  qinv[2] = -q[2]*sqleninv;
  qinv[3] = -q[3]*sqleninv;
}


/**
 * @brief Set Identity in the quaternion
 */
template <typename T>
DFM2_INLINE void Quat_Identity(
    T q[4]);

template <typename REAL>
DFM2_INLINE void Quat_Bryant(
    REAL q[4],
    REAL x, REAL y, REAL z);

/**
 * @brief Quaternion for cartesian rotation angle (3D axis with magnitude of rotation angle)
 * @tparam REAL float and double
 * @param q (out)
 * @param a (in)
 */
template <typename REAL>
DFM2_INLINE void Quat_CartesianAngle(
    REAL q[4],
    const REAL a[3]);


/**
 * @func copy quaternion
 * @tparam REAL float or double
 */
template <typename REAL>
DFM2_INLINE void Copy_Quat(
    REAL r[],
    const REAL p[]);

/**
 * @brief multiply two quaternion
 * @tparam REAL float or double
 * @param r (out)
 * @param p (in) lhs quaternion as 4D array (r,x,y,z)
 * @param q (in) rhs quaternion as 4D array (r,z,y,z)
 * @details quaternions don't commute (qp!=pq)
 */
template <typename REAL>
DFM2_INLINE void QuatQuat(
    REAL r[],
    const REAL p[],
    const REAL q[]);

template <typename REAL>
DFM2_INLINE REAL Dot_Quat(
    const REAL p[],
    const REAL q[]);

template <typename T>
DFM2_INLINE T Length_Quat(
    const T q[]);

// -----------------------------
// below: quaternion and vector

/**
 * @func transform a 3D vector with quaternion vo  = q*vi*adj(q)
 * @tparam REAL float or double
 */
template <typename REAL>
DFM2_INLINE void QuatVec(
    REAL vo[],
    const REAL q[],
    const REAL vi[]);

DFM2_INLINE void QuatConjVec(
    double vo[],
    const double q[],
    const double vi[]);

// -------------------------------------------------------

template <typename T>
class CQuat;

template <typename T>
CQuat<T> operator+(const CQuat<T>&, const CQuat<T>&);

template <typename T>
CQuat<T> operator-(const CQuat<T>&, const CQuat<T>&);
  
template <typename T>
CQuat<T> operator*(const CQuat<T>&, T);	//!< multiply scalar
  
template <typename T>
CQuat<T> operator/(const CQuat<T>&, T);	//!< divide by scalar
  
template <typename T>
CQuat<T> operator*(const CQuat<T>&, const CQuat<T>&);

template <typename T>
CQuat<T> SphericalLinearInterp(const CQuat<T>&, const CQuat<T>&, T);


template <typename T>
CQuat<T> operator*(T, const CQuat<T>&);  //!< multiply scalar

  
/**
 * @class class of Quaternion
 **/
template <typename T>
class CQuat  
{
public:
  CQuat() : q{1,0,0,0} {}
  explicit CQuat(const T rhs[4]) : q{rhs[0], rhs[1], rhs[2], rhs[3]} {};
  CQuat(T r, T v0, T v1, T v2) : q{r, v0, v1, v2} {};
  ~CQuat()= default;
  // -----------
  static CQuat Random(T a){
    CQuat<T> q;
    q.q[0] = 1.0;
    q.q[1] = 2*a*rand()/(RAND_MAX+1.0)-a;
    q.q[2] = 2*a*rand()/(RAND_MAX+1.0)-a;
    q.q[3] = 2*a*rand()/(RAND_MAX+1.0)-a;
    Normalize_Quat(q.q);
    return q;
  }
  void CopyTo(float* q1) const {
    q1[0] = q[0];
    q1[1] = q[1];
    q1[2] = q[2];
    q1[3] = q[3];
  }
  void CopyTo(double* q1) const {
    q1[0] = q[0];
    q1[1] = q[1];
    q1[2] = q[2];
    q1[3] = q[3];
  }
  CQuat<double> Double() const {
    return CQuat<double>((double)q[0], (double)q[1], (double)q[2], (double)q[3]);
  }
  CQuat<T> Conjugate() const {
    return CQuat<T>(q[0], -q[1], -q[2], -q[3]);
  }
  
  void SetNormalized();
  void SetSmallerRotation();
  
  static CQuat<T> Identity() {
    return CQuat<T>(1,0,0,0);
  }
public:
  T q[4]; // w,x,y,z
};
using CQuatd = CQuat<double>;
using CQuatf = CQuat<float>;
  
}

#ifdef DFM2_HEADER_ONLY
  #include "delfem2/quat.cpp"
#endif

#endif // !defined(DFM2_QUAT_H)
