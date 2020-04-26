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
 * @func transform a 3D vector with quaternion vo  = q*vi*adj(q)
 * @tparam REAL float or double
 */
template <typename REAL>
DFM2_INLINE void QuatVec(
    REAL vo[],
    const REAL q[],
    const REAL vi[]);

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

DFM2_INLINE void QuatConjVec(
    double vo[],
    const double q[],
    const double vi[]);

// ----------


/*
DFM2_INLINE void Mat4_Quat(
    double r[],
    const double q[]);

DFM2_INLINE void Mat4_QuatConj(
    double r[],
    const double q[]);
 */

/**
 * @brief applying transformation in the order of scale, rotation and translation
 */
 /*
void Mat4_ScaleRotTrans(double m[16],
                        double scale, const double quat[4], const double trans[3]);
void MatMat4(double m01[16],
             const double m0[16], const double m1[16]);
void Copy_Mat4(double m1[16],
               const double m0[16]);
               */
  
// -------------------------------------------------------

template <typename T>
class CQuat;

template <typename T>
CQuat<T> operator+(const CQuat<T>&, const CQuat<T>&);

template <typename T>
CQuat<T> operator-(const CQuat<T>&, const CQuat<T>&);

template <typename T>
CQuat<T> operator*(double, const CQuat<T>&);	//!< multiply scalar
  
template <typename T>
CQuat<T> operator*(const CQuat<T>&, T);	//!< multiply scalar
  
template <typename T>
CQuat<T> operator/(const CQuat<T>&, T);	//!< divide by scalar
  
template <typename T>
CQuat<T> operator*(const CQuat<T>&, const CQuat<T>&);

template <typename T>
CQuat<T> SphericalLinearInterp(const CQuat<T>&, const CQuat<T>&, T);

  
/**
 * @class class of Quaternion
 **/
template <typename T>
class CQuat  
{
public:
  CQuat() : q{1,0,0,0} {}
  CQuat(const T rhs[4]) : q{rhs[0], rhs[1], rhs[2], rhs[3]} {};
  CQuat(T r, T v0, T v1, T v2) : q{r, v0, v1, v2} {};
  ~CQuat(){}
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
  /*
	CQuaternion(const CQuaternion& rhs)
		:vector( rhs.vector ){
		real = rhs.real;
	}
	CQuaternion(const CVector3D& axis);
	CQuaternion(double real, const CVector3D& vector);
	CQuaternion(const CVector3D& a_vector, const CVector3D& b_vector);
	~CQuaternion();

	CQuaternion GetConjugate() const;	//!< get conjugate quaternion
  CQuaternion GetInverse() const;
	double GetReal() const{ return real; }	//!< get real part
	CVector3D GetVector(){ return vector; }	//!< get imaginary part

	//! normalization
	void Normalize();
	//! set unit quaternion
	void SetUnit(){ real = 1.0; vector.SetZero(); }
	//! initialize from axial rotation vector
	void AxisToQuat(const CVector3D& axis);
	void VectorTrans(const CVector3D& a_vector, const CVector3D& b_vector);
  void RotMatrix33(double* m) const;
  void RotMatrix44(double* m) const;
  

	friend bool operator==(const CQuaternion&, const CQuaternion&);
	friend bool operator!=(const CQuaternion&, const CQuaternion&);
  friend CQuaternion SphericalLinearInterp(const CQuaternion&, const CQuaternion&, double);

	CQuaternion& operator=(const CQuaternion&);
	CQuaternion& operator+=(const CQuaternion&);
	CQuaternion& operator-=(const CQuaternion&);
	CQuaternion& operator*=(const CQuaternion&);
	CQuaternion& operator*=(double);
	CQuaternion& operator/=(const CQuaternion&);
  
	double Length() const;
  double SquareLength() const;
	*/
public:
  T q[4]; // w,x,y,z
//	CVector3D vector;	//!< imaginary part
//	double real;	//!< real part
};
using CQuatd = CQuat<double>;
using CQuatf = CQuat<float>;
  
}

#ifdef DFM2_HEADER_ONLY
#  include "delfem2/quat.cpp"
#endif

#endif // !defined(DFM2_QUAT_H)
