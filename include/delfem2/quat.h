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

namespace delfem2 {

template <typename T>
void Normalize_Quat(T q[4]);

template <typename T>
void SetIdentity_Quat(T q[4]);
  
  
void QuatCopy(double r[], const double p[]);
void QuatQuat(double r[], const double p[], const double q[]);
void QuatVec(double vo[], const double q[], const double vi[]);
void QuatConjVec(double vo[], const double q[], const double vi[]);
void Quat_Bryant(double q[4],
                 double x, double y, double z);
void Mat4_Quat(double r[],
               const double q[]);
void Mat4_QuatConj(double r[],
                   const double q[]);
void Mat4_ScaleRotTrans(double m[16],
                        double scale, const double quat[4], const double trans[3]);
void MatMat4(double m01[16],
             const double m0[16], const double m1[16]);
void Copy_Mat4(double m1[16],
               const double m0[16]);
  
// -------------------------------------------------------

template <typename T>
class CQuaternion;

template <typename T>
CQuaternion<T> operator+(const CQuaternion<T>&, const CQuaternion<T>&);

template <typename T>
CQuaternion<T> operator-(const CQuaternion<T>&, const CQuaternion<T>&);

template <typename T>
CQuaternion<T> operator*(double, const CQuaternion<T>&);	//!< multiply scalar
  
template <typename T>
CQuaternion<T> operator*(const CQuaternion<T>&, T);	//!< multiply scalar
  
template <typename T>
CQuaternion<T> operator/(const CQuaternion<T>&, T);	//!< divide by scalar
  
template <typename T>
CQuaternion<T> operator*(const CQuaternion<T>&, const CQuaternion<T>&);

template <typename T>
CQuaternion<T> SphericalLinearInterp(const CQuaternion<T>&, const CQuaternion<T>&, T);

  
/**
 * @brief class of Quaternion
 **/
template <typename T>
class CQuaternion  
{
public:
  CQuaternion() : q{1,0,0,0} {}
  CQuaternion(const T rhs[4]) : q{rhs[0], rhs[1], rhs[2], rhs[3]} {};
  CQuaternion(T r, T v0, T v1, T v2) : q{r, v0, v1, v2} {};
  ~CQuaternion(){}
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
  
  
}

#endif // !defined(QUATERNION_H)
