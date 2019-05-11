/*
 * Copyright (c) 2019 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */


#ifndef QUAT_H
#define QUAT_H

// stand alone implementation of the quaternion function and class

void QuatSetIdentity(double q[]);
void QuatCopy(double r[], const double p[]);
void QuatMult(double r[], const double p[], const double q[]);
void QuatVec(double vo[], const double q[], const double vi[]);
void QuatConjVec(double vo[], const double q[], const double vi[]);
void QuatGetAffineMatrix(double r[], const double q[]);
void QuatGetAffineMatrixTrans(double r[], const double q[]);
void QuatNormalize(double q[]);

class CQuaternion;

CQuaternion operator+(const CQuaternion&, const CQuaternion&);
CQuaternion operator-(const CQuaternion&, const CQuaternion&);
CQuaternion operator*(double, const CQuaternion&);	//!< multiply scalar
CQuaternion operator*(const CQuaternion&, double);	//!< multiply scalar
CQuaternion operator/(const CQuaternion&, double);	//!< divide by scalar
CQuaternion operator*(const CQuaternion&, const CQuaternion&);
  
CQuaternion SphericalLinearInterp(const CQuaternion&, const CQuaternion&, double);

//! class of Quaternion
class CQuaternion  
{
public:
  CQuaternion(){ q[0]=1.0; q[1]=0.0; q[2]=0.0; q[3]=0.0; }
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
private:
  double q[4];
//	CVector3D vector;	//!< imaginary part
//	double real;	//!< real part
};

#endif // !defined(QUATERNION_H)
