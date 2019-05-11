/*
 * Copyright (c) 2019 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */


#include <cmath>

#include "delfem2/quat.h"

#ifndef M_PI
#define M_PI 3.141592
#endif


// multiply two quaternion
void QuatMult(double r[], const double p[], const double q[])
{
  r[0] = p[0] * q[0] - p[1] * q[1] - p[2] * q[2] - p[3] * q[3];
  r[1] = p[0] * q[1] + p[1] * q[0] + p[2] * q[3] - p[3] * q[2];
  r[2] = p[0] * q[2] - p[1] * q[3] + p[2] * q[0] + p[3] * q[1];
  r[3] = p[0] * q[3] + p[1] * q[2] - p[2] * q[1] + p[3] * q[0];
}

void QuatSetIdentity(double q[]){
  q[0] = 1;
  q[1] = 0;
  q[2] = 0;
  q[3] = 0;
}

// transform vector with quaternion
void QuatVec(double vo[], const double q[], const double vi[])
{
  double x2 = q[1] * q[1] * 2.0;
  double y2 = q[2] * q[2] * 2.0;
  double z2 = q[3] * q[3] * 2.0;
  double xy = q[1] * q[2] * 2.0;
  double yz = q[2] * q[3] * 2.0;
  double zx = q[3] * q[1] * 2.0;
  double xw = q[1] * q[0] * 2.0;
  double yw = q[2] * q[0] * 2.0;
  double zw = q[3] * q[0] * 2.0;
  
  vo[0] = (1.0 - y2 - z2)*vi[0] + (xy + zw      )*vi[1] + (zx - yw      )*vi[2];
  vo[1] = (xy - zw      )*vi[0] + (1.0 - z2 - x2)*vi[1] + (yz + xw      )*vi[2];
  vo[2] = (zx + yw      )*vi[0] + (yz - xw      )*vi[1] + (1.0 - x2 - y2)*vi[2];
}

// transform vector with conjugate of quaternion
void QuatConjVec(double vo[], const double q[], const double vi[])
{
  double x2 = q[1] * q[1] * 2.0;
  double y2 = q[2] * q[2] * 2.0;
  double z2 = q[3] * q[3] * 2.0;
  double xy = q[1] * q[2] * 2.0;
  double yz = q[2] * q[3] * 2.0;
  double zx = q[3] * q[1] * 2.0;
  double xw = q[1] * q[0] * 2.0;
  double yw = q[2] * q[0] * 2.0;
  double zw = q[3] * q[0] * 2.0;
  
  vo[0] = (1.0 - y2 - z2)*vi[0] + (xy - zw      )*vi[1] + (zx + yw      )*vi[2];
  vo[1] = (xy + zw      )*vi[0] + (1.0 - z2 - x2)*vi[1] + (yz - xw      )*vi[2];
  vo[2] = (zx - yw      )*vi[0] + (yz + xw      )*vi[1] + (1.0 - x2 - y2)*vi[2];
}

// copy quaternion
void QuatCopy(double r[], const double p[])
{
  r[0] = p[0];
  r[1] = p[1];
  r[2] = p[2];
  r[3] = p[3];
}

void QuatGetAffineMatrix(double r[], const double q[])
{
  double x2 = q[1] * q[1] * 2.0;
  double y2 = q[2] * q[2] * 2.0;
  double z2 = q[3] * q[3] * 2.0;
  double xy = q[1] * q[2] * 2.0;
  double yz = q[2] * q[3] * 2.0;
  double zx = q[3] * q[1] * 2.0;
  double xw = q[1] * q[0] * 2.0;
  double yw = q[2] * q[0] * 2.0;
  double zw = q[3] * q[0] * 2.0;
  
  r[ 0] = 1.0 - y2 - z2;
  r[ 1] = xy - zw;
  r[ 2] = zx + yw;
  r[ 4] = xy + zw;
  r[ 5] = 1.0 - z2 - x2;
  r[ 6] = yz - xw;
  r[ 8] = zx - yw;
  r[ 9] = yz + xw;
  r[10] = 1.0 - x2 - y2;
  r[ 3] = r[ 7] = r[11] = r[12] = r[13] = r[14] = 0.0;
  r[15] = 1.0;
}

void QuatGetAffineMatrixTrans(double r[], const double q[])
{
  double x2 = q[1] * q[1] * 2.0;
  double y2 = q[2] * q[2] * 2.0;
  double z2 = q[3] * q[3] * 2.0;
  double xy = q[1] * q[2] * 2.0;
  double yz = q[2] * q[3] * 2.0;
  double zx = q[3] * q[1] * 2.0;
  double xw = q[1] * q[0] * 2.0;
  double yw = q[2] * q[0] * 2.0;
  double zw = q[3] * q[0] * 2.0;
  
  r[ 0] = 1.0 - y2 - z2;
  r[ 1] = xy + zw;
  r[ 2] = zx - yw;
  r[ 4] = xy - zw;
  r[ 5] = 1.0 - z2 - x2;
  r[ 6] = yz + xw;
  r[ 8] = zx + yw;
  r[ 9] = yz - xw;
  r[10] = 1.0 - x2 - y2;
  r[ 3] = r[ 7] = r[11] = r[12] = r[13] = r[14] = 0.0;
  r[15] = 1.0;
}

void QuatNormalize(double q[])
{
  const double len = sqrt(q[0]*q[0]+q[1]*q[1]+q[2]*q[2]+q[3]*q[3]);
  double invlen = 1.0/len;
  q[0] *= invlen;
  q[1] *= invlen;
  q[2] *= invlen;
  q[3] *= invlen;
}

//////////////////////////////////////////////////////////////////////
// �\�z/����
//////////////////////////////////////////////////////////////////////

/*
CQuaternion::CQuaternion(const CVector3D& axis ){
	AxisToQuat( axis );
}

CQuaternion::CQuaternion(double real, const CVector3D& vector){
	this->real = real;
	this->vector = vector;
}

CQuaternion::CQuaternion(const CVector3D& a_vector, const CVector3D& b_vector){
	VectorTrans(a_vector, b_vector);
}
 */


//////////////////////////////////////////////////////////////////////
// �񃁃��o�֐��̔�t�����h�֐�
//////////////////////////////////////////////////////////////////////

/*
CQuaternion operator+(const CQuaternion& lhs, const CQuaternion& rhs){
    CQuaternion temp = lhs;
	temp += rhs;
    return temp;
}

CQuaternion operator-(const CQuaternion& lhs, const CQuaternion& rhs){
    CQuaternion temp = lhs;
	temp -= rhs;
    return temp;
}

CQuaternion operator*(const CQuaternion& quat, double d){
    CQuaternion temp = quat;
	temp *= d;
    return temp;
}

CQuaternion operator*(double d, const CQuaternion& quat){
    CQuaternion temp = quat;
    temp *= d;
	return temp;
}

CQuaternion operator*(const CQuaternion& lhs, const CQuaternion& rhs){
    CQuaternion temp = lhs;
    temp *= rhs;
	return temp;
}

CVector3D Rotate(const CQuaternion& quat, const CVector3D& vec){
	CQuaternion tmp(0, vec);
	tmp = quat *  (tmp * quat.GetConjugate());  
//	tmp = quat.GetConjugate() *  tmp * quat ;   
	return tmp.GetVector();
}

CVector3D UnRotate(const CQuaternion& quat, const CVector3D& vec){
	CQuaternion tmp(0, vec);
	tmp = (quat.GetConjugate() *  tmp) * quat ; 
//	tmp = quat *  tmp * quat.GetConjugate();    
	return tmp.GetVector();
}


//////////////////////////////////////////////////////////////////////
// �����o�֐��̃t�����h�֐�
//////////////////////////////////////////////////////////////////////


bool operator==(const CQuaternion& lhs, const CQuaternion& rhs){
	if( fabs(lhs.real - rhs.real) < NEARLY_ZERO && lhs == rhs )	return true;
	else return false;
}

bool operator!=(const CQuaternion& lhs, const CQuaternion& rhs){
	if( lhs == rhs ) return false;
	else return true;
}

//////////////////////////////////////////////////////////////////////
// �����o�֐��̔�t�����h�֐��@
//////////////////////////////////////////////////////////////////////

CQuaternion& CQuaternion::operator = (const CQuaternion& rhs){
	if( this != &rhs ){
		real = rhs.real;
		vector = rhs.vector;
	}
	return *this;
}

CQuaternion& CQuaternion::operator *= (const CQuaternion& rhs){
	const double last_real = real;
	real = real * rhs.real - Dot( vector, rhs.vector );
	vector = ( (rhs.vector * last_real) + (vector * rhs.real) ) + Cross( vector, rhs.vector );
	return *this;
}

CQuaternion& CQuaternion::operator *= (double d){
	real *= d;
	vector *= d;
	return *this;
}

CQuaternion& CQuaternion::operator += (const CQuaternion& rhs){
	real += rhs.real;
	vector += rhs.vector;
	return *this;
}

CQuaternion& CQuaternion::operator -= (const CQuaternion& rhs){
	real -= rhs.real;
	vector -= rhs.vector;
	return *this;
}

CQuaternion CQuaternion::GetInverse() const
{
  double invslen = 1.0/SquareLength();
	CQuaternion tmp(real*invslen, vector*(-invslen));
	return tmp;
}

CQuaternion CQuaternion::GetConjugate() const{
	CQuaternion tmp(real, vector*(-1.0));
	return tmp;
}

void CQuaternion::RotMatrix44(double* m) const
{
  
	double vx=vector.x, vy=vector.y, vz=vector.z;
  
	m[ 0] = 1.0 - 2.0 * ( vy * vy + vz * vz );
	m[ 4] = 2.0 * ( vx * vy - vz * real );
	m[ 8] = 2.0 * ( vz * vx + vy * real );
	m[12] = 0;
  
	m[ 1] = 2.0 * ( vx * vy + vz * real );
	m[ 5] = 1.0 - 2.0 * ( vz * vz + vx * vx );
	m[ 9] = 2.0 * ( vy * vz - vx * real );
	m[13] = 0.0;
  
	m[ 2] = 2.0 * ( vz * vx - vy * real );
	m[ 6] = 2.0 * ( vy * vz + vx * real );
	m[10] = 1.0 - 2.0 * ( vy * vy + vx * vx );    
	m[14] = 0.0;  
  
  m[ 3] = 0.0;
  m[ 7] = 0.0;
  m[11] = 0.0;
  m[15] = 1.0;  
}



void CQuaternion::RotMatrix33(double* m) const
{
	double vx=vector.x, vy=vector.y, vz=vector.z;

	m[0] = 1.0 - 2.0 * ( vy * vy + vz * vz );
	m[1] = 2.0 * ( vx * vy - vz * real );
	m[2] = 2.0 * ( vz * vx + vy * real );

	m[3] = 2.0 * ( vx * vy + vz * real );
	m[4] = 1.0 - 2.0 * ( vz * vz + vx * vx );
	m[5] = 2.0 * ( vy * vz - vx * real );

	m[6] = 2.0 * ( vz * vx - vy * real );
	m[7] = 2.0 * ( vy * vz + vx * real );
	m[8] = 1.0 - 2.0 * ( vy * vy + vx * vx );
}

void CQuaternion::Normalize()
{	
	const double len = ( real * real + vector.DLength() );
	real /= len;
	vector /= len;
}

void CQuaternion::AxisToQuat(const CVector3D &axis )
{
	const double phi = axis.Length();
	if( phi < 1.0e-30 ){
		vector = CVector3D(0,0,0);
		real = 1.0;
		return;
	}
	vector = axis;
	vector.SetNormalizedVector();
	vector *= sin( phi * 0.5 );
	real = cos( phi * 0.5 );
}

void CQuaternion::VectorTrans(const CVector3D& a_vector, const CVector3D& b_vector){
	vector = Cross(a_vector, b_vector);
	if( vector.DLength() < NEARLY_ZERO ){
		real = 1.0;
		vector.SetZero();
		return;
	}
	vector.SetNormalizedVector();
	double cos_theta = Dot(a_vector, b_vector) / ( a_vector.Length() * b_vector.Length() );
	real = sqrt( 0.5*(1+cos_theta) );
	vector *= sqrt( 0.5*(1-cos_theta) );
}

double CQuaternion::SquareLength() const
{
	return real*real + vector.DLength();
}


double CQuaternion::Length() const
{
	return sqrt( real*real + vector.DLength() );
}

CQuaternion SphericalLinearInterp(const CQuaternion& q0, const CQuaternion& q1, double t)
{  
  
  double qr = q0.real * q1.real + Dot(q0.vector,q1.vector);
  double ss = 1.0 - qr * qr;
  
  if (ss == 0.0) { return q0; }
  double sp = sqrt(ss);
  double ph = acos(qr);
  double pt = ph * t;
  double t1 = sin(pt) / sp;
  double t0 = sin(ph - pt) / sp;    
  CQuaternion q;
  q.real = t0*q0.real + t1*q1.real;
  q.vector = t0*q0.vector + t1*q1.vector;
  return q;
}

*/
