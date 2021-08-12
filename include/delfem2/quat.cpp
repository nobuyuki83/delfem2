/*
 * Copyright (c) 2019 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */


#include <cmath>
#include <cstring>

#include "delfem2/quat.h"

#ifndef M_PI
#define M_PI 3.141592
#endif

// ----------------------------------

template<typename REAL>
DFM2_INLINE REAL delfem2::Dot_Quat(
    const REAL p[],
    const REAL q[]) {
  return p[0] * q[0] + p[1] * q[1] + p[2] * q[2] + p[3] * q[3];
}

// ---------------------

template<typename T>
DFM2_INLINE T delfem2::Length_Quat(const T q[]) {
  return std::sqrt(q[0] * q[0] + q[1] * q[1] + q[2] * q[2] + q[3] * q[3]);
}
#ifdef DFM2_STATIC_LIBRARY
template float delfem2::Length_Quat(const float q[]);
template double delfem2::Length_Quat(const double q[]);
#endif

// -----------------------------


template<typename T>
DFM2_INLINE void delfem2::Normalize_Quat(T q[]) {
  const T len = std::sqrt(
      q[0] * q[0] + q[1] * q[1] + q[2] * q[2] + q[3] * q[3]);
  const T invlen = 1 / len;
  q[0] *= invlen;
  q[1] *= invlen;
  q[2] *= invlen;
  q[3] *= invlen;
}
#ifdef DFM2_STATIC_LIBRARY
template void delfem2::Normalize_Quat(float q[]);
template void delfem2::Normalize_Quat(double q[]);
#endif

// -----------------------------------

template<typename T>
DFM2_INLINE void delfem2::Quat_Identity(T q[4]) {
  q[0] = 0;
  q[1] = 0;
  q[2] = 0;
  q[3] = 1;
}
#ifdef DFM2_STATIC_LIBRARY
template void delfem2::Quat_Identity(float q[4]);
template void delfem2::Quat_Identity(double q[4]);
#endif

// ----------------------------------

/**
 * @details for the relationship between quaternion and rotation matrix,
 * take a look at https://en.wikipedia.org/wiki/Quaternions_and_spatial_rotation
 */
template<typename T>
DFM2_INLINE void delfem2::QuatVec(
    T vo[3],
    const T q[4],
    const T vi[3]) {
  const T x2 = q[0] * q[0] * 2;
  const T y2 = q[1] * q[1] * 2;
  const T z2 = q[2] * q[2] * 2;
  const T xy = q[0] * q[1] * 2;
  const T yz = q[1] * q[2] * 2;
  const T zx = q[2] * q[0] * 2;
  const T xw = q[0] * q[3] * 2;
  const T yw = q[1] * q[3] * 2;
  const T zw = q[2] * q[3] * 2;
  vo[0] = (1 - y2 - z2) * vi[0] + (xy - zw) * vi[1] + (zx + yw) * vi[2];
  vo[1] = (xy + zw) * vi[0] + (1 - z2 - x2) * vi[1] + (yz - xw) * vi[2];
  vo[2] = (zx - yw) * vi[0] + (yz + xw) * vi[1] + (1 - x2 - y2) * vi[2];
}
#ifdef DFM2_STATIC_LIBRARY
template void delfem2::QuatVec(
    float vo[3], const float q[4], const float vi[3]);
template void delfem2::QuatVec(
    double vo[3], const double q[4], const double vi[3]);
#endif

// -------------------------------------

/**
 * multiply two quaternions (x,y,z,w)
 */
template<typename REAL>
DFM2_INLINE void delfem2::QuatQuat(
    REAL r[],
    const REAL p[],
    const REAL q[]) {
  r[0] = p[3] * q[0] + p[0] * q[3] + p[1] * q[2] - p[2] * q[1];
  r[1] = p[3] * q[1] - p[0] * q[2] + p[1] * q[3] + p[2] * q[0];
  r[2] = p[3] * q[2] + p[0] * q[1] - p[1] * q[0] + p[2] * q[3];
  r[3] = p[3] * q[3] - p[0] * q[0] - p[1] * q[1] - p[2] * q[2];
}
#ifdef DFM2_STATIC_LIBRARY
template void delfem2::QuatQuat(float r[], const float p[], const float q[]);
template void delfem2::QuatQuat(double r[], const double p[], const double q[]);
#endif


// ----------------------------------------

// transform vector with conjugate of quaternion
template<typename T>
DFM2_INLINE void delfem2::QuatConjVec(
    T vo[],
    const T q[],
    const T vi[]) {
  const T x2 = q[0] * q[0] * 2;
  const T y2 = q[1] * q[1] * 2;
  const T z2 = q[2] * q[2] * 2;
  const T xy = q[0] * q[1] * 2;
  const T yz = q[1] * q[2] * 2;
  const T zx = q[2] * q[0] * 2;
  const T xw = q[0] * q[3] * 2;
  const T yw = q[1] * q[3] * 2;
  const T zw = q[2] * q[3] * 2;
  vo[0] = (1 - y2 - z2) * vi[0] + (xy + zw) * vi[1] + (zx - yw) * vi[2];
  vo[1] = (xy - zw) * vi[0] + (1 - z2 - x2) * vi[1] + (yz + xw) * vi[2];
  vo[2] = (zx + yw) * vi[0] + (yz - xw) * vi[1] + (1 - x2 - y2) * vi[2];
}

// -------------------------

// copy quaternion
template<typename REAL>
DFM2_INLINE void delfem2::Copy_Quat(
    REAL r[],
    const REAL p[]) {
  std::memcpy(r, p, sizeof(REAL) * 4);
}
#ifdef DFM2_STATIC_LIBRARY
template void delfem2::Copy_Quat(float r[], const float p[]);
template void delfem2::Copy_Quat(double r[], const double p[]);
#endif

// -------------------------

template<typename T>
DFM2_INLINE void delfem2::Quat_Bryant(
    T q[4],
    T x,
    T y,
    T z) {
  constexpr T half = static_cast<T>(0.5);
  const T dqx[4] = {std::sin(x * half), 0.0, 0.0, std::cos(x * half)};
  const T dqy[4] = {0.0, std::sin(y * half), 0.0, std::cos(y * half)};
  const T dqz[4] = {0.0, 0.0, std::sin(z * half), std::cos(z * half)};
  T qtmp_yx[4];
  delfem2::QuatQuat(qtmp_yx, dqy, dqx);
  delfem2::QuatQuat(q, dqz, qtmp_yx);
}
#ifdef DFM2_STATIC_LIBRARY
template DFM2_INLINE void delfem2::Quat_Bryant(
    float q[4], float x, float y, float z);
template DFM2_INLINE void delfem2::Quat_Bryant(
    double q[4], double x, double y, double z);
#endif

// ------------------------

template<typename T>
DFM2_INLINE void delfem2::Quat_CartesianAngle(
    T q[4],
    const T a[3]) {
  constexpr T half = static_cast<T>(0.5);
  const T sqlen = a[0] * a[0] + a[1] * a[1] + a[2] * a[2];
  if (sqlen < 1.0e-10) {
    q[0] = half * a[0];
    q[1] = half * a[1];
    q[2] = half * a[2];
    q[3] = 1 - 0.125 * sqlen;
    return;
  }
  const T lena = std::sqrt(sqlen);
  q[0] = std::sin(lena * half) * a[0] / lena;
  q[1] = std::sin(lena * half) * a[1] / lena;
  q[2] = std::sin(lena * half) * a[2] / lena;
  q[3] = std::cos(lena * half);
}
#ifdef DFM2_STATIC_LIBRARY
template DFM2_INLINE void delfem2::Quat_CartesianAngle(float q[4], const float a[3]);
template DFM2_INLINE void delfem2::Quat_CartesianAngle(double q[4], const double a[3]);
#endif


// ---------------------------------------------------------------------

// functions for Quaternion above
// --------------------------------------------------------------------------------------------
// CQuat from here

namespace delfem2 {

template<typename T>
CQuat<T> operator+(const CQuat<T> &lhs, const CQuat<T> &rhs) {
  return CQuat<T>(lhs.w + rhs.w,
                  lhs.x + rhs.x,
                  lhs.y + rhs.y,
                  lhs.z + rhs.z);
}
#ifdef DFM2_STATIC_LIBRARY
template CQuat<double> operator+(
    const CQuat<double> &lhs, const CQuat<double> &rhs);
template CQuat<float> operator+(
    const CQuat<float> &lhs, const CQuat<float> &rhs);
#endif

// -----------------------

template<typename T>
CQuat<T> operator-(const CQuat<T> &lhs, const CQuat<T> &rhs) {
  return CQuat<T>(lhs.w - rhs.w,
                  lhs.x - rhs.x,
                  lhs.y - rhs.y,
                  lhs.z - rhs.z);
}
#ifdef DFM2_STATIC_LIBRARY
template CQuat<double> operator-(
    const CQuat<double> &lhs, const CQuat<double> &rhs);
template CQuat<float> operator-(
    const CQuat<float> &lhs, const CQuat<float> &rhs);
#endif

// -----------------------

template<typename T>
CQuat<T> operator*(const CQuat<T> &lhs, const CQuat<T> &rhs) {
  CQuat<T> q;
  QuatQuat(q.p, lhs.p, rhs.p);
  return q;
}
#ifdef DFM2_STATIC_LIBRARY
template CQuat<double> operator*(
    const CQuat<double> &lhs, const CQuat<double> &rhs);
template CQuat<float> operator*(
    const CQuat<float> &lhs, const CQuat<float> &rhs);
#endif

// ---------------------

template<typename T>
CQuat<T> operator*(const CQuat<T> &lhs, const T rhs) {
  CQuat<T> q(lhs.w * rhs,
             lhs.x * rhs,
             lhs.y * rhs,
             lhs.z * rhs);
  return q;
}
#ifdef DFM2_STATIC_LIBRARY
template CQuat<double> operator*(const CQuat<double> &lhs, double rhs);
template CQuat<float> operator*(const CQuat<float> &lhs, float rhs);
#endif

// ---------------------

template<typename T>
std::ostream &operator<<(std::ostream &output, const CQuat<T> &q) {
  output << q.w << " " << q.x << " " << q.y << " " << q.z;
  return output;
}

} // end namespace delfem2

// -------------------------

template<typename T>
void delfem2::CQuat<T>::normalize() {
  const T len = std::sqrt(
      p[0] * p[0] + p[1] * p[1] + p[2] * p[2] + p[3] * p[3]);
  const T leninv = 1 / len;
  p[0] *= leninv;
  p[1] *= leninv;
  p[2] *= leninv;
  p[3] *= leninv;
}
#ifdef DFM2_STATIC_LIBRARY
template void delfem2::CQuat<double>::normalize();
template void delfem2::CQuat<float>::normalize();
#endif

// ------------------------

template<typename T>
void delfem2::CQuat<T>::SetSmallerRotation() {
  if (w > 0) { return; }
  x = -x;
  y = -y;
  z = -z;
  w = -w;
}
#ifdef DFM2_STATIC_LIBRARY
template void delfem2::CQuat<double>::SetSmallerRotation();
template void delfem2::CQuat<float>::SetSmallerRotation();
#endif

// member function of CQuat
// ======================================================================

template<typename T>
delfem2::CQuat<T> delfem2::SphericalLinearInterp(
    const delfem2::CQuat<T> &q0,
    const delfem2::CQuat<T> &q1,
    T t) {
  const T qr = Dot_Quat(q0.p, q1.p);
  const T ss = 1 - qr * qr;

  if (ss == 0.0) { return q0; }
  const T sp = std::sqrt(ss);
  const T ph = std::acos(qr);
  const T pt = ph * t;
  const T t1 = std::sin(pt) / sp;
  const T t0 = std::sin(ph - pt) / sp;
  CQuat<T> q(
      t0 * q0.w + t1 * q1.w,
      t0 * q0.x + t1 * q1.x,
      t0 * q0.y + t1 * q1.y,
      t0 * q0.z + t1 * q1.z);
  return q;
}
#ifdef DFM2_STATIC_LIBRARY
template delfem2::CQuat<double> delfem2::SphericalLinearInterp(
    const delfem2::CQuat<double> &q0,
    const delfem2::CQuat<double> &q1,
    double t);
template delfem2::CQuat<float> delfem2::SphericalLinearInterp(
    const delfem2::CQuat<float> &q0,
    const delfem2::CQuat<float> &q1,
    float t);
#endif

/*
CQuat operator-(const CQuat& lhs, const CQuat& rhs){
    CQuat temp = lhs;
	temp -= rhs;
    return temp;
}

CQuat operator*(const CQuat& quat, double d){
    CQuat temp = quat;
	temp *= d;
    return temp;
}

CQuat operator*(double d, const CQuat& quat){
    CQuat temp = quat;
    temp *= d;
	return temp;
}

CQuat operator*(const CQuat& lhs, const CQuat& rhs){
    CQuat temp = lhs;
    temp *= rhs;
	return temp;
}

CVector3D Rotate(const CQuat& quat, const CVector3D& vec){
	CQuat tmp(0, vec);
	tmp = quat *  (tmp * quat.GetConjugate());  
//	tmp = quat.GetConjugate() *  tmp * quat ;   
	return tmp.GetVector();
}

CVector3D UnRotate(const CQuat& quat, const CVector3D& vec){
	CQuat tmp(0, vec);
	tmp = (quat.GetConjugate() *  tmp) * quat ; 
//	tmp = quat *  tmp * quat.GetConjugate();    
	return tmp.GetVector();
}


//////////////////////////////////////////////////////////////////////
// �����o�֐��̃t�����h�֐�
//////////////////////////////////////////////////////////////////////


bool operator==(const CQuat& lhs, const CQuat& rhs){
	if( fabs(lhs.real - rhs.real) < NEARLY_ZERO && lhs == rhs )	return true;
	else return false;
}

bool operator!=(const CQuat& lhs, const CQuat& rhs){
	if( lhs == rhs ) return false;
	else return true;
}

//////////////////////////////////////////////////////////////////////
// �����o�֐��̔�t�����h�֐��@
//////////////////////////////////////////////////////////////////////

CQuat& CQuat::operator = (const CQuat& rhs){
	if( this != &rhs ){
		real = rhs.real;
		vector = rhs.vector;
	}
	return *this;
}

CQuat& CQuat::operator *= (const CQuat& rhs){
	const double last_real = real;
	real = real * rhs.real - Dot( vector, rhs.vector );
	vector = ( (rhs.vector * last_real) + (vector * rhs.real) ) + Cross( vector, rhs.vector );
	return *this;
}

CQuat& CQuat::operator *= (double d){
	real *= d;
	vector *= d;
	return *this;
}

CQuat& CQuat::operator += (const CQuat& rhs){
	real += rhs.real;
	vector += rhs.vector;
	return *this;
}

CQuat& CQuat::operator -= (const CQuat& rhs){
	real -= rhs.real;
	vector -= rhs.vector;
	return *this;
}

CQuat CQuat::GetInverse() const
{
  double invslen = 1.0/SquareLength();
	CQuat tmp(real*invslen, vector*(-invslen));
	return tmp;
}

CQuat CQuat::GetConjugate() const{
	CQuat tmp(real, vector*(-1.0));
	return tmp;
}

void CQuat::RotMatrix44(double* m) const
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



void CQuat::RotMatrix33(double* m) const
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



void CQuat::AxisToQuat(const CVector3D &axis )
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

void CQuat::VectorTrans(const CVector3D& a_vector, const CVector3D& b_vector){
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

double CQuat::SquareLength() const
{
	return real*real + vector.DLength();
}


double CQuat::Length() const
{
	return sqrt( real*real + vector.DLength() );
}


*/


// CQuat
// -----------------------------
// std::vector


