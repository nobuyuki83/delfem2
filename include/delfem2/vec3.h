 /*
 * Copyright (c) 2019 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */


#ifndef DFM2_VEC3_H
#define DFM2_VEC3_H

#include <cassert>
#include <cmath>
#include <iostream>
#include <vector>
#include "delfem2/dfm2_inline.h"

#define NEARLY_ZERO 1.e-16

namespace delfem2 {

template <typename T>
DFM2_INLINE T Distance3(const T p0[3], const T p1[3]);

template <typename T>
DFM2_INLINE T SquareDistance3(const T p0[3], const T p1[3]);

template <typename T>
DFM2_INLINE T Length3(const T v[3]);

template <typename T>
DFM2_INLINE T SquareLength3(const T v[3]);
  
template <typename T>
DFM2_INLINE T Dot3(const T a[3], const T b[3]);

template <typename T>
DFM2_INLINE void Normalize3(T v[3]);

template <typename REAL>
void AverageTwo3(REAL po[3],
    const REAL p0[3], const REAL p1[3]);

template <typename REAL>
void AverageFour3(REAL po[3],
    const REAL p0[3], const REAL p1[3], const REAL p2[3], const REAL p3[3]);

/**
 * @func add values for 3-array (vo += vi)
 * @tparam REAL float and double
 * @param vo (out)
 * @param vi (in)
 */
template <typename REAL>
DFM2_INLINE void Add3(
    REAL vo[3],
    const REAL vi[3]);

// above: functions general for any dimensions
// ------------------------------------------------
// below: functions specific to 3 dimension
  
template <typename T>
T Volume_Tet3(const T v1[3], const T v2[3], const T v3[3], const T v4[3]);
  
template <typename T>
void Cross3(T r[3],
    const T v1[3], const T v2[3]);

template <typename T>
T Area_Tri3(const T v1[3], const T v2[3], const T v3[3]);
  
template <typename T>
T ScalarTripleProduct3(const T a[], const T b[], const T c[]);

template <typename T>
void NormalTri3(T n[3],
    const T v1[3], const T v2[3], const T v3[3]);


template <typename REAL>
void UnitNormalAreaTri3(
    REAL n[3], REAL& a,
    const REAL v1[3], const REAL v2[3], const REAL v3[3]);

// --------------------------------
  
DFM2_INLINE void GetVertical2Vector3D(
    const double vec_n[3],
    double vec_x[3],
    double vec_y[3]);

// -------------------------------------------------------------

template <typename T>
class CVec3;
  
template <typename T>
CVec3<T> operator + (const CVec3<T>& lhs, const CVec3<T>& rhs);
  
template <typename T>
CVec3<T> operator - (const CVec3<T>& lhs, const CVec3<T>& rhs);
  
template <typename T>
CVec3<T> operator * (T d, const CVec3<T>& rhs);
  
template <typename T>
CVec3<T> operator * (const CVec3<T>& vec, T d);
  
template <typename T>
T operator * (const CVec3<T>& lhs, const CVec3<T>& rhs);
  
template <typename T>
CVec3<T> operator / (const CVec3<T>& vec, T d);
  
template <typename T>
CVec3<T> operator ^ (const CVec3<T>& lhs, const CVec3<T>& rhs);
  
template <typename T>
std::ostream &operator<<(std::ostream &output, const CVec3<T>& v);
  
template <typename T>
std::istream &operator>>(std::istream &input, CVec3<T>& v);
  
template <typename T>
std::ostream &operator<<(std::ostream &output, const std::vector<CVec3<T> >& aV);

template <typename T>
std::istream &operator>>(std::istream &input, std::vector<CVec3<T> >& aV);
  
/**
 * @class 3 dimentional vector class
 */
template <typename T>
class CVec3
{
public:
  CVec3(T vx, T vy, T vz) : p{vx,vy,vz} {}
  CVec3(): p{0.0, 0.0, 0.0} {}
	CVec3(const CVec3& rhs){ p[0] = rhs.p[0]; p[1] = rhs.p[1]; p[2] = rhs.p[2]; }
  CVec3(const double* prhs){ p[0] = prhs[0]; p[1] = prhs[1]; p[2] = prhs[2]; }
  CVec3(const float* prhs){ p[0] = prhs[0]; p[1] = prhs[1]; p[2] = prhs[2]; }
  CVec3(const std::vector<double>& v){ p[0] = v[0]; p[1] = v[1]; p[2] = v[2]; }
	virtual ~CVec3(){}

  std::vector<double> stlvec() const {
    std::vector<double> d = {p[0], p[1], p[2]};
    return d;
  }
	void SetVector(double vx, double vy, double vz){ p[0] = vx; p[1] = vy; p[2] = vz; }
  void CopyTo(double* v) const { v[0]=p[0]; v[1]=p[1]; v[2]=p[2]; }
  void CopyToScale(double* v, double s) const { v[0]=p[0]*s; v[1]=p[1]*s; v[2]=p[2]*s; }
  void AddToScale(double* v, double s) const { v[0]+=p[0]*s; v[1]+=p[1]*s; v[2]+=p[2]*s; }

	inline const CVec3 operator-() const{ return ((T)(-1))*(*this); }
	inline const CVec3 operator+() const{ return *this; }
	inline CVec3& operator=(const CVec3& rhs){
		if( this != &rhs ){ p[0]= rhs.p[0]; p[1] = rhs.p[1]; p[2] = rhs.p[2]; }
		return *this;
	}
	inline CVec3& operator+=(const CVec3& rhs){
		p[0] += rhs.p[0];
    p[1] += rhs.p[1];
    p[2] += rhs.p[2];
		return *this;
	}
	inline CVec3& operator-=(const CVec3& rhs){
		p[0] -= rhs.p[0];
    p[1] -= rhs.p[1];
    p[2] -= rhs.p[2];
		return *this;
	}
	inline CVec3& operator*=(double d){
		p[0] *= d;
    p[1] *= d;
    p[2] *= d;
		return *this;
	}
	inline CVec3& operator/=(double d){
		if( fabs(d) < NEARLY_ZERO ){ return *this; }
		p[0] /= d;
    p[1] /= d;
    p[2] /= d;
		return *this;
	}
  inline double operator[](int i) const{
    if( i == 0 ) return p[0];
    if( i == 1 ) return p[1];
    if( i == 2 ) return p[2];
    return 0;
  }
  inline T& operator[](int i){
    if( i == 0 ) return p[0];
    if( i == 1 ) return p[1];
    if( i == 2 ) return p[2];
    assert(0);
    return p[0];
  }  
	inline CVec3 operator+(){ return *this; }
	inline CVec3 operator-(){ return CVec3(-p[0],-p[1],-p[2]); }
  inline double x() const { return p[0]; }
  inline double y() const { return p[1]; }
  inline double z() const { return p[2]; }

  CVec3 Normalize() const {
    CVec3 r = (*this);
    r.SetNormalizedVector();
    return r;
  }
  CVec3<float> Float() const {
    return CVec3<float>((float)p[0], (float)p[1], (float)p[2]);
  }
  CVec3<double> Double() const {
    return CVec3<double>((double)p[0], (double)p[1], (double)p[2]);
  }
	inline double Length()  const{ return sqrt( p[0]*p[0]+p[1]*p[1]+p[2]*p[2] ); }
	inline double DLength() const{ return p[0]*p[0]+p[1]*p[1]+p[2]*p[2]; }
	void SetNormalizedVector();
	void SetZero();
  void Print() const {
    std::cout <<p[0]<< " " << p[1] << " " << p[2] << std::endl;
  }
  bool isNaN() const{
    double s=p[0]+p[1]+p[2];
    return !(s > s-1.0);
  }
  static CVec3 Axis(unsigned int idim){
    CVec3 r(0,0,0);
    if( idim < 3) { r[idim] = 1; }
    return r;
  }
  T* data() { return p; }
  const T* data() const { return p; }
  //
  template <typename DIST, typename ENG>
  void SetRandom(DIST& dist, ENG& eng){
    p[0] = dist(eng);
    p[1] = dist(eng);
    p[2] = dist(eng);
  }
  //
  template <typename DIST, typename ENG>
  static CVec3 Random(DIST& dist, ENG& eng){
    return CVec3(dist(eng),dist(eng),dist(eng));
  }
//  void SetRandom();
public:
  T p[3];
};
using CVec3d = CVec3<double>;
using CVec3f = CVec3<float>;

// ------------------------------
  
template <typename T>
T Dot(const CVec3<T>& arg1, const CVec3<T>& arg2);

template <typename T>
CVec3<T> Cross(const CVec3<T>& arg1, const CVec3<T>& arg2);

  
// --------------------------------------------------------------------------
// rule about naming, the method starts "Set" change it self (not const)
  
template <typename T>
CVec3<T> mult_GlAffineMatrix(const float* m,
                             const CVec3<T>& p);

template <typename T>
CVec3<T> solve_GlAffineMatrix(const float* m,
                              const CVec3<T>& p);
  
template <typename T>
CVec3<T> solve_GlAffineMatrixDirection(const float* m,
                                       const CVec3<T>& v);

template <typename T>
CVec3<T> Mat3Vec(const T M[ 9], const CVec3<T>& v);
  
template <typename T>
CVec3<T> Mat4Vec(const T M[16], const CVec3<T>& v);
  
template <typename T>
DFM2_INLINE CVec3<T> QuatVec
 (const T quat[4],
  const CVec3<T>& v0);
  
template <typename REAL>
CVec3<REAL> QuatConjVec(
    const REAL quat[4],
    const CVec3<REAL>& v0);
  
template <typename T>
CVec3<T> screenProjection(
    const CVec3<T>& v,
    const float* mMV,
    const float* mPj);
  
// opposite to the screen normal
template <typename T>
CVec3<T> screenDepthDirection(
    const CVec3<T>& v,
    const float* mMV,
    const float* mPj);
  
template <typename T>
CVec3<T> screenUnProjection(
    const CVec3<T>& v,
    const float* mMV,
    const float* mPj);
  
template <typename T>
CVec3<T> screenUnProjectionDirection(
    const CVec3<T>& v,
    const float* mMV,
    const float* mPj);
  
template <typename T>
CVec3<T> screenDepthDirection(
    const CVec3<T>& v,
    const float* mMV,
    const float* mPj);

template <typename T>
double ScalarTripleProduct(const CVec3<T>& a, const CVec3<T>& b, const CVec3<T>& c);
  
template <typename T>
bool operator == (const CVec3<T>& lhs, const CVec3<T>& rhs);
  
template <typename T>
bool operator != (const CVec3<T>& lhs, const CVec3<T>& rhs);
  
template <typename T>
double Height(const CVec3<T>& v1, const CVec3<T>& v2, const CVec3<T>& v3, const CVec3<T>& v4);
  
template <typename T>
void GetVertical2Vector (const CVec3<T>& vec_n, CVec3<T>& vec_x, CVec3<T>& vec_y);

// ---------------------------------------------------------

template <typename T>
bool barycentricCoord_Origin_Tet(
    double& r0,
    double& r1,
    double& r2,
    const CVec3<T>& p0,
    const CVec3<T>& p1,
    const CVec3<T>& p2,
    const CVec3<T>& p3);
  
template <typename T>
bool barycentricCoord_Origin_Pyramid(
    double& r0,
    double& r1,
    double& r2,
    const CVec3<T>& p0,
    const CVec3<T>& p1,
    const CVec3<T>& p2,
    const CVec3<T>& p3,
    const CVec3<T>& p4);
  
template <typename T>
bool barycentricCoord_Origin_Wedge(
    double& r0,
    double& r1,
    double& r2,
    const CVec3<T>& p0,
    const CVec3<T>& p1,
    const CVec3<T>& p2,
    const CVec3<T>& p3,
    const CVec3<T>& p4,
    const CVec3<T>& p5);

template <typename T>
CVec3<T> positionBarycentricCoord_Pyramid(
    double r0,
    double r1,
    double r2,
    const CVec3<T>& p0,
    const CVec3<T>& p1,
    const CVec3<T>& p2,
    const CVec3<T>& p3,
    const CVec3<T>& p4);
  
template <typename T>
CVec3<T> positionBarycentricCoord_Wedge(
    double r0,
    double r1,
    double r2,
    const CVec3<T>& p0,
    const CVec3<T>& p1,
    const CVec3<T>& p2,
    const CVec3<T>& p3,
    const CVec3<T>& p4,
    const CVec3<T>& p5);
  
template <typename T>
void iteration_barycentricCoord_Origin_Solid(
    double& r0,
    double& r1,
    double& r2,
    const CVec3<T>& q, // q=positionBarycentricCoord_Wedge
    const CVec3<T>& dpdr0,
    const CVec3<T>& dpdr1,
    const CVec3<T>& dpdr2,
    double damp=1.0);

  
// -----------------------------------------------

template <typename T>
T Volume_OrgTet(
    const CVec3<T>& v1,
    const CVec3<T>& v2,
    const CVec3<T>& v3 );
  
template <typename T>
T Volume_Tet(
    const CVec3<T>& v0,
    const CVec3<T>& v1,
    const CVec3<T>& v2,
    const CVec3<T>& v3 );

template <typename T>
double Volume_Pyramid(
    const CVec3<T>& p0,
    const CVec3<T>& p1,
    const CVec3<T>& p2,
    const CVec3<T>& p3,
    const CVec3<T>& p4);

template <typename T>
double Volume_Wedge(
    const CVec3<T>& p0,
    const CVec3<T>& p1,
    const CVec3<T>& p2,
    const CVec3<T>& p3,
    const CVec3<T>& p4,
    const CVec3<T>& p5);

// ---------------------------------------------

template <typename T>
double SolidAngleTri(const CVec3<T>& v1,
                     const CVec3<T>& v2,
                     const CVec3<T>& v3);

template <typename T>
void Cross( CVec3<T>& lhs, const CVec3<T>& v1, const CVec3<T>& v2 );
  
template <typename T>
double Area_Tri(const CVec3<T>& v1, const CVec3<T>& v2, const CVec3<T>& v3);
  
template <typename T>
double SquareTriArea(const CVec3<T>& v1, const CVec3<T>& v2, const CVec3<T>& v3);
  
template <typename T>
double SquareDistance(const CVec3<T>& ipo0, const CVec3<T>& ipo1);

template <typename T>
double SquareLength(const CVec3<T>& point);
  
template <typename T>
double Length(
    const CVec3<T>& point);
  
template <typename T>
double Distance(
    const CVec3<T>& ipo0,
    const CVec3<T>& ipo1);

template <typename T>
double SqareLongestEdgeLength(
    const CVec3<T>& ipo0,
    const CVec3<T>& ipo1,
    const CVec3<T>& ipo2,
    const CVec3<T>& ipo3 );

template <typename T>
double LongestEdgeLength(
    const CVec3<T>& ipo0,
    const CVec3<T>& ipo1,
    const CVec3<T>& ipo2,
    const CVec3<T>& ipo3 );

template <typename T>
double SqareShortestEdgeLength(
    const CVec3<T>& ipo0,
    const CVec3<T>& ipo1,
    const CVec3<T>& ipo2,
    const CVec3<T>& ipo3 );

template <typename T>
double ShortestEdgeLength(
    const CVec3<T>& ipo0,
    const CVec3<T>& ipo1,
    const CVec3<T>& ipo2,
    const CVec3<T>& ipo3 );

template <typename T>
void Normal(CVec3<T>& vnorm,
            const CVec3<T>& v1,
            const CVec3<T>& v2,
            const CVec3<T>& v3);

template <typename T>
CVec3<T> Normal(const CVec3<T>& v1,
                const CVec3<T>& v2,
                const CVec3<T>& v3);

template <typename T>
void UnitNormal(CVec3<T>& vnorm,
                const CVec3<T>& v1,
                const CVec3<T>& v2,
                const CVec3<T>& v3);
  
template <typename T>
CVec3<T> UnitNormal(const CVec3<T>& v1,
                    const CVec3<T>& v2,
                    const CVec3<T>& v3);
  
/**
 * @function check if Delaunay condition satisfied
 * @return
 * 0 : p3 is inside circum circle on the p0,p1,p2
 * 1 :       on
 * 2 :       outsdie
 */
template <typename T>
int DetDelaunay(const CVec3<T>& p0,
                const CVec3<T>& p1,
                const CVec3<T>& p2,
                const CVec3<T>& p3);

template <typename T>
double SquareCircumradius(const CVec3<T>& ipo0,
                          const CVec3<T>& ipo1,
                          const CVec3<T>& ipo2,
                          const CVec3<T>& ipo3);

template <typename T>
CVec3<T> CircumCenter(const CVec3<T>& ipo0,
                      const CVec3<T>& ipo1,
                      const CVec3<T>& ipo2,
                      const CVec3<T>& ipo3);

template <typename T>
double Circumradius(const CVec3<T>& ipo0,
                    const CVec3<T>& ipo1,
                    const CVec3<T>& ipo2,
                    const CVec3<T>& ipo3);

template <typename T>
CVec3<T> RotateVector(const CVec3<T>& vec0, const CVec3<T>& rot );
  
template <typename T>
CVec3<T> RandVector();
  
template <typename T>
CVec3<T> RandUnitVector();
  
template <typename T>
CVec3<T> RandGaussVector();

template <typename T>
void MeanValueCoordinate(
    double w[3],
    const CVec3<T>& v0,
    const CVec3<T>& v1,
    const CVec3<T>& v2);

// ----------------------------------------------------------
// here starts std::vector<CVector3>
  
template <typename T>
double Volume_Tet(
    int iv1,
    int iv2,
    int iv3,
    int iv4,
    const std::vector<CVec3<T> >& aPoint);

template <typename T>
double Area_Tri(
    const int iv1,
    const int iv2,
    const int iv3,
    const std::vector<CVec3<T> >& node );

template <typename T>
CVec3<T> CG_Tri3(
    unsigned int itri,
    const std::vector<unsigned int>& aTri,
    const std::vector<T>& aXYZ)
{
  const unsigned int i0 = aTri[itri*3+0];
  const unsigned int i1 = aTri[itri*3+1];
  const unsigned int i2 = aTri[itri*3+2];
  CVec3<T> p;
  p.p[0] = (aXYZ[i0*3+0]+aXYZ[i1*3+0]+aXYZ[i2*3+0])/3.0;
  p.p[1] = (aXYZ[i0*3+1]+aXYZ[i1*3+1]+aXYZ[i2*3+1])/3.0;
  p.p[2] = (aXYZ[i0*3+2]+aXYZ[i1*3+2]+aXYZ[i2*3+2])/3.0;
  return p;
}

template <typename T>
inline CVec3<T> Normal_Tri3(
    unsigned int itri,
    const std::vector<unsigned int>& aTri,
    const std::vector<T>& aXYZ)
{
  const unsigned int i0 = aTri[itri*3+0];
  const unsigned int i1 = aTri[itri*3+1];
  const unsigned int i2 = aTri[itri*3+2];
  const CVec3<T> p0(aXYZ[i0*3+0],aXYZ[i0*3+1],aXYZ[i0*3+2]);
  const CVec3<T> p1(aXYZ[i1*3+0],aXYZ[i1*3+1],aXYZ[i1*3+2]);
  const CVec3<T> p2(aXYZ[i2*3+0],aXYZ[i2*3+1],aXYZ[i2*3+2]);
  return (p1-p0)^(p2-p0);
}

} // end namespace delfem2

#ifdef DFM2_HEADER_ONLY
#  include "delfem2/vec3.cpp"
#endif


#endif // VEC3_H
