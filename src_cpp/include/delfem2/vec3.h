/*
 * Copyright (c) 2019 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */


#ifndef VEC3_H
#define VEC3_H

#include <cassert>
#include <math.h>
#include <iostream>
#include <vector>

#define NEARLY_ZERO 1.e-16

double ScalarTripleProduct3D(const double a[], const double b[], const double c[]);
double Dot3D(const double a[], const double b[]);
double Length3D(const double v[3]);
void Normalize3D(double v[3]);
double SquareLength3D(const double v[3]);
double SquareDistance3D(const double p0[3], const double p1[3]);
double Distance3D(const double p0[3], const double p1[3]);
double TriArea3D(const double v1[3], const double v2[3], const double v3[3]);
double TetVolume3D(const double v1[3], const double v2[3], const double v3[3], const double v4[3] );
void UnitNormalAreaTri3D(double n[3], double& a,
                         const double v1[3], const double v2[3], const double v3[3]);
void NormalTri3D(double n[3],
                 const double v1[3], const double v2[3], const double v3[3]);
void Cross3D(double r[3],
             const double v1[3], const double v2[3]);
void InverseMat3(double Ainv[],
                 const double A[]);
void transposeMat3(double At[],
                   const double A[]);
void GetVertical2Vector3D(const double vec_n[3], double vec_x[3], double vec_y[3]);
void MatVec3(double y[3],
             const double m[9], const double x[3]);
void VecMat3D(const double x[3], const double m[9],  double y[3]);
void GetRotMatrix_Rodrigues3D(double rot[9],
                              const double n[3], double theta);

/////////////////////////////////////////////////////////
class CVector3;
double Dot(const CVector3 &arg1, const CVector3 &arg2);
CVector3 Cross(const CVector3& arg1, const CVector3& arg2);
CVector3 operator + (const CVector3& lhs, const CVector3& rhs);
CVector3 operator - (const CVector3& lhs, const CVector3& rhs);
CVector3 operator * (double d, const CVector3& rhs);
CVector3 operator * (const CVector3& vec, double d);
double operator * (const CVector3& lhs, const CVector3& rhs);
CVector3 operator / (const CVector3& vec, double d);
CVector3 operator ^ (const CVector3& lhs, const CVector3& rhs);
std::ostream &operator<<(std::ostream &output, const CVector3& v);
std::istream &operator>>(std::istream &input, CVector3& v);

//! 3 dimentional vector class
class CVector3  
{
public:
	CVector3(double vx, double vy, double vz) : x(vx), y(vy), z(vz){}
	CVector3(): x(0.0), y(0.0), z(0.0){}
	CVector3(const CVector3& rhs){ x = rhs.x; y = rhs.y; z = rhs.z; }
  CVector3(const double* prhs){ x = prhs[0]; y = prhs[1]; z = prhs[2]; }
  CVector3(const std::vector<double>& v){ x = v[0]; y = v[1]; z = v[2]; }
	virtual ~CVector3(){}

  std::vector<double> stlvec() const {
    std::vector<double> d(3);
    d[0] = x; d[1] = y; d[2] = z;
    return d;
  }
	void SetVector(double vx, double vy, double vz){ x = vx; y = vy; z = vz; }
  void CopyValueTo(double* v) const { v[0]=x; v[1]=y; v[2]=z; }
  void CopyValueToScale(double* v, double s) const { v[0]=x*s; v[1]=y*s; v[2]=z*s; }

	inline const CVector3 operator-() const{ return -1.0*(*this); }
	inline const CVector3 operator+() const{ return *this; }  
	inline CVector3& operator=(const CVector3& rhs){
		if( this != &rhs ){ x = rhs.x; y = rhs.y; z = rhs.z; }
		return *this;
	}
	inline CVector3& operator+=(const CVector3& rhs){
		x += rhs.x; y += rhs.y; z += rhs.z;
		return *this;
	}
	inline CVector3& operator-=(const CVector3& rhs){
		x -= rhs.x; y -= rhs.y; z -= rhs.z;
		return *this;
	}
	inline CVector3& operator*=(double d){
		x *= d; y *= d; z *= d;
		return *this;
	}
	inline CVector3& operator/=(double d){
		if( fabs(d) < NEARLY_ZERO ){ return *this; }
		x /= d; y /= d; z /= d;
		return *this;
	}
  inline double operator[](int i) const{
    if( i == 0 ) return x;
    if( i == 1 ) return y;
    if( i == 2 ) return z;
    return 0;
  }
  inline double& operator[](int i){
    if( i == 0 ) return x;
    if( i == 1 ) return y;
    if( i == 2 ) return z;
    assert(0);
    return x;
  }  
	inline CVector3 operator+(){ return *this; }
	inline CVector3 operator-(){ return CVector3(-x,-y,-z); }

	friend bool operator==(const CVector3&, const CVector3&);
	friend bool operator!=(const CVector3&, const CVector3&);

	friend CVector3 Cross(const CVector3&, const CVector3&);
	friend double Dot(const CVector3&, const CVector3&);
  CVector3 Normalize() const {
    CVector3 r = (*this);
    r.SetNormalizedVector();
    return r;
  }
	inline double Length()  const{ return sqrt( x*x+y*y+z*z ); }
	inline double DLength() const{ return x*x+y*y+z*z; }
	void SetNormalizedVector();
	void SetZero();
  void Print() const {
    std::cout << x << " " << y << " " << z << std::endl;
  }
  bool isNaN() const{
    double s=x+y+z;
    return !(s > s-1.0);
  }
  static CVector3 Axis(int idim){
    CVector3 r(0,0,0);
    if( idim >= 0 && idim < 3) { r[idim] = 1; }
    return r;
  }
public:
	double x;	//!< x axis coordinate
	double y;	//!< y axis coordinate
	double z;	//!< z axis coordinate
};
//////////////////////////////////////////////////////////////////////////////


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// rule about naming, the method starts "Set" change it self (not const)

CVector3 MatVec(double M[9], const CVector3& v);
CVector3 screenProjection(const CVector3& v,
                          const float* mMV, const float* mPj);
// opposite to the screen normal
CVector3 screenDepthDirection(const CVector3& v,
                              const float* mMV,
                              const float* mPj);
CVector3 screenUnProjection(const CVector3& v,
                            const float* mMV,
                            const float* mPj);
CVector3 screenUnProjectionDirection(const CVector3& v,
                                     const float* mMV,
                                     const float* mPj);
CVector3 screenDepthDirection(const CVector3& v,
                              const float* mMV,
                              const float* mPj);

double ScalarTripleProduct(const CVector3& a, const CVector3& b, const CVector3& c);
bool operator == (const CVector3& lhs, const CVector3& rhs);
bool operator != (const CVector3& lhs, const CVector3& rhs);
double Height(const CVector3& v1, const CVector3& v2, const CVector3& v3, const CVector3& v4);
void GetVertical2Vector (const CVector3& vec_n, CVector3& vec_x, CVector3& vec_y);

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

CVector3 nearest_Line_Point(const CVector3& p, // point
                              const CVector3& s, // source
                              const CVector3& d);
CVector3 nearest_Line_Point(double& t,
                              const CVector3& p, // point
                              const CVector3& s, // source
                              const CVector3& d); // direction
CVector3 nearest_LineSeg_Point(const CVector3& p, // point
                                 const CVector3& s, // source
                                 const CVector3& e); // end
CVector3 nearest_LineSeg_Point(double& t,
                                 const CVector3& p, // point
                                 const CVector3& s,
                                 const CVector3& e); // direction
void nearest_Line_Line(double& D, CVector3& Da, CVector3& Db,
                         const CVector3& pa_, const CVector3& va,
                         const CVector3& pb_, const CVector3& vb);
void nearest_Line_Line(double& D, CVector3& Da, CVector3& Db,
                         double& ta, double& tb,
                         const CVector3& pa_, const CVector3& va,
                         const CVector3& pb_, const CVector3& vb);
CVector3 nearest_Plane_Point(const CVector3& p, // point
                             const CVector3& o, // origin
                             const CVector3& n); // normal
CVector3 nearest_Origin_Tri(double& r0,
                            double& r1,
                            const CVector3& q0,
                            const CVector3& q1,
                            const CVector3& q2);
CVector3 nearst_Origin_Quad(double& s0, double& s1,
                            const CVector3& p,
                            const CVector3& q0,
                            const CVector3& q1,
                            const CVector3& q2,
                            const CVector3& q3);
void Nearest_Line_Circle(CVector3& p0,
                         CVector3& q0,
                         const CVector3& src,
                         const CVector3& dir,
                         const CVector3& org,
                         const CVector3& normal,
                         double rad);

////////////////////////////////////////////////////////////////

bool intersection_Plane_Line(CVector3& p0, double& r0, double& r1, double& r2,
                             double eps,
                             const CVector3& src, const CVector3& dir,
                             const CVector3& q0, const CVector3& q1, const CVector3& q2);
bool intersection_Point_Quad(CVector3& psec, double& s0, double& s1,
                             const CVector3& src, const CVector3& dir,
                             const CVector3& q0, const CVector3& q1, const CVector3& q2, const CVector3& q3);
CVector3 intersection_Plane_Line(const CVector3& o, // one point on plane
                                 const CVector3& n, // plane normal
                                 const CVector3& s, // one point on line
                                 const CVector3& d); // direction of line

////////////////////////////////////////////////////////////////

bool barycentricCoord_Origin_Tet(double& r0,
                                 double& r1,
                                 double& r2,
                                 const CVector3& p0,
                                 const CVector3& p1,
                                 const CVector3& p2,
                                 const CVector3& p3);
bool barycentricCoord_Origin_Pyramid(double& r0,
                                     double& r1,
                                     double& r2,
                                     const CVector3& p0,
                                     const CVector3& p1,
                                     const CVector3& p2,
                                     const CVector3& p3,
                                     const CVector3& p4);
bool barycentricCoord_Origin_Wedge(double& r0,
                                   double& r1,
                                   double& r2,
                                   const CVector3& p0,
                                   const CVector3& p1,
                                   const CVector3& p2,
                                   const CVector3& p3,
                                   const CVector3& p4,
                                   const CVector3& p5);

////////////////////////////////////////////////////////////////

bool IsInside_Orgin_BoundingBoxPoint4(const CVector3& p0,
                                    const CVector3& p1,
                                    const CVector3& p2,
                                    const CVector3& p3);

bool IsInside_Orgin_BoundingBoxPoint5(const CVector3& p0,
                                    const CVector3& p1,
                                    const CVector3& p2,
                                    const CVector3& p3,
                                    const CVector3& p4);

bool IsInside_Orgin_BoundingBoxPoint6(const CVector3& p0,
                                    const CVector3& p1,
                                    const CVector3& p2,
                                    const CVector3& p3,
                                    const CVector3& p4,
                                    const CVector3& p5);

////////////////////////////////////////////////////////////////

double volume_OrgTet(const CVector3& v1,
                     const CVector3& v2,
                     const CVector3& v3 );
double volume_Tet(const CVector3& v0,
                  const CVector3& v1,
                  const CVector3& v2,
                  const CVector3& v3 );
double volume_Pyramid(const CVector3& p0,
                      const CVector3& p1,
                      const CVector3& p2,
                      const CVector3& p3,
                      const CVector3& p4);
double volume_Wedge(const CVector3& p0,
                    const CVector3& p1,
                    const CVector3& p2,
                    const CVector3& p3,
                    const CVector3& p4,
                    const CVector3& p5);

//////////////////////////////////////////////////////////

double SolidAngleTri(const CVector3& v1,
                     const CVector3& v2,
                     const CVector3& v3);

void Cross( CVector3& lhs, const CVector3& v1, const CVector3& v2 );
double TriArea(const CVector3& v1, const CVector3& v2, const CVector3& v3);
double SquareTriArea(const CVector3& v1, const CVector3& v2, const CVector3& v3);
double SquareDistance(const CVector3& ipo0, const CVector3& ipo1);
double SquareLength(const CVector3& point);
double Length(const CVector3& point);
double Distance(const CVector3& ipo0, const CVector3& ipo1);

double SqareLongestEdgeLength(const CVector3& ipo0,
                              const CVector3& ipo1,
                              const CVector3& ipo2,
                              const CVector3& ipo3 );

double LongestEdgeLength(const CVector3& ipo0,
                         const CVector3& ipo1,
                         const CVector3& ipo2,
                         const CVector3& ipo3 );

double SqareShortestEdgeLength(const CVector3& ipo0,
                               const CVector3& ipo1,
                               const CVector3& ipo2,
                               const CVector3& ipo3 );

double ShortestEdgeLength(const CVector3& ipo0,
                          const CVector3& ipo1,
                          const CVector3& ipo2,
                          const CVector3& ipo3 );

void Normal(CVector3& vnorm,
            const CVector3& v1,
            const CVector3& v2,
            const CVector3& v3);

CVector3 Normal(const CVector3& v1,
                const CVector3& v2,
                const CVector3& v3);

void UnitNormal(CVector3& vnorm,
                const CVector3& v1,
                const CVector3& v2,
                const CVector3& v3);

//! check if Delaunay condition satisfied
// 0 : p3 is inside circum circle on the p0,p1,p2
// 1 :       on
// 2 :       outsdie
int DetDelaunay(const CVector3& p0,
                const CVector3& p1,
                const CVector3& p2,
                const CVector3& p3);

double SquareCircumradius(const CVector3& ipo0,
                          const CVector3& ipo1,
                          const CVector3& ipo2,
                          const CVector3& ipo3);

CVector3 CircumCenter(const CVector3& ipo0,
                      const CVector3& ipo1,
                      const CVector3& ipo2,
                      const CVector3& ipo3);

double Circumradius(const CVector3& ipo0,
                    const CVector3& ipo1,
                    const CVector3& ipo2,
                    const CVector3& ipo3);
CVector3 RotateVector(const CVector3& vec0, const CVector3& rot );
CVector3 RandVector();
CVector3 RandUnitVector();
CVector3 RandGaussVector();

void MeanValueCoordinate(double w[3],
                         const CVector3& v0,
                         const CVector3& v1,
                         const CVector3& v2);

CVector3 ProjectPointOnTriangle(const CVector3 &p0,
                                const CVector3 &tri_p1, const CVector3 &tri_p2, const CVector3 &tri_p3);
bool isPointInsideTriangle(const CVector3 &p0,
                           const CVector3 &tri_p1, const CVector3 &tri_p2, const CVector3 &tri_p3);
bool isPointSameSide(const CVector3 &p0, const CVector3 &p1,
                     const CVector3 &line_p0, const CVector3 &line_p1);
bool isRayIntersectingTriangle(const CVector3 &line0, const CVector3 &line1,
                               const CVector3 &tri0, const CVector3 &tri1, const CVector3 &tri2,
                               CVector3 &intersectionPoint);


///////////////////
// here starts std::vector<CVector3>


double TetVolume( int iv1, int iv2, int iv3, int iv4,
                 const std::vector<CVector3>& node);

double TriArea(const int iv1, const int iv2, const int iv3,
               const std::vector<CVector3>& node );

CVector3 QuadBilinear(int iq, double r0, double r1,
                      std::vector<int>& aQuad,
                      std::vector<CVector3>& aPoint);

void ConvexHull(std::vector<int>& aTri,
                const std::vector<CVector3>& aXYZ);

inline CVector3 cg_Tri(int itri,
                       const std::vector<unsigned int>& aTri,
                       const std::vector<double>& aXYZ)
{
  CVector3 p;
  int i0 = aTri[itri*3+0];
  int i1 = aTri[itri*3+1];
  int i2 = aTri[itri*3+2];
  p.x = (aXYZ[i0*3+0]+aXYZ[i1*3+0]+aXYZ[i2*3+0])/3.0;
  p.y = (aXYZ[i0*3+1]+aXYZ[i1*3+1]+aXYZ[i2*3+1])/3.0;
  p.z = (aXYZ[i0*3+2]+aXYZ[i1*3+2]+aXYZ[i2*3+2])/3.0;
  return p;
}

inline CVector3 normalTri(int itri,
                          const std::vector<int>& aTri,
                          const std::vector<double>& aXYZ)
{
  int i0 = aTri[itri*3+0];
  int i1 = aTri[itri*3+1];
  int i2 = aTri[itri*3+2];
  CVector3 p0(aXYZ[i0*3+0],aXYZ[i0*3+1],aXYZ[i0*3+2]);
  CVector3 p1(aXYZ[i1*3+0],aXYZ[i1*3+1],aXYZ[i1*3+2]);
  CVector3 p2(aXYZ[i2*3+0],aXYZ[i2*3+1],aXYZ[i2*3+2]);
  return (p1-p0)^(p2-p0);
}

std::ostream &operator<<(std::ostream &output, const std::vector<CVector3>& aV);
std::istream &operator>>(std::istream &input, std::vector<CVector3>& aV);




#endif // PHYSICS_VECTOR_3D_H
