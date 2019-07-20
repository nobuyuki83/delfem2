/*
 * Copyright (c) 2019 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */

#ifndef VEC2_H
#define VEC2_H

#include <cstdlib>
#include <cassert>
#include <iostream>
#include <cmath>
#include <vector>

double TriArea2D(const double v1[2], const double v2[2], const double v3[2]);
double SqDistance2D(const double v1[2], const double v2[2]);
double Distance2D(const double v1[2], const double v2[2]);
void noise2D(double noise[2]);
bool InverseMat2(double invB[4], const double B[4]);
void matMat2(double AB[4], const double A[4], const double B[4]);
void MatVec2(double w[2], const double A[4], const double v[2]);
void setNormalized2(double w[2]);
double dot2(const double w[2], const double v[2]);
double len2(const double v[2]);
double sqLen2(const double v[2]);
void gramian2(double AtA[3], const double A[4]);
void VLVt2(double A[4], double l0, double l1, const double V[4]);
void RotationalComponentOfMatrix2(double R[4], const double M[4]);

////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////

class CVector2;

CVector2 operator*(double, const CVector2&);
CVector2 operator*(const CVector2&, double);
std::ostream &operator<<(std::ostream &output, const CVector2& v);
std::istream &operator>>(std::istream &input, CVector2& v);

//! 2 dimensional vector class
class CVector2{
public:
	CVector2(){
    x = 0; y = 0;
  }
	CVector2( const CVector2& rhs ){
		this->x = rhs.x;
		this->y = rhs.y;
	}
	CVector2(double x, double y){
		this->x = x;
		this->y = y;
	}
	
	friend double Dot(const CVector2&, const CVector2&);
  CVector2 operator-() const{
    return CVector2(-x,-y);
  }
	inline CVector2& operator+=(const CVector2& rhs){
		x += rhs.x;
		y += rhs.y;
		return *this;
	}
	inline CVector2& operator-=(const CVector2& rhs){
		x -= rhs.x;
		y -= rhs.y;
		return *this;
	}
	inline CVector2& operator*=(double scale){
		x *= scale;
		y *= scale;
		return *this;
  }
  inline CVector2& operator/=(double d){
    if (fabs(d) < 1.0e-6){ assert(0); return *this; }
    x /= d; y /= d;
    return *this;
  }
	inline CVector2 operator+(const CVector2& rhs) const {
		CVector2 v = *this;
		return v += rhs;
	}
	inline CVector2 operator-(const CVector2& rhs) const {
		CVector2 v = *this;
		return v -= rhs;
	}
  inline double operator[](int i) const{
    if (i==0) return x;
    if (i==1) return y;
    return 0;
  }
  inline double& operator[](int i){
    if (i==0) return x;
    if (i==1) return y;
    assert(0);
    return x;
  }
	//! normalize length
	inline void SetNormalizedVector(){
		const double mag = Length();
		x /= mag;
		y /= mag;
	}
  CVector2 Normalize() const {
    CVector2 r(*this);
    r.SetNormalizedVector();
    return r;
  }
	//! set zero vector
	inline void SetZero(){
		x = 0.0;
		y = 0.0;
	}
  double Length() const{
		return sqrt( x*x+y*y );
	}
  double SqLength() const{
		return x*x+y*y;
	}
public:
	double x;	//!< x coordinate value
	double y;	//!< y coordinate value
};

CVector2 operator*(double c, const CVector2& v0);
CVector2 operator*(const CVector2& v0, double c);
double operator * (const CVector2& lhs, const CVector2& rhs);
double operator ^ (const CVector2& lhs, const CVector2& rhs);
CVector2 operator/ (const CVector2& vec, double d); //! divide by real number
CVector2 rotate(const CVector2& p0, double theta);
CVector2 rotate90(const CVector2& p0);

////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////

CVector2 matVec(const double A[4], const CVector2& v);

//! Area of the Triangle
double TriArea(const CVector2& v1,
               const CVector2& v2,
               const CVector2& v3);
double Cross(const CVector2& v1, const CVector2& v2);
double SquareLength(const CVector2& ipo0, const CVector2& ipo1);
double SquareLength(const CVector2& point);
double Length(const CVector2& point);

//! Length between two points
double Distance(const CVector2& ipo0, const CVector2& ipo1);

//! Length between two points
double SquareDistance(const CVector2& ipo0, const CVector2& ipo1);

//! Hight of a triangle : between v1 and line of v2-v3
double TriHeight(const CVector2& v1, const CVector2& v2, const CVector2& v3);

//! compute dot product
inline double Dot(const CVector2& ipo0, const CVector2& ipo1);

// get parameter 't' of the line against point. t=0 is po_s, t=1 is po_e
// this one has implementation in header because GetDist_LineSeg_Point below refers this
double FindNearestPointParameter_Line_Point(const CVector2& po_c,
                                            const CVector2& po_s, const CVector2& po_e);

// get parameter 't' of the line against point. t=0 is po_s, t=1 is po_e
// this one has implementation in header because utility_glvector.h refers this
CVector2 GetNearest_LineSeg_Point(const CVector2& po_c,
                                  const CVector2& po_s, const CVector2& po_e);

// get parameter 't' of the line against point. t=0 is po_s, t=1 is po_e
// this one has implementation in header because utility_glvector.h refers this
double GetDist_LineSeg_Point(const CVector2& po_c,
                             const CVector2& po_s, const CVector2& po_e);
  
bool IsCross_LineSeg_LineSeg(const CVector2& po_s0, const CVector2& po_e0,
                             const CVector2& po_s1, const CVector2& po_e1 );
  
double GetDist_LineSeg_LineSeg(const CVector2& po_s0, const CVector2& po_e0,
                               const CVector2& po_s1, const CVector2& po_e1);
  
//! square root of circumradius
double SquareCircumradius(const CVector2& p0,
                          const CVector2& p1,
                          const CVector2& p2 );

//! center of the circumcircle
bool CenterCircumcircle(const CVector2& p0,
                        const CVector2& p1,
                        const CVector2& p2,
                        CVector2& center);

//! check if Delaunay condition satisfied
// 0 : p3 is inside circum circle on the p0,p1,p2
// 1 :       on         
// 2 :       outsdie 
int DetDelaunay(const CVector2& p0,
                const CVector2& p1,
                const CVector2& p2,
                const CVector2& p3);

// move to paramgeo2d?
CVector2 pointCurve_BezierCubic(double t,
                                const CVector2& p1, const CVector2& p2, const CVector2& p3, const CVector2& p4);

/////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////

//! Area of the Triangle (3 indexes and vertex array)
double TriArea(const int iv1, const int iv2, const int iv3,
               const std::vector<CVector2>& point );

void Polyline_CubicBezierCurve(std::vector<CVector2>& aP,
                               const int n,
                               const std::vector<CVector2>& aCP);

std::vector<CVector2> Polygon_Resample_Polygon(const std::vector<CVector2>& stroke0,
                                               double l);
std::vector<CVector2> Polyline_Resample_Polyline(const std::vector<CVector2>& stroke0,
                                                 double l);
void SecondMomentOfArea_Polygon(CVector2& cg,  double& area,
                               CVector2& pa1, double& I1,
                               CVector2& pa2, double& I2,
                              const std::vector<CVector2>& aVec2D);
double Length_Polygon(const std::vector<CVector2>& aP);

void MeanValueCoordinate2D(double* aW,
                           double px, double py,
                           const double* aXY, int nXY);
void MeanValueCoordinate(std::vector<double>& aW,
                         CVector2& p,
                         std::vector<CVector2>& aVtx);

void makeRandomLoop(unsigned int nCV,
                    std::vector<double>& aCV);

void makeSplineLoop(const std::vector<double>& aCV,
                    std::vector<double>& aVecCurve);

void Translate(std::vector<CVector2>& aP,
               double dx, double dy);
void Rotate(std::vector<CVector2>& aP,
            double dt);
  
//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////

//! 2D bounding box class
class CBoundingBox2D
{
public:
  CBoundingBox2D(){
    x_min=0;	x_max=0;
    y_min=0;	y_max=0;
    isnt_empty = false;
  }
  CBoundingBox2D(double x_min0,double x_max0,  double y_min0,double y_max0)
  : x_min(x_min0),x_max(x_max0),  y_min(y_min0),y_max(y_max0)
  {
    assert( x_min <= x_max );
    assert( y_min <= y_max );
    isnt_empty = true;
  }
  CBoundingBox2D( const CBoundingBox2D& bb )
  : x_min(bb.x_min),x_max(bb.x_max), y_min(bb.y_min),y_max(bb.y_max),
  isnt_empty(bb.isnt_empty){}
  
  CBoundingBox2D& operator+=(const CBoundingBox2D& bb)
  {
    if( !bb.isnt_empty ) return *this;
    if( !isnt_empty ){
      x_max = bb.x_max;	x_min = bb.x_min;
      y_max = bb.y_max;	y_min = bb.y_min;
      this->isnt_empty = bb.isnt_empty;
      return *this;
    }
    x_max = ( x_max > bb.x_max ) ? x_max : bb.x_max;
    x_min = ( x_min < bb.x_min ) ? x_min : bb.x_min;
    y_max = ( y_max > bb.y_max ) ? y_max : bb.y_max;
    y_min = ( y_min < bb.y_min ) ? y_min : bb.y_min;
    return *this;
  }
  bool IsInside(const CVector2& vec)
  {
    if( !isnt_empty ) return false;
    if(   vec.x >= x_min && vec.x <= x_max
       && vec.y >= y_min && vec.y <= y_max ) return true;
    return false;
  }
  bool IsIntersectSphere(const CVector2& vec, const double radius ) const
  {
    if( !isnt_empty ) return false;
    if( vec.x < x_min-radius || vec.x > x_max+radius ||
       vec.y < y_min-radius || vec.y > y_max+radius ) return false;
    return true;
  }
  bool IsIntersect(const CBoundingBox2D& bb_j, double clearance) const
  {
    if( bb_j.x_min > x_max+clearance || bb_j.x_max < x_min-clearance ) return false;
    if( bb_j.y_min > y_max+clearance || bb_j.y_max < y_min-clearance ) return false;
    return true;
  }
public:
  double x_min,x_max,  y_min,y_max;
  bool isnt_empty;	//!< false if there is nothing inside
};

#endif // VECTOR_2D_H


