/*
 * Copyright (c) 2019 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */

#ifndef DFM2_VEC2_H
#define DFM2_VEC2_H

#include <cstdlib>
#include <cassert>
#include <iostream>
#include <cmath>
#include <vector>

// -----------------------------------------------------

namespace delfem2 {


template <typename T>
T Dot2(const T w[2], const T v[2]);

template <typename T>
T Length2(const T v[2]);

template <typename T>
T Distance2(const T v1[2], const T v2[2]);

template <typename T>
void MatVec2(T w[2], const T A[4], const T v[2]);

template <typename T>
void MatMat2(T AB[4], const T A[4], const T B[4]);

template <typename T>
T SquareLength2(const T v[2]);

template <typename T>
T SquareDistance2(const T v1[2], const T v2[2]);

template <typename T>
void Normalize2(T w[2]);

// ----------------------------------------------------
  
template <typename T>
T Area_Tri2(const T v1[2], const T v2[2], const T v3[2]);
  
template <typename T>
void GaussianDistribution2(T noise[2]);
  
void MeanValueCoordinate2D(double *aW,
                           double px, double py,
                           const double *aXY, unsigned int nXY);
  
} // end namespace delfem2

bool InverseMat2(double invB[4], const double B[4]);
void gramian2(double AtA[3], const double A[4]);
void VLVt2(double A[4], double l0, double l1, const double V[4]);
void RotationalComponentOfMatrix2(double R[4], const double M[4]);

// -----------------------------------------------------

namespace delfem2 {

template <typename T>
class CVec2;

template <typename T>
CVec2<T> operator*(double, const CVec2<T>&);
  
template <typename T>
CVec2<T> operator*(const CVec2<T>&, double);
  
template <typename T>
std::ostream &operator<<(std::ostream &output, const CVec2<T>& v);
  
template <typename T>
std::istream &operator>>(std::istream &input, CVec2<T>& v);
  
template <typename T>
CVec2<T> operator*(double c, const CVec2<T>& v0);

template <typename T>
CVec2<T> operator*(const CVec2<T>& v0, double c);

template <typename T>
double operator * (const CVec2<T>& lhs, const CVec2<T>& rhs);

template <typename T>
double operator ^ (const CVec2<T>& lhs, const CVec2<T>& rhs);

template <typename T>
CVec2<T> operator/ (const CVec2<T>& vec, double d); //! divide by real number

/**
 * @brief 2 dimensional vector class
 * @todo use template for this class
 */
template <typename T>
class CVec2{
public:
  CVec2(): p{0.0, 0.0} {}
	CVec2( const CVec2& rhs ){
		this->p[0] = rhs.p[0];
		this->p[1] = rhs.p[1];
	}
	CVec2(double x, double y){
		this->p[0] = x;
		this->p[1] = y;
	}	
  
  CVec2 operator-() const{
    return CVec2(-p[0],-p[1]);
  }
	inline CVec2& operator+=(const CVec2& rhs){
		p[0] += rhs.p[0];
		p[1] += rhs.p[1];
		return *this;
	}
	inline CVec2& operator-=(const CVec2& rhs){
		p[0] -= rhs.p[0];
		p[1] -= rhs.p[1];
		return *this;
	}
	inline CVec2& operator*=(double scale){
		p[0] *= scale;
		p[1] *= scale;
		return *this;
  }
  inline CVec2& operator/=(double d){
    if (fabs(d) < 1.0e-6){ assert(0); return *this; }
    p[0] /= d; p[1] /= d;
    return *this;
  }
	inline CVec2 operator+(const CVec2& rhs) const {
		CVec2 v = *this;
		return v += rhs;
	}
	inline CVec2 operator-(const CVec2& rhs) const {
		CVec2 v = *this;
		return v -= rhs;
	}
  inline double operator[](int i) const{
    if (i==0) return p[0];
    if (i==1) return p[1];
    return 0;
  }
  inline double& operator[](int i){
    if (i==0) return p[0];
    if (i==1) return p[1];
    assert(0);
    return p[0];
  }
	//! @brief normalize length
	inline void SetNormalizedVector(){
		const double mag = Length();
		p[0] /= mag;
		p[1] /= mag;
	}
  CVec2 Normalize() const {
    CVec2 r(*this);
    r.SetNormalizedVector();
    return r;
  }
	//! set zero vector
	inline void SetZero(){
		p[0] = 0.0;
		p[1] = 0.0;
	}
  inline double x() const { return p[0]; }
  inline double y() const { return p[1]; }
  double Length() const{
		return sqrt( p[0]*p[0]+p[1]*p[1] );
	}
  double SqLength() const{
		return p[0]*p[0]+p[1]*p[1];
	}
public:
  double p[2];
};
using CVec2d = CVec2<double>;
using CVec2f = CVec2<float>;
  
template <typename T>
CVec2<T> rotate(const CVec2<T>& p0, double theta);

template <typename T>
CVec2<T> rotate90(const CVec2<T>& p0);

// --------------------------------------------

template <typename T>
CVec2<T> Mat2Vec(const double A[4], const CVec2<T>& v);

//! @brief Area of the Triangle
template <typename T>
double Area_Tri(const CVec2<T>& v1,
               const CVec2<T>& v2,
               const CVec2<T>& v3);

template <typename T>
double Cross(const CVec2<T>& v1, const CVec2<T>& v2);
  
template <typename T>
double SquareDistance(const CVec2<T>& ipo0, const CVec2<T>& ipo1);

template <typename T>
double SquareLength(const CVec2<T>& point);
  
template <typename T>
double Length(const CVec2<T>& point);

//! @brief Length between two points
template <typename T>
double Distance(const CVec2<T>& ipo0, const CVec2<T>& ipo1);

//! @brief Length between two points
template <typename T>
double SquareDistance(const CVec2<T>& ipo0, const CVec2<T>& ipo1);

//! @brief Hight of a triangle : between v1 and line of v2-v3
template <typename T>
double TriHeight(const CVec2<T>& v1, const CVec2<T>& v2, const CVec2<T>& v3);

//! @brief compute dot product
template <typename T>
double Dot(const CVec2<T>& ipo0, const CVec2<T>& ipo1);

//! @details get parameter 't' of the line against point. t=0 is po_s, t=1 is po_e
template <typename T>
double FindNearestPointParameter_Line_Point(const CVec2<T>& po_c,
                                            const CVec2<T>& po_s, const CVec2<T>& po_e);

//! @details  get parameter 't' of the line against point. t=0 is po_s, t=1 is po_e
template <typename T>
CVec2<T> GetNearest_LineSeg_Point(const CVec2<T>& po_c,
                                  const CVec2<T>& po_s, const CVec2<T>& po_e);

//! @details  get parameter 't' of the line against point. t=0 is po_s, t=1 is po_e
template <typename T>
double GetDist_LineSeg_Point(const CVec2<T>& po_c,
                             const CVec2<T>& po_s, const CVec2<T>& po_e);
  
template <typename T>
bool IsCross_LineSeg_LineSeg(const CVec2<T>& po_s0, const CVec2<T>& po_e0,
                             const CVec2<T>& po_s1, const CVec2<T>& po_e1 );
  
template <typename T>
double GetDist_LineSeg_LineSeg(const CVec2<T>& po_s0, const CVec2<T>& po_e0,
                               const CVec2<T>& po_s1, const CVec2<T>& po_e1);
  
/**
 * @brief square root of circumradius
 */
template <typename T>
double SquareCircumradius(const CVec2<T>& p0,
                          const CVec2<T>& p1,
                          const CVec2<T>& p2 );

/**
 * @brief center of the circumcircle
 */
template <typename T>
bool CenterCircumcircle(const CVec2<T>& p0,
                        const CVec2<T>& p1,
                        const CVec2<T>& p2,
                        CVec2<T>& center);

/**
 * @brief check if Delaunay condition satisfied
 * @return
 * 0 : p3 is inside circum circle on the p0,p1,p2
 * 1 :       on
 * 2 :       outsdie
 */
template <typename T>
int DetDelaunay(const CVec2<T>& p0,
                const CVec2<T>& p1,
                const CVec2<T>& p2,
                const CVec2<T>& p3);

// move to paramgeo2d?
template <typename T>
CVec2<T> pointCurve_BezierCubic
(double t,
 const CVec2<T>& p1, const CVec2<T>& p2, const CVec2<T>& p3, const CVec2<T>& p4);

template <typename T>
CVec2<T> pointCurve_BezierQuadratic
(double t,
 const CVec2<T>& p1, const CVec2<T>& p2, const CVec2<T>& p3);


// ---------------------------------------------------------------


//! @brief translate all the points
template <typename T>
void Translate(std::vector<CVec2<T> >& aP,
               double dx, double dy);

template <typename T>
void Rotate(std::vector<CVec2<T> >& aP,
            double dt);

//! @brief Area of the Triangle (3 indexes and vertex array)
template <typename T>
double Area_Tri(const int iv1, const int iv2, const int iv3,
               const std::vector<CVec2<T> >& point );

template <typename T>
void Polyline_CubicBezierCurve(std::vector<CVec2<T> >& aP,
                               const int n,
                               const std::vector<CVec2<T> >& aCP);

template <typename T>
void Polyline_BezierCubic(std::vector<CVec2<T> >& aP,
                          const unsigned int n,
                          const CVec2<T>& p1, const CVec2<T>& p2,
                          const CVec2<T>& p3, const CVec2<T>& p4);

template <typename T>
void Polyline_BezierQuadratic(std::vector<CVec2<T>>& aP,
                              const unsigned int n,
                              const CVec2<T>& p1, const CVec2<T>& p2, const CVec2<T>& p3);

template <typename T>
std::vector<CVec2<T> > Polygon_Resample_Polygon(const std::vector<CVec2<T> >& stroke0,
                                               double l);
  
template <typename T>
std::vector<CVec2<T> > Polyline_Resample_Polyline(const std::vector<CVec2<T> >& stroke0,
                                                 double l);
  
template <typename T>
void SecondMomentOfArea_Polygon(CVec2<T>& cg,  double& area,
                               CVec2<T>& pa1, double& I1,
                               CVec2<T>& pa2, double& I2,
                              const std::vector<CVec2<T> >& aVec2D);
  
template <typename T>
double Length_Polygon(const std::vector<CVec2<T> >& aP);

template <typename T>
double Area_Polygon(const std::vector<CVec2<T> >& aP);

template <typename T>
void MeanValueCoordinate(std::vector<double>& aW,
                         CVec2<T>& p,
                         std::vector<CVec2<T> >& aVtx);

template <typename T>
void makeRandomLoop(unsigned int nCV,
                    std::vector<double>& aCV);

template <typename T>
void makeSplineLoop(const std::vector<double>& aCV,
                    std::vector<double>& aVecCurve);

template <typename T>
void FixLoopOrientation(std::vector<int>& loopIP,
                        const std::vector<int>& loopIP_ind,
                        const std::vector<CVec2<T> >& aXY);

template <typename T>
std::vector<CVec2<T> > Polygon_Invert(const std::vector<CVec2<T> >& aP);

template <typename T>
std::vector<double> XY_Polygon(const std::vector<CVec2<T> >& aP);

template <typename T>
void ResamplingLoop(std::vector<int>& loopIP1_ind,
                    std::vector<int>& loopIP1,
                    std::vector<CVec2<T> >& aXY,
                    double max_edge_length);

template <typename T>
void JArray_FromVecVec_XY(std::vector<int>& aIndXYs,
                          std::vector<int>& loopIP0,
                          std::vector<CVec2<T> >& aXY,
                          const std::vector< std::vector<double> >& aaXY);

template <typename T>
void MakeMassMatrixTri(double M[9],
                       double rho,
                       const unsigned int aIP[3],
                       const std::vector<CVec2<T> >& aVec2);
  
template <typename T>
bool IsInclude_Loop(const double co[],
                    const int ixy_stt, const int ixy_end,
                    const std::vector<CVec2<T> >& aXY);

template <typename T>
bool CheckInputBoundaryForTriangulation (const std::vector<int>& loopIP_ind,
                                         const std::vector<CVec2<T> >& aXY);
  

// -------------------------------------------------------------

/**
 * @brief 2D bounding box
 * @details inactive if x_min > x_max
 */
class CBoundingBox2D
{
public:
  CBoundingBox2D(){
    x_min=1;	x_max=-1;
    y_min=0;	y_max=0;
  }
  CBoundingBox2D(double x_min0,double x_max0,  double y_min0,double y_max0)
  : x_min(x_min0),x_max(x_max0),  y_min(y_min0),y_max(y_max0)
  {}
  CBoundingBox2D( const CBoundingBox2D& bb )
  : x_min(bb.x_min),x_max(bb.x_max), y_min(bb.y_min),y_max(bb.y_max) {}
  
  // -------------------------
  // const functions from here
  
  bool isActive() const { return x_min <= x_max; }
  bool IsIntersectSphere(const CVec2<double>& vec, const double radius ) const
  {
    if( !isActive() ) return false;
    if( vec.p[0] < x_min-radius || vec.p[0] > x_max+radius ||
       vec.p[1] < y_min-radius || vec.p[1] > y_max+radius ) return false;
    return true;
  }
  bool IsIntersect(const CBoundingBox2D& bb_j, double clearance) const
  {
    if( bb_j.x_min > x_max+clearance || bb_j.x_max < x_min-clearance ) return false;
    if( bb_j.y_min > y_max+clearance || bb_j.y_max < y_min-clearance ) return false;
    return true;
  }
  std::vector<double> MinMaxXYZ() const {
    const double tmp[6] = {x_min,x_max, y_min,y_max, 0.0, 0.0};
    std::vector<double> bb(tmp,tmp+6);
    return bb;
  }
  double LengthDiagonal() const {
    return sqrt( (x_max-x_min)*(x_max-x_min) + (y_max-y_min)*(y_max-y_min) );
  }
  
  // -------------------------------
  // non const functions from here
  
  void Add(double x0, double y0){
    if( !isActive() ){
      x_min = x_max = x0;
      y_min = y_max = y0;
      return;
    }
    x_max = ( x_max > x0 ) ? x_max : x0;
    x_min = ( x_min < x0 ) ? x_min : x0;
    y_max = ( y_max > y0 ) ? y_max : y0;
    y_min = ( y_min < y0 ) ? y_min : y0;
  }
  CBoundingBox2D& operator+=(const CBoundingBox2D& bb)
  {
    if( !bb.isActive() ){ return *this; }
    if( !isActive() ){
      x_max = bb.x_max;	x_min = bb.x_min;
      y_max = bb.y_max;	y_min = bb.y_min;
      return *this;
    }
    x_max = ( x_max > bb.x_max ) ? x_max : bb.x_max;
    x_min = ( x_min < bb.x_min ) ? x_min : bb.x_min;
    y_max = ( y_max > bb.y_max ) ? y_max : bb.y_max;
    y_min = ( y_min < bb.y_min ) ? y_min : bb.y_min;
    return *this;
  }
  bool IsInside(const CVec2<double>& vec)
  {
    if( !isActive() ) return false;
    if(   vec.p[0] >= x_min && vec.p[0] <= x_max
       && vec.p[1] >= y_min && vec.p[1] <= y_max ) return true;
    return false;
  }
public:
  double x_min,x_max,  y_min,y_max;
};
  
  
}

#endif // VEC_2


