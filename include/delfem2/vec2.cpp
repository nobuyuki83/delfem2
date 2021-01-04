/*
 * Copyright (c) 2019 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */

#include <cstdlib>
#include "delfem2/vec2.h"

DFM2_INLINE bool delfem2::IsCrossLines(
    const double po_s0[2],
    const double po_e0[2],
    const double po_s1[2],
    const double po_e1[2])
{
  const double area1 = Area_Tri2(po_s0, po_e0, po_s1);
  const double area2 = Area_Tri2(po_s0, po_e0, po_e1);
  if (area1 * area2 > 0.0) return false;
  const double area3 = Area_Tri2(po_s1, po_e1, po_s0);
  const double area4 = Area_Tri2(po_s1, po_e1, po_e0);
  if (area3 * area4 > 0.0) return false;
  return true;
}

// ================================================

template <typename T>
T delfem2::Area_Tri2(
    const T v1[2],
    const T v2[2],
    const T v3[2])
{
  T v0 = ( v2[0] - v1[0] )*( v3[1] - v1[1] ) - ( v3[0] - v1[0] )*( v2[1] - v1[1] );
  return 0.5*v0;
}
#ifndef DFM2_HEADER_ONLY
template float delfem2::Area_Tri2(const float v1[2], const float v2[2], const float v3[2]);
template double delfem2::Area_Tri2(const double v1[2], const double v2[2], const double v3[2]);
#endif

// --------------------------------------

template <typename T>
T delfem2::Dot2(const T w[2], const T v[2])
{
  return w[0]*v[0]+w[1]*v[1];
}
#ifndef DFM2_HEADER_ONLY
template float delfem2::Dot2(const float v1[2], const float v2[2]);
template double delfem2::Dot2(const double v1[2], const double v2[2]);
#endif

// --------------------------------------

// template specialization need to be done in the namespace
namespace delfem2 {
  
template <>
DFM2_INLINE double Length2
 (const double v[2])
{
  return sqrt(v[0]*v[0]+v[1]*v[1]);
}

template <>
DFM2_INLINE float Length2
 (const float v[2])
{
  return sqrtf(v[0]*v[0]+v[1]*v[1]);
}
  
}

// --------------------------------------

// template specialization need to be done in the namespace
namespace delfem2 {

template <>
DFM2_INLINE double Distance2
 (const double v1[2],
  const double v2[2])
{
  return sqrt( (v1[0]-v2[0])*(v1[0]-v2[0]) + (v1[1]-v2[1])*(v1[1]-v2[1]) );
}

template <>
DFM2_INLINE float Distance2
 (const float v1[2],
  const float v2[2])
{
  return sqrtf( (v1[0]-v2[0])*(v1[0]-v2[0]) + (v1[1]-v2[1])*(v1[1]-v2[1]) );
}
  
}

// -------------------------------------

template <typename T>
void delfem2::MatVec2(T w[2], const T A[4], const T v[2])
{
  w[0] = A[0*2+0]*v[0]+A[0*2+1]*v[1];
  w[1] = A[1*2+0]*v[0]+A[1*2+1]*v[1];
}
#ifndef DFM2_HEADER_ONLY
template void delfem2::MatVec2(float w[2], const float A[4], const float v[2]);
template void delfem2::MatVec2(double w[2], const double A[4], const double v[2]);
#endif

// -------------------------------------

template <typename T>
void delfem2::MatMat2(T AB[4], const T A[4], const T B[4])
{
  AB[0*2+0] = A[0*2+0]*B[0*2+0]+A[0*2+1]*B[1*2+0];
  AB[0*2+1] = A[0*2+0]*B[0*2+1]+A[0*2+1]*B[1*2+1];
  AB[1*2+0] = A[1*2+0]*B[0*2+0]+A[1*2+1]*B[1*2+0];
  AB[1*2+1] = A[1*2+0]*B[0*2+1]+A[1*2+1]*B[1*2+1];
}
#ifndef DFM2_HEADER_ONLY
template void delfem2::MatMat2(float AB[4], const float A[4], const float B[4]);
template void delfem2::MatMat2(double AB[4], const double A[4], const double B[4]);
#endif

// --------------------------------------------------

template <typename T>
T delfem2::SquareLength2(const T v[2]){
  return v[0]*v[0]+v[1]*v[1];
}
#ifndef DFM2_HEADER_ONLY
template float delfem2::SquareLength2(const float v[2]);
template double delfem2::SquareLength2(const double v[2]);
#endif

// ------------------------------------------------

template <typename T>
T delfem2::SquareDistance2(const T v1[2], const T v2[2]){
  return (v1[0]-v2[0])*(v1[0]-v2[0]) + (v1[1]-v2[1])*(v1[1]-v2[1]);
}
#ifndef DFM2_HEADER_ONLY
template float delfem2::SquareDistance2(const float v1[2], const float v2[2]);
template double delfem2::SquareDistance2(const double v1[2], const double v2[2]);
#endif

// ------------------------------------------------

// specification of template requires namespace
namespace delfem2 {
  
template <>
DFM2_INLINE void GaussianDistribution2(float noise[2])
{
  float a0 = (float)(rand()/(RAND_MAX+1.0));
  float a1 = (float)(rand()/(RAND_MAX+1.0));
  noise[0] = (float)(sqrtf(-2.f*logf(a0))*cosf(3.1415f*2.f*a1));
  noise[1] = (float)(sqrtf(-2.f*logf(a0))*sinf(3.1415f*2.f*a1));
}

template <>
DFM2_INLINE void GaussianDistribution2(double noise[2])
{
  double a0 = rand()/(RAND_MAX+1.0);
  double a1 = rand()/(RAND_MAX+1.0);
  noise[0] = sqrt(-2.0*log(a0))*cos(3.1415*2*a1);
  noise[1] = sqrt(-2.0*log(a0))*sin(3.1415*2*a1);
}
  
}

// ------------------------------------------------

template <typename T>
void delfem2::Normalize2(T w[2])
{
  double l = Length2(w);
  double invl = 1.0/l;
  w[0] *= invl;
  w[1] *= invl;
}
#ifndef DFM2_HEADER_ONLY
template void delfem2::Normalize2(float w[2]);
template void delfem2::Normalize2(double w[2]);
#endif

// --------------------------------------------------------------

DFM2_INLINE bool delfem2::InverseMat2
 (double invB[4],
  const double B[4])
{
  double det = B[0]*B[3]-B[1]*B[2];
  if (fabs(det)<1.0e-10) return false;
  double invdet = 1.0/det;
  invB[0] = +invdet*B[3];
  invB[1] = -invdet*B[1];
  invB[2] = -invdet*B[2];
  invB[3] = +invdet*B[0];
  return true;
}

DFM2_INLINE void delfem2::gramian2
 (double AtA[3],
  const double A[4])
{
  AtA[0] = A[0*2+0]*A[0*2+0]+A[1*2+0]*A[1*2+0];
  AtA[1] = A[0*2+0]*A[0*2+1]+A[1*2+0]*A[1*2+1];
  AtA[2] = AtA[1];
  AtA[3] = A[0*2+1]*A[0*2+1]+A[1*2+1]*A[1*2+1];
}

DFM2_INLINE void delfem2::VLVt2
 (double A[4],
  double l0,
  double l1,
  const double V[4])
{
  A[0] = l0*V[0]*V[0]+l1*V[1]*V[1];
  A[1] = l0*V[2]*V[0]+l1*V[3]*V[1];
  A[2] = l0*V[0]*V[2]+l1*V[1]*V[3];
  A[3] = l0*V[2]*V[2]+l1*V[3]*V[3];
}


DFM2_INLINE void delfem2::RotationalComponentOfMatrix2
 (double R[4], const double M[4])
{
  const double eps = 1.0e-20;
  double A[4];
  {
    double s = fabs(M[0])+fabs(M[1])+fabs(M[2])+fabs(M[3]);
    if (s<1.0e-10){
      R[0] = 1.0;  R[1] = 0.0; R[2] = 0.0; R[3] = 1.0;
      return;
    }
    double invs = 1.0/s;
    A[0] = invs*M[0];
    A[1] = invs*M[1];
    A[2] = invs*M[2];
    A[3] = invs*M[3];
  }
  double G[4]; gramian2(G, A);
  double l0, l1;
  double v0[2], v1[2];
  {
    double b = G[0]+G[3];
    double c = G[0]*G[3]-G[1]*G[2];
    double d = b*b-4*c;
    if (d<eps){
      l0 = 0.5*b;
      l1 = 0.5*b;
      v0[0] = 0;
      v0[1] = 1;
      v1[0] = 1;
      v1[1] = 0;
    }
    else{
      d = sqrt(d);
      l0 = 0.5*(b+d);
      l1 = 0.5*(b-d);
      v0[0] = G[1];
      v0[1] = G[3]-l1;
      if (SquareLength2(v0)>eps){ Normalize2(v0); }
      v1[0] = G[0]-l0;
      v1[1] = G[2];
      if (SquareLength2(v1)>eps){ Normalize2(v1); }
    }
  }
  double V[4] = { v0[0], v1[0], v0[1], v1[1] };
  if (l0<eps){ l0 = 1; }
  if (l1<eps){ l1 = 1; }
  double il0 = 1.0/sqrt(l0);
  double il1 = 1.0/sqrt(l1);
  double invS[4]; VLVt2(invS, il0, il1, V);
  MatMat2(R, A, invS);
}


// ------------------------------------------------------------------------

// definining operator requires defining namespace
namespace delfem2 {
  
template <typename T>
std::ostream &operator<<(std::ostream &output, const CVec2<T>& v)
{
  output.setf(std::ios::scientific);
  output<<v.p[0]<<" "<<v.p[1];
  return output;
}

template <typename T>
std::istream &operator>>(std::istream &input, CVec2<T>& v)
{
  input>>v.p[0]>>v.p[1];
  return input;
}

template <typename T>
delfem2::CVec2<T> operator*(double c, const CVec2<T>& v0)
{
  return CVec2<T>(v0.p[0]*c, v0.p[1]*c);
}
#ifndef DFM2_HEADER_ONLY
template CVec2d operator*(double c, const CVec2d& v0);
#endif
  
//  ---------------------

template <typename T>
delfem2::CVec2<T> operator*(const CVec2<T>& v0, double c)
{
  return CVec2<T>(v0.p[0]*c, v0.p[1]*c);
}
#ifndef DFM2_HEADER_ONLY
template CVec2d operator*(const CVec2d& v0, double c);
#endif

//  ---------------------

template <typename T>
T operator * (const CVec2<T>& lhs, const CVec2<T>& rhs)
{
  return lhs.p[0]*rhs.p[0] + lhs.p[1]*rhs.p[1];
}
#ifndef DFM2_HEADER_ONLY
template double operator * (const CVec2d& lhs, const CVec2d& rhs);
template float  operator * (const CVec2f& lhs, const CVec2f& rhs);
#endif
  
//  ---------------------

template <typename T>
T operator ^ (const CVec2<T>& lhs, const CVec2<T>& rhs)
{
  return lhs.p[0]*rhs.p[1] - lhs.p[1]*rhs.p[0];
}
#ifndef DFM2_HEADER_ONLY
template double operator ^ (const CVec2d& lhs, const CVec2d& rhs);
template float  operator ^ (const CVec2f& lhs, const CVec2f& rhs);
#endif
  
// -------------------

//! divide by real number
template <typename T>
CVec2<T> operator / (const CVec2<T>& vec, double d)
{
  CVec2<T> temp = vec;
  temp /= d;
  return temp;
}
#ifndef DFM2_HEADER_ONLY
template CVec2d operator / (const CVec2d& vec, double d);
#endif
  
} // namespace delfem2

// -----------------------------------------------------------------

template <typename T>
delfem2::CVec2<T> delfem2::rotate
 (const CVec2<T>& p0, double theta)
{
  CVec2<T> p1;
  double c = cos(theta), s = sin(theta);
  p1.p[0] = c*p0.p[0] - s*p0.p[1];
  p1.p[1] = s*p0.p[0] + c*p0.p[1];
  return p1;
}

template <typename T>
delfem2::CVec2<T> delfem2::rotate90
 (const CVec2<T>& p0)
{
  return CVec2<T>(-p0.p[1], p0.p[0]);
}

// ---------------------------------------------------------------

template <typename T>
delfem2::CVec2<T> delfem2::Mat2Vec(const double A[4], const CVec2<T>& v)
{
  CVec2<T> w;
  MatVec2(w.p, A, v.p);
  return w;
  //  w.p[0] = A[0]*v.p[0]+A[1]*v.p[1];
  //  w.p[1] = A[2]*v.p[0]+A[3]*v.p[1];
}
#ifndef DFM2_HEADER_ONLY
template delfem2::CVec2d delfem2::Mat2Vec(const double A[4], const CVec2d& v);
#endif

// --------------------------------

//! Area of the Triangle
template <typename T>
double delfem2::Area_Tri
(const CVec2<T>& v1,
 const CVec2<T>& v2,
 const CVec2<T>& v3)
{
  return 0.5*( (v2.p[0]-v1.p[0])*(v3.p[1]-v1.p[1]) - (v3.p[0]-v1.p[0])*(v2.p[1]-v1.p[1]) );
}
#ifndef DFM2_HEADER_ONLY
template double delfem2::Area_Tri(const CVec2d& v1, const CVec2d& v2, const CVec2d& v3);
#endif


template <typename T>
double delfem2::Cross(const CVec2<T>& v1, const CVec2<T>& v2){
  return v1.p[0]*v2.p[1] - v2.p[0]*v1.p[1];
}
#ifndef DFM2_HEADER_ONLY
template double delfem2::Cross(const CVec2d& v1, const CVec2d& v2);
#endif


template <typename T>
double delfem2::SquareLength(const CVec2<T>& point)
{
  return  point.p[0]*point.p[0] + point.p[1]*point.p[1];
}
#ifndef DFM2_HEADER_ONLY
template double delfem2::SquareLength(const CVec2d& point);
#endif

// --------------------

template <typename T>
double delfem2::Length(const CVec2<T>& point)
{
  return Length2(point.p);
}
#ifndef DFM2_HEADER_ONLY
template double delfem2::Length(const CVec2d& point);
#endif

// ---------------------

// Distance between two points
template <typename T>
double delfem2::Distance
(const CVec2<T>& ipo0,
 const CVec2<T>& ipo1)
{
  return Distance2(ipo0.p, ipo1.p);
}
#ifndef DFM2_HEADER_ONLY
template double delfem2::Distance(const CVec2d& ipo0, const CVec2d& ipo1);
#endif

// ---------------------

// Distance between two points
template <typename T>
double delfem2::SquareDistance
(const CVec2<T>& ipo0,
 const CVec2<T>& ipo1)
{
  return delfem2::SquareDistance2(ipo0.p, ipo1.p);
}

// Hight of a triangle : between v1 and line of v2-v3
template <typename T>
double delfem2::TriHeight(const CVec2<T>& v1, const CVec2<T>& v2, const CVec2<T>& v3){
  const double area = Area_Tri(v1,v2,v3);
  const double len = sqrt( SquareDistance(v2,v3) );
  return area*2.0/len;
}

// compute dot product
template <typename T>
double delfem2::Dot(const CVec2<T>& ipo0, const CVec2<T>& ipo1)
{
  return Dot2(ipo0.p, ipo1.p);
}
#ifndef DFM2_HEADER_ONLY
template double delfem2::Dot(const CVec2d& ipo0, const CVec2d& ipo1);
#endif

// get parameter 't' of the line against point. t=0 is po_s, t=1 is po_e
// this one has implementation in header because GetDist_LineSeg_Point below refers this
template <typename T>
double delfem2::FindNearestPointParameter_Line_Point
(const CVec2<T>& po_c,
 const CVec2<T>& po_s, const CVec2<T>& po_e)
{
  const CVec2<T>& es = po_e-po_s;
  const CVec2<T>& sc = po_s-po_c;
  const double a = SquareLength(es);
  const double b = Dot(es, sc);
  return -b/a;
}

// --------------------------------------

template <typename T>
delfem2::CVec2<T> delfem2::GetNearest_LineSeg_Point(
    const CVec2<T>& po_c,
    const CVec2<T>& po_s,
    const CVec2<T>& po_e)
{
  double t = FindNearestPointParameter_Line_Point(po_c, po_s, po_e);
  if (t < 0){ return po_s; }
  if (t > 1){ return po_e; }
  return po_s+t*(po_e-po_s);
}
#ifndef DFM2_HEADER_ONLY
template delfem2::CVec2d delfem2::GetNearest_LineSeg_Point(
    const CVec2d& po_c, const CVec2d& po_s, const CVec2d& po_e);
#endif

// --------------------------------------

template <typename T>
double delfem2::GetDist_LineSeg_Point
(const CVec2<T>& po_c,
 const CVec2<T>& po_s, const CVec2<T>& po_e)
{
  CVec2<T> p = GetNearest_LineSeg_Point(po_c, po_s,po_e);
  return Distance(p, po_c);
}
#ifndef DFM2_HEADER_ONLY
template double delfem2::GetDist_LineSeg_Point(
    const CVec2d& po_c, const CVec2d& po_s, const CVec2d& po_e);
#endif


template <typename T>
bool delfem2::IsCross_LineSeg_LineSeg
(const CVec2<T>& po_s0, const CVec2<T>& po_e0,
 const CVec2<T>& po_s1, const CVec2<T>& po_e1 )
{
  {
    const double min0x = ( po_s0.p[0] < po_e0.p[0] ) ? po_s0.p[0] : po_e0.p[0];
    const double max0x = ( po_s0.p[0] > po_e0.p[0] ) ? po_s0.p[0] : po_e0.p[0];
    const double max1x = ( po_s1.p[0] > po_e1.p[0] ) ? po_s1.p[0] : po_e1.p[0];
    const double min1x = ( po_s1.p[0] < po_e1.p[0] ) ? po_s1.p[0] : po_e1.p[0];
    const double min0y = ( po_s0.p[1] < po_e0.p[1] ) ? po_s0.p[1] : po_e0.p[1];
    const double max0y = ( po_s0.p[1] > po_e0.p[1] ) ? po_s0.p[1] : po_e0.p[1];
    const double max1y = ( po_s1.p[1] > po_e1.p[1] ) ? po_s1.p[1] : po_e1.p[1];
    const double min1y = ( po_s1.p[1] < po_e1.p[1] ) ? po_s1.p[1] : po_e1.p[1];
    const double len = ((max0x-min0x)+(max0y-min0y)+(max1x-min1x)+(max1y-min1y))*0.0001;
    //		std::cout << len << std::endl;
    if( max1x+len < min0x ) return false;
    if( max0x+len < min1x ) return false;
    if( max1y+len < min0y ) return false;
    if( max0y+len < min1y ) return false;
  }
  const double area1 = Area_Tri(po_s0,po_e0,po_s1);
  const double area2 = Area_Tri(po_s0,po_e0,po_e1);
  const double area3 = Area_Tri(po_s1,po_e1,po_s0);
  const double area4 = Area_Tri(po_s1,po_e1,po_e0);  
  //	std::cout << area1 << " " << area2 << " " << area3 << " " << area4 << std::endl;
  const double a12 = area1*area2; if( a12 > 0 ) return false;
  const double a34 = area3*area4; if( a34 > 0 ) return false;
  return true;
}

template <typename T>
double delfem2::GetDist_LineSeg_LineSeg
(const CVec2<T>& po_s0, const CVec2<T>& po_e0,
 const CVec2<T>& po_s1, const CVec2<T>& po_e1)
{
  if( IsCross_LineSeg_LineSeg(po_s0,po_e0, po_s1,po_e1) ) return -1;
  const double ds1 = GetDist_LineSeg_Point(po_s0,po_s1,po_e1);
  const double de1 = GetDist_LineSeg_Point(po_e0,po_s1,po_e1);
  const double ds0 = GetDist_LineSeg_Point(po_s1,po_s0,po_e0);
  const double de0 = GetDist_LineSeg_Point(po_e1,po_s0,po_e0);
  double min_dist = ds1;
  min_dist = ( de1 < min_dist ) ? de1 : min_dist;
  min_dist = ( ds0 < min_dist ) ? ds0 : min_dist;
  min_dist = ( de0 < min_dist ) ? de0 : min_dist;    
  return min_dist;  
}

// square root of circumradius
template <typename T>
double delfem2::SquareCircumradius
(const CVec2<T>& p0,
 const CVec2<T>& p1,
 const CVec2<T>& p2 )
{
	const double area = Area_Tri(p0,p1,p2);
  
	const double dtmp0 = SquareDistance(p1,p2);
	const double dtmp1 = SquareDistance(p0,p2);
	const double dtmp2 = SquareDistance(p0,p1);
  
	return dtmp0*dtmp1*dtmp2/(16.0*area*area);
}

//! center of the circumcircle
template <typename T>
bool delfem2::CenterCircumcircle
(const CVec2<T>& p0,
 const CVec2<T>& p1,
 const CVec2<T>& p2,
 CVec2<T>& center)
{
  const double area = Area_Tri(p0,p1,p2);
  if( fabs(area) < 1.0e-10 ){ return false; }
  const double tmp_val = 1.0/(area*area*16.0);
  
  const double dtmp0 = SquareDistance(p1,p2);
  const double dtmp1 = SquareDistance(p0,p2);
  const double dtmp2 = SquareDistance(p0,p1);
  
  const double etmp0 = tmp_val*dtmp0*(dtmp1+dtmp2-dtmp0);
  const double etmp1 = tmp_val*dtmp1*(dtmp0+dtmp2-dtmp1);
  const double etmp2 = tmp_val*dtmp2*(dtmp0+dtmp1-dtmp2);
  
  center.p[0] = etmp0*p0.p[0] + etmp1*p1.p[0] + etmp2*p2.p[0];
  center.p[1] = etmp0*p0.p[1] + etmp1*p1.p[1] + etmp2*p2.p[1];
  return true;
}


// check if Delaunay condition satisfied
// 0 : p3 is inside circum circle on the p0,p1,p2
// 1 :       on
// 2 :       outsdie
template <typename T>
int delfem2::DetDelaunay
(const CVec2<T>& p0,
 const CVec2<T>& p1,
 const CVec2<T>& p2,
 const CVec2<T>& p3)
{
	const double area = Area_Tri(p0,p1,p2);
	if( fabs(area) < 1.0e-10 ){
		return 3;
	}
	const double tmp_val = 1.0/(area*area*16.0);
  
	const double dtmp0 = SquareDistance(p1,p2);
	const double dtmp1 = SquareDistance(p0,p2);
	const double dtmp2 = SquareDistance(p0,p1);
  
	const double etmp0 = tmp_val*dtmp0*(dtmp1+dtmp2-dtmp0);
	const double etmp1 = tmp_val*dtmp1*(dtmp0+dtmp2-dtmp1);
	const double etmp2 = tmp_val*dtmp2*(dtmp0+dtmp1-dtmp2);
  
	const CVec2<T> out_center(etmp0*p0.p[0] + etmp1*p1.p[0] + etmp2*p2.p[0],
                            etmp0*p0.p[1] + etmp1*p1.p[1] + etmp2*p2.p[1] );
  
	const double qradius = SquareDistance(out_center,p0);
	const double qdistance = SquareDistance(out_center,p3);
  
  //	assert( fabs( qradius - SquareLength(out_center,p1) ) < 1.0e-10*qradius );
  //	assert( fabs( qradius - SquareLength(out_center,p2) ) < 1.0e-10*qradius );
  
	const double tol = 1.0e-20;
	if( qdistance > qradius*(1.0+tol) ){ return 2; }	// outside the circumcircle
	else{
		if( qdistance < qradius*(1.0-tol) ){ return 0; }	// inside the circumcircle
		else{ return 1;	}	// on the circumcircle
	}
	return 0;
}
#ifndef DFM2_HEADER_ONLY
template int delfem2::DetDelaunay(const CVec2d& p0, const CVec2d& p1, const CVec2d& p2, const CVec2d& p3);
#endif



template <typename T>
delfem2::CVec2<T> delfem2::pointCurve_BezierCubic
(double t,
 const CVec2<T>& p1, const CVec2<T>& p2, const CVec2<T>& p3, const CVec2<T>& p4)
{
  double tp = 1.0-t;
  return t*t*t*p4 + 3*t*t*tp*p3 + 3*t*tp*tp*p2 + tp*tp*tp*p1;
}
#ifndef DFM2_HEADER_ONLY
template delfem2::CVec2d delfem2::pointCurve_BezierCubic(double t,
                                                   const CVec2d& p1,
                                                   const CVec2d& p2,
                                                   const CVec2d& p3,
                                                   const CVec2d& p4);
#endif

template <typename T>
delfem2::CVec2<T> delfem2::pointCurve_BezierQuadratic
 (double t,
  const CVec2<T>& p1, const CVec2<T>& p2, const CVec2<T>& p3)
{
  double tp = 1.0-t;
  return (t*t)*p3 + (2*t*tp)*p2 + (tp*tp)*p1;
}
#ifndef DFM2_HEADER_ONLY
template delfem2::CVec2d delfem2::pointCurve_BezierQuadratic(double t,
                                                       const CVec2d& p1,
                                                       const CVec2d& p2,
                                                       const CVec2d& p3);
#endif


// ----------------------------------------------------------------------------------
// std::vector starts from here

//! Area of the Triangle (3 indexes and vertex array)
template <typename T>
double delfem2::Area_Tri
(const int iv1, const int iv2, const int iv3,
 const std::vector<CVec2<T>>& point )
{
  return Area_Tri(point[iv1],point[iv2],point[iv3]);
}

template <typename T>
void delfem2::Polyline_CubicBezierCurve
(std::vector<CVec2<T>>& aP,
 const int n, 
 const std::vector<CVec2<T>>& aCP)
{
  int ns = (int)(aCP.size()/3);
  aP.resize(ns*n+1);
  for(int is=0;is<ns;is++){
    for(int i=0;i<n;i++){
      double t = (double)i/n;
      aP[is*n+i] = pointCurve_BezierCubic(t,aCP[is*3+0],aCP[is*3+1],aCP[is*3+2],aCP[is*3+3]);
    }
  }
  aP[ns*n] = aCP[ns*3];
}
#ifndef DFM2_HEADER_ONLY
template void delfem2::Polyline_CubicBezierCurve(std::vector<CVec2d>& aP, const int n, const std::vector<CVec2d>& aCP);
#endif


// ------------------------------

template <typename T>
void delfem2::Polyline_BezierCubic
(std::vector<CVec2<T>>& aP,
 const unsigned int n,
 const CVec2<T>& p1, const CVec2<T>& p2, const CVec2<T>& p3, const CVec2<T>& p4)
{
  aP.resize(n);
  for(unsigned int i=0;i<n;++i){
    const double t = (double)i/(n-1.0);
    aP[i] = pointCurve_BezierCubic(t,
                                   p1, p2, p3, p4);
  }
}
template void delfem2::Polyline_BezierCubic
 (std::vector<CVec2d>& aP,
  const unsigned int n,
  const CVec2d& p1, const CVec2d& p2, const CVec2d& p3, const CVec2d& p4);

// --------------

template <typename T>
void delfem2::Polyline_BezierQuadratic(
    std::vector<CVec2<T>>& aP,
    const unsigned int n,
    const CVec2<T>& p1,
    const CVec2<T>& p2,
    const CVec2<T>& p3)
{
  aP.resize(n);
  for(unsigned int i=0;i<n;++i){
    const double t = (double)i/(n-1.0);
    aP[i] = pointCurve_BezierQuadratic(t,
        p1, p2, p3);
  }
}
#ifndef DFM2_HEADER_ONLY
template void delfem2::Polyline_BezierQuadratic(
    std::vector<CVec2d>& aP,
    const unsigned int n,
    const CVec2d& p1,
    const CVec2d& p2,
    const CVec2d& p3);
#endif

// ---------------------------------------

template <typename T>
void delfem2::Translate(
    std::vector<CVec2<T>>& aP,
    double dx, double dy)
{
  for(auto & ip : aP){
    ip.p[0] += dx;
    ip.p[1] += dy;
  }
}
#ifndef DFM2_HEADER_ONLY
template void delfem2::Translate(std::vector<CVec2d>& aP, double dx, double dy);
#endif

// ----------------------------

namespace delfem2 {

template<>
DFM2_INLINE void Rotate(std::vector<CVec2d> &aP, double dt) {
  for (auto &ip : aP) {
    double x0 = ip.p[0];
    double y0 = ip.p[1];
    ip.p[0] = cos(dt) * x0 - sin(dt) * y0;
    ip.p[1] = sin(dt) * x0 + cos(dt) * y0;
  }
}

}

// ----------------------------

template <typename T>
std::vector<delfem2::CVec2<T>> delfem2::Polyline_Resample_Polyline
(const std::vector<CVec2<T>>& stroke0,
 double l)
{
  if( stroke0.empty() ){
    std::vector<CVec2<T>> a;
    return a;
  }
  std::vector<CVec2<T>> stroke;
  stroke.push_back( stroke0[0] );
  int jcur = 0;
  double rcur = 0;
  double lcur = l;
  for(;;){
    if( jcur >= (int)stroke0.size()-1 ) break;
    double lenj  = (stroke0[jcur+1]-stroke0[jcur]).Length();
    double lenjr = lenj*(1.0-rcur);
    if( lenjr > lcur ){ // put point in this segment
      rcur += lcur/lenj;
      stroke.push_back( (1-rcur)*stroke0[jcur] + rcur*stroke0[jcur+1] );
      lcur = l;
    }
    else{ // next segment
      lcur -= lenjr;
      rcur =0;
      jcur++;
    }
  }
  //  stroke.push_back( stroke0.back() );
  return stroke;
}
#ifndef DFM2_HEADER_ONLY
template std::vector<delfem2::CVec2d> delfem2::Polyline_Resample_Polyline(
    const std::vector<CVec2d>& stroke0, double l);
#endif

// -----------------------------------------------------------

template <typename T>
void delfem2::MakeMassMatrixTri(
    double M[9],
    double rho,
    const unsigned int aIP[3],
    const std::vector<CVec2<T>>& aVec2)
{
  assert( aIP[0]<aVec2.size() );
  assert( aIP[1]<aVec2.size() );
  assert( aIP[2]<aVec2.size() );
  const double P[3][2] = {
    {aVec2[aIP[0]].p[0],aVec2[aIP[0]].p[1]},
    {aVec2[aIP[1]].p[0],aVec2[aIP[1]].p[1]},
    {aVec2[aIP[2]].p[0],aVec2[aIP[2]].p[1]} };
  const double Area = delfem2::Area_Tri2(P[0], P[1], P[2]);
  {
    const double tmp = rho*Area/3.0;
    M[0] = M[4] = M[8] = tmp;
    M[1] = M[2] = M[3] = M[5] = M[6] = M[7] = 0.0;
  }
  {
    const double tmp = rho*Area/12.0;
    M[0] = M[4] = M[8] = tmp*2.0;
    M[1] = M[2] = M[3] = M[5] = M[6] = M[7] = tmp;
  }
}
