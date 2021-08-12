/*
 * Copyright (c) 2019 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */

#include <cstdlib>
#include <cmath>
#include <stack>

#include "delfem2/geosolidelm_v3.h"

#ifndef M_PI
#  define M_PI 3.14159265358979323846
#endif

// =====================================
// below: unexposed 

namespace delfem2 {
namespace solidelm {

DFM2_INLINE bool MyIsnan(double x) { return x!=x; }

// evaluate cubic function
template <typename REAL>
DFM2_INLINE REAL EvaluateCubic(
    REAL x,
    REAL k0, REAL k1, REAL k2, REAL k3) // coefficient of cubic function
{
  return k0 + k1*x + k2*x*x + k3*x*x*x;
}
#ifdef DFM2_STATIC_LIBRARY
template float EvaluateCubic(float r2, float k0, float k1, float k2, float k3);
template double EvaluateCubic(double r2, double k0, double k1, double k2, double k3);
#endif


// find root of cubic function using bisection method
DFM2_INLINE double FindRootCubic_Bisect(
    double r0, double r1,
    double v0, double v1,
    double k0, double k1, double k2, double k3)
{
  assert( v0*v1 <= 0 );
  if( v0*v1 == 0 ){
    if( v0 == 0 ){ return r0; }
    else{ return r1; }
  }
  for(unsigned int itr=0;itr<15;itr++){
    const double r2 = 0.5*(r0+r1);
    const double v2 = EvaluateCubic(r2, k0,k1,k2,k3);
    if( v2 == 0 ){ return r2; }
    if( v0*v2 < 0 ){ r1 = r2; v1 = v2; }
    else{            r0 = r2; v0 = v2; }
  }
  return 0.5*(r0+r1);
}


// there is another impelemntation in quat.h so this is "static function"
// transform vector with quaternion
template <typename REAL>
DFM2_INLINE void MyQuatVec(
  REAL vo[],
  const REAL q[],
  const REAL vi[])
{
  REAL x2 = q[1] * q[1] * 2;
  REAL y2 = q[2] * q[2] * 2;
  REAL z2 = q[3] * q[3] * 2;
  REAL xy = q[1] * q[2] * 2;
  REAL yz = q[2] * q[3] * 2;
  REAL zx = q[3] * q[1] * 2;
  REAL xw = q[1] * q[0] * 2;
  REAL yw = q[2] * q[0] * 2;
  REAL zw = q[3] * q[0] * 2;
  vo[0] = (1 - y2 - z2)*vi[0] + (xy - zw    )*vi[1] + (zx + yw    )*vi[2];
  vo[1] = (xy + zw    )*vi[0] + (1 - z2 - x2)*vi[1] + (yz - xw    )*vi[2];
  vo[2] = (zx - yw    )*vi[0] + (yz + xw    )*vi[1] + (1 - x2 - y2)*vi[2];
}
#ifdef DFM2_STATIC_LIBRARY
template void MyQuatVec(float vo[], const float q[], const float vi[]);
template void MyQuatVec(double vo[], const double q[], const double vi[]);
#endif

// ----------------------

template <typename REAL>
DFM2_INLINE void MyMat4Vec3
 (REAL vo[3],
  const REAL M[16], const REAL vi[3])
{
  vo[0] = M[0*4+0]*vi[0] + M[0*4+1]*vi[1] + M[0*4+2]*vi[2];
  vo[1] = M[1*4+0]*vi[0] + M[1*4+1]*vi[1] + M[1*4+2]*vi[2];
  vo[2] = M[2*4+0]*vi[0] + M[2*4+1]*vi[1] + M[2*4+2]*vi[2];
}
#ifdef DFM2_STATIC_LIBRARY
template void MyMat4Vec3(float vo[3],
    const float M[16], const float vi[3]);
template void MyMat4Vec3(double vo[3],
    const double M[16], const double vi[3]);
#endif

// ----------------------

// there is formal implementation in quat.cpp so this is static to avoid dumplicated
template <typename REAL>
DFM2_INLINE void MyQuatConjVec(
    REAL vo[3],
    const REAL q[4],
    const REAL vi[3])
{
  REAL x2 = q[1] * q[1] * 2;
  REAL y2 = q[2] * q[2] * 2;
  REAL z2 = q[3] * q[3] * 2;
  REAL xy = q[1] * q[2] * 2;
  REAL yz = q[2] * q[3] * 2;
  REAL zx = q[3] * q[1] * 2;
  REAL xw = q[1] * q[0] * 2;
  REAL yw = q[2] * q[0] * 2;
  REAL zw = q[3] * q[0] * 2;
  vo[0] = (1 - y2 - z2)*vi[0] + (xy + zw    )*vi[1] + (zx - yw    )*vi[2];
  vo[1] = (xy - zw    )*vi[0] + (1 - z2 - x2)*vi[1] + (yz + xw    )*vi[2];
  vo[2] = (zx + yw    )*vi[0] + (yz - xw    )*vi[1] + (1 - x2 - y2)*vi[2];
//  vo[0] = (1.0 - y2 - z2)*vi[0] + (xy - zw      )*vi[1] + (zx + yw      )*vi[2];
//  vo[1] = (xy + zw      )*vi[0] + (1.0 - z2 - x2)*vi[1] + (yz - xw      )*vi[2];
//  vo[2] = (zx - yw      )*vi[0] + (yz + xw      )*vi[1] + (1.0 - x2 - y2)*vi[2];
}
#ifdef DFM2_STATIC_LIBRARY
template void MyQuatConjVec(float vo[3], const float q[4], const float vi[3]);
template void MyQuatConjVec(double vo[3], const double q[4], const double vi[3]);
#endif

// --------------------------

template <typename REAL>
DFM2_INLINE void MyInverse_Mat3
    (REAL Ainv[9],
     const REAL A[9])
{
  const REAL det =
      + A[0]*A[4]*A[8] + A[3]*A[7]*A[2] + A[6]*A[1]*A[5]
      - A[0]*A[7]*A[5] - A[6]*A[4]*A[2] - A[3]*A[1]*A[8];
  const REAL inv_det = 1.0/det;
  Ainv[0] = inv_det*(A[4]*A[8]-A[5]*A[7]);
  Ainv[1] = inv_det*(A[2]*A[7]-A[1]*A[8]);
  Ainv[2] = inv_det*(A[1]*A[5]-A[2]*A[4]);
  Ainv[3] = inv_det*(A[5]*A[6]-A[3]*A[8]);
  Ainv[4] = inv_det*(A[0]*A[8]-A[2]*A[6]);
  Ainv[5] = inv_det*(A[2]*A[3]-A[0]*A[5]);
  Ainv[6] = inv_det*(A[3]*A[7]-A[4]*A[6]);
  Ainv[7] = inv_det*(A[1]*A[6]-A[0]*A[7]);
  Ainv[8] = inv_det*(A[0]*A[4]-A[1]*A[3]);
}

template <typename T>
DFM2_INLINE void MyMatVec3
    (T y[3],
     const T m[9], const T x[3])
{
  y[0] = m[0]*x[0] + m[1]*x[1] + m[2]*x[2];
  y[1] = m[3]*x[0] + m[4]*x[1] + m[5]*x[2];
  y[2] = m[6]*x[0] + m[7]*x[1] + m[8]*x[2];
}

}
}

// ---------------------------------

template <typename T>
T delfem2::Volume_Tet3(
    const T v1[3],
    const T v2[3],
    const T v3[3],
    const T v4[3])
{
  return
  ((v2[0]-v1[0])*((v3[1]-v1[1])*(v4[2]-v1[2])-(v4[1]-v1[1])*(v3[2]-v1[2]))
   -(v2[1]-v1[1])*((v3[0]-v1[0])*(v4[2]-v1[2])-(v4[0]-v1[0])*(v3[2]-v1[2]))
   +(v2[2]-v1[2])*((v3[0]-v1[0])*(v4[1]-v1[1])-(v4[0]-v1[0])*(v3[1]-v1[1]))
   ) * static_cast<T>(1.0/6.0);
}
#ifdef DFM2_STATIC_LIBRARY
template float delfem2::Volume_Tet3(const float v1[3],
                                    const float v2[3],
                                    const float v3[3],
                                    const float v4[3]);
template double delfem2::Volume_Tet3(const double v1[3],
                                     const double v2[3],
                                     const double v3[3],
                                     const double v4[3]);
#endif
  
// --------------------------------

//! Hight of a tetrahedra
template <typename T>
double delfem2::Height(
    const CVec3<T>& v1,
    const CVec3<T>& v2,
    const CVec3<T>& v3,
    const CVec3<T>& v4)
{
  T n[3];
  NormalTri3(n,
      v1.p,v2.p,v3.p);
  Normalize3(n);
  return (v4.p[0]-v1.p[0])*n[0]+(v4.p[1]-v1.p[1])*n[1]+(v4.p[2]-v1.p[2])*n[2];
}
#ifdef DFM2_STATIC_LIBRARY
template double delfem2::Height(const CVec3f& v1,
                                const CVec3f& v2,
                                const CVec3f& v3,
                                const CVec3f& v4);
template double delfem2::Height(const CVec3d& v1,
                                const CVec3d& v2,
                                const CVec3d& v3,
                                const CVec3d& v4);
#endif

// ------------------------------------------------------------------------------

template <typename T>
delfem2::CVec3<T> delfem2::positionBarycentricCoord_Pyramid
(double r0,
 double r1,
 double r2,
 const CVec3<T>& p0,
 const CVec3<T>& p1,
 const CVec3<T>& p2,
 const CVec3<T>& p3,
 const CVec3<T>& p4)
{
  return (1.0-r2)*(1.0-r0)*(1.0-r1)*p0
  + (1.0-r2)*r0*(1.0-r1)*p1
  + (1.0-r2)*r0*r1*p2
  + (1.0-r2)*(1.0-r0)*r1*p3
  + r2*p4;
}

template <typename T>
delfem2::CVec3<T> delfem2::positionBarycentricCoord_Wedge
(double r0,
 double r1,
 double r2,
 const CVec3<T>& p0,
 const CVec3<T>& p1,
 const CVec3<T>& p2,
 const CVec3<T>& p3,
 const CVec3<T>& p4,
 const CVec3<T>& p5)
{
  return (1.0-r2)*r0*p0
  + (1.0-r2)*r1*p1
  + (1.0-r2)*(1.0-r0-r1)*p2
  + r2*r0*p3
  + r2*r1*p4
  + r2*(1.0-r0-r1)*p5;
}

template <typename T>
void delfem2::iteration_barycentricCoord_Origin_Solid
(double& r0,
 double& r1,
 double& r2,
 const CVec3<T>& q, // q=positionBarycentricCoord_Wedge
 const CVec3<T>& dpdr0,
 const CVec3<T>& dpdr1,
 const CVec3<T>& dpdr2,
 double damp)
{
  namespace lcl = delfem2::solidelm;
  const double cxx = dpdr0.p[0]*dpdr0.p[0] + dpdr1.p[0]*dpdr1.p[0] + dpdr2.p[0]*dpdr2.p[0];
  const double cxy = dpdr0.p[0]*dpdr0.p[1] + dpdr1.p[0]*dpdr1.p[1] + dpdr2.p[0]*dpdr2.p[1];
  const double cxz = dpdr0.p[0]*dpdr0.p[2] + dpdr1.p[0]*dpdr1.p[2] + dpdr2.p[0]*dpdr2.p[2];
  const double cyy = dpdr0.p[1]*dpdr0.p[1] + dpdr1.p[1]*dpdr1.p[1] + dpdr2.p[1]*dpdr2.p[1];
  const double cyz = dpdr0.p[1]*dpdr0.p[2] + dpdr1.p[1]*dpdr1.p[2] + dpdr2.p[1]*dpdr2.p[2];
  const double czz = dpdr0.p[2]*dpdr0.p[2] + dpdr1.p[2]*dpdr1.p[2] + dpdr2.p[2]*dpdr2.p[2];
  double C[9] = {cxx,cxy,cxz, cxy,cyy,cyz, cxz,cyz,czz};
  double Cinv[9]; lcl::MyInverse_Mat3(Cinv, C);
  const CVec3<T> d = damp*Mat3Vec(Cinv,q);
  r0 -= dpdr0*d;
  r1 -= dpdr1*d;
  r2 -= dpdr2*d;
}

template <typename T>
bool delfem2::barycentricCoord_Origin_Tet(
    double& r0,
    double& r1,
    double& r2,
    const CVec3<T>& p0,
    const CVec3<T>& p1,
    const CVec3<T>& p2,
    const CVec3<T>& p3)
{
  double v0 = Volume_OrgTet(p1, p2, p3);
  double v1 = Volume_OrgTet(p2, p0, p3);
  double v2 = Volume_OrgTet(p1, p3, p0);
  double v3 = Volume_OrgTet(p1, p0, p2);
  double vt_inv = 1.0/(v0+v1+v2+v3);
  r0 = v0*vt_inv;
  r1 = v1*vt_inv;
  r2 = v2*vt_inv;
  return true;
}

template <typename T>
bool delfem2::barycentricCoord_Origin_Pyramid(
    double& r0,
    double& r1,
    double& r2,
    const CVec3<T>& p0,
    const CVec3<T>& p1,
    const CVec3<T>& p2,
    const CVec3<T>& p3,
    const CVec3<T>& p4)
{
  CVec3<T> q = positionBarycentricCoord_Pyramid(r0,r1,r2, p0,p1,p2,p3,p4);
  for(int itr=0;itr<5;++itr){
    const CVec3<T> dpdr0 = -(1.0-r2)*(1.0-r1)*p0 + (1.0-r2)*(1.0-r1)*p1 + (1.0-r2)*r1*p2 - (1.0-r2)*r1*p3;
    const CVec3<T> dpdr1 = -(1.0-r2)*(1.0-r0)*p0 - (1.0-r2)*r0*p1 + (1.0-r2)*r0*p2 + (1.0-r2)*(1.0-r0)*p3;
    const CVec3<T> dpdr2 = -(1.0-r0)*(1.0-r1)*p0 - r0*(1.0-r1)*p1 - r0*r1*p2 - (1.0-r0)*r1*p3 + p4;
    iteration_barycentricCoord_Origin_Solid(r0,r1,r2, q,
                                            dpdr0,dpdr1,dpdr2,
                                            1.0);
    q = positionBarycentricCoord_Pyramid(r0,r1,r2, p0,p1,p2,p3,p4);
  }
  return true;
}

template <typename T>
bool delfem2::barycentricCoord_Origin_Wedge(
    double& r0,
    double& r1,
    double& r2,
    const CVec3<T>& p0,
    const CVec3<T>& p1,
    const CVec3<T>& p2,
    const CVec3<T>& p3,
    const CVec3<T>& p4,
    const CVec3<T>& p5)
{
  CVec3<T> q = positionBarycentricCoord_Wedge(r0,r1,r2, p0,p1,p2,p3,p4,p5);
  for(int itr=0;itr<5;++itr){
    const CVec3<T> dpdr0 = (1.0-r2)*(p0-p2)+r2*(p3-p5);
    const CVec3<T> dpdr1 = (1.0-r2)*(p1-p2)+r2*(p4-p5);
    const CVec3<T> dpdr2 = r0*(p3-p0)+r1*(p4-p1)+(1.0-r0-r1)*(p5-p2);
    iteration_barycentricCoord_Origin_Solid(r0,r1,r2, q,
                                            dpdr0,dpdr1,dpdr2,
                                            1.0);
    q = positionBarycentricCoord_Wedge(r0,r1,r2, p0,p1,p2,p3,p4,p5);
  }
  return true;
}


// ----------------------------------------------------------------------

//! Volume of a tetrahedra
template <typename T>
T delfem2::Volume_Tet(
    const CVec3<T>& v0,
    const CVec3<T>& v1,
    const CVec3<T>& v2,
    const CVec3<T>& v3 )
{
  return delfem2::Volume_Tet3(v0.p, v1.p, v2.p, v3.p);
  /*
  double v = (v1.p[0]-v0.p[0])*( (v2.p[1]-v0.p[1])*(v3.p[2]-v0.p[2]) - (v3.p[1]-v0.p[1])*(v2.p[2]-v0.p[2]) )
           + (v1.p[1]-v0.p[1])*( (v2.p[2]-v0.p[2])*(v3.p[0]-v0.p[0]) - (v3.p[2]-v0.p[2])*(v2.p[0]-v0.p[0]) )
           + (v1.p[2]-v0.p[2])*( (v2.p[0]-v0.p[0])*(v3.p[1]-v0.p[1]) - (v3.p[0]-v0.p[0])*(v2.p[1]-v0.p[1]) );
  return v*0.16666666666666666666666666666667;
   */
}
#ifdef DFM2_STATIC_LIBRARY
template double delfem2::Volume_Tet(const CVec3d& v0, const CVec3d& v1, const CVec3d& v2, const CVec3d& v3 );
template float delfem2::Volume_Tet(const CVec3f& v0, const CVec3f& v1, const CVec3f& v2, const CVec3f& v3 );
#endif


//! Volume of a tetrahedra v0 is orgin
template <typename T>
T delfem2::Volume_OrgTet(
    const CVec3<T>& v1,
    const CVec3<T>& v2,
    const CVec3<T>& v3 )
{
  T v = v1.p[0]*(v2.p[1]*v3.p[2]-v3.p[1]*v2.p[2])
      + v1.p[1]*(v2.p[2]*v3.p[0]-v3.p[2]*v2.p[0])
      + v1.p[2]*(v2.p[0]*v3.p[1]-v3.p[0]*v2.p[1]);
  return v * static_cast<T>(1.0 / 6.0); 
}
#ifdef DFM2_STATIC_LIBRARY
template double delfem2::Volume_OrgTet(
    const CVec3d& v1,
    const CVec3d& v2,
    const CVec3d& v3 );
template float delfem2::Volume_OrgTet(
    const CVec3f& v1,
    const CVec3f& v2,
    const CVec3f& v3 );
#endif

template <typename T>
double delfem2::Volume_Pyramid
(const CVec3<T>& p0,
 const CVec3<T>& p1,
 const CVec3<T>& p2,
 const CVec3<T>& p3,
 const CVec3<T>& p4)
{
  double v0124 = Volume_Tet(p0, p1, p2, p4);
  double v0234 = Volume_Tet(p0, p2, p3, p4);
  double v0134 = Volume_Tet(p0, p1, p3, p4);
  double v2314 = Volume_Tet(p2, p3, p1, p4);
  double v0 = v0124+v0234;
  double v1 = v0134+v2314;
  return (v0+v1)*0.5;
}

template <typename T>
double delfem2::Volume_Wedge(
    const CVec3<T>& p0,
    const CVec3<T>& p1,
    const CVec3<T>& p2,
    const CVec3<T>& p3,
    const CVec3<T>& p4,
    const CVec3<T>& p5)
{
  CVec3<T> pm = (p0+p1+p2+p3+p4+p5)/6.0;
  double vm012 = Volume_Tet(pm,p0,p1,p2);
  double vm435 = Volume_Tet(pm,p4,p3,p5);
  double vp0143 = Volume_Pyramid(p0,p1,p4,p3,pm);
  double vp1254 = Volume_Pyramid(p1,p2,p5,p4,pm);
  double vp2035 = Volume_Pyramid(p2,p2,p3,p5,pm);
  return vm012+vm435+vp0143+vp1254+vp2035;
}

// ---------------------------------------------------------------

template <typename T>
double delfem2::SolidAngleTri(
    const CVec3<T>& v1,
    const CVec3<T>& v2,
    const CVec3<T>& v3)
{
  double l1 = v1.norm();
  double l2 = v2.norm();
  double l3 = v3.norm();
  double den = (v1^v2).dot(v3);
  double num = l1*l2*l3+(v1.dot(v2))*l3+(v2.dot(v3))*l1+(v3.dot(v1))*l2;
  double tho = den/num;
  double v = atan(tho);
  if (v<0){ v += 2*M_PI; }
  v *= 2;
  return v;
}
#ifdef DFM2_STATIC_LIBRARY
template double delfem2::SolidAngleTri(const CVec3d& v1, const CVec3d& v2, const CVec3d& v3);
#endif



// ------------------------------------------------------------------

template <typename T>
double delfem2::SqareLongestEdgeLength(
    const CVec3<T>& ipo0,
    const CVec3<T>& ipo1,
    const CVec3<T>& ipo2,
    const CVec3<T>& ipo3 )
{
  double edge1, edge2;
  edge1 = SquareDistance( ipo0, ipo1 );
  edge2 = SquareDistance( ipo0, ipo2 );
  if( edge2 > edge1 ) edge1 = edge2;
  edge2 = SquareDistance( ipo0, ipo3 );
  if( edge2 > edge1 ) edge1 = edge2;
  edge2 = SquareDistance( ipo1, ipo2 );
  if( edge2 > edge1 ) edge1 = edge2;
  edge2 = SquareDistance( ipo1, ipo3 );
  if( edge2 > edge1 ) edge1 = edge2;
  edge2 = SquareDistance( ipo2, ipo3 );
  if( edge2 > edge1 ) edge1 = edge2;
  return edge1;
}

// --------------------------------------

template <typename T>
double delfem2::LongestEdgeLength(
    const CVec3<T>& ipo0,
    const CVec3<T>& ipo1,
    const CVec3<T>& ipo2,
    const CVec3<T>& ipo3 )
{
  return sqrt( SqareLongestEdgeLength(ipo0,ipo1,ipo2,ipo3) );
}

// --------------------------------------

template <typename T>
double delfem2::SqareShortestEdgeLength(
    const CVec3<T>& ipo0,
    const CVec3<T>& ipo1,
    const CVec3<T>& ipo2,
    const CVec3<T>& ipo3 )
{
  double edge1, edge2;
  edge1 = SquareDistance( ipo0, ipo1 );
  edge2 = SquareDistance( ipo0, ipo2 );
  if( edge2 < edge1 ) edge1 = edge2;
  edge2 = SquareDistance( ipo0, ipo3 );
  if( edge2 < edge1 ) edge1 = edge2;
  edge2 = SquareDistance( ipo1, ipo2 );
  if( edge2 < edge1 ) edge1 = edge2;
  edge2 = SquareDistance( ipo1, ipo3 );
  if( edge2 < edge1 ) edge1 = edge2;
  edge2 = SquareDistance( ipo2, ipo3 );
  if( edge2 < edge1 ) edge1 = edge2;
  return edge1;
}

// -----------------------------------------

template <typename T>
double delfem2::ShortestEdgeLength
(const CVec3<T>& ipo0,
 const CVec3<T>& ipo1,
 const CVec3<T>& ipo2,
 const CVec3<T>& ipo3 )
{
  return sqrt( SqareShortestEdgeLength(ipo0,ipo1,ipo2,ipo3) );
}


template <typename T>
double delfem2::Volume_Tet(
    int iv1,
    int iv2,
    int iv3,
    int iv4,
    const std::vector<CVec3<T>>& aPoint)
{
  return Volume_Tet(aPoint[iv1],aPoint[iv2],aPoint[iv3],aPoint[iv4]);
}
#ifdef DFM2_STATIC_LIBRARY
template double delfem2::Volume_Tet(
    int iv1,
    int iv2,
    int iv3,
    int iv4,
   const std::vector<CVec3d>& aPoint);
#endif


