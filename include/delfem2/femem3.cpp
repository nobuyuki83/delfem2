/*
 * Copyright (c) 2019 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */


#include <math.h>
#include <assert.h>
#include <complex>

#include "delfem2/femem3.h"

// --------------------------------------------------------

namespace delfem2 {
namespace femem3 {


const static unsigned int NIntLineGauss[4] = {
  1, 2, 3, 4
};
const static double LineGauss[4][4][2] =
{
  {
    { 0.0, 2.0 },
    { 0.0, 0.0 },
    { 0.0, 0.0 },
    { 0.0, 0.0 },
  },
  {
    { -0.577350269189626, 1.0 },
    {  0.577350269189626, 1.0 },
    {  0.0,               0.0 },
    {  0.0,               0.0 },
  },
  {
    { -0.774596669241483, 0.555555555555556 },
    {  0.0,               0.888888888888889 },
    {  0.774596669241483, 0.555555555555556 },
    {  0.0,               0.0               },
  },
  {
    { -0.861136311594053, 0.347854845137454 },
    { -0.339981043584856, 0.652145154862546 },
    {  0.339981043584856, 0.652145154862546 },
    {  0.861136311594053, 0.347854845137454 },
  }
};


DFM2_INLINE void ShapeFunc_Hex8
(const double& r0, const double& r1,	const double& r2,
 const double coords[][3],
 double& detjac,
 double dndx[][3],
 double an[] )
{
  an[0] = 0.125*(1.0-r0)*(1.0-r1)*(1.0-r2);
  an[1] = 0.125*(1.0+r0)*(1.0-r1)*(1.0-r2);
  an[2] = 0.125*(1.0-r0)*(1.0+r1)*(1.0-r2);
  an[3] = 0.125*(1.0+r0)*(1.0+r1)*(1.0-r2);
  an[4] = 0.125*(1.0-r0)*(1.0-r1)*(1.0+r2);
  an[5] = 0.125*(1.0+r0)*(1.0-r1)*(1.0+r2);
  an[6] = 0.125*(1.0-r0)*(1.0+r1)*(1.0+r2);
  an[7] = 0.125*(1.0+r0)*(1.0+r1)*(1.0+r2);
  
  double dndr[8][3];
  dndr[0][0] = -0.125*(1.0-r1)*(1.0-r2);
  dndr[1][0] = -dndr[0][0];
  dndr[2][0] = -0.125*(1.0+r1)*(1.0-r2);
  dndr[3][0] = -dndr[2][0];
  dndr[4][0] = -0.125*(1.0-r1)*(1.0+r2);
  dndr[5][0] = -dndr[4][0];
  dndr[6][0] = -0.125*(1.0+r1)*(1.0+r2);
  dndr[7][0] = -dndr[6][0];
  
  dndr[0][1] = -0.125*(1.0-r0)*(1.0-r2);
  dndr[1][1] = -0.125*(1.0+r0)*(1.0-r2);
  dndr[2][1] = -dndr[0][1];
  dndr[3][1] = -dndr[1][1];
  dndr[4][1] = -0.125*(1.0-r0)*(1.0+r2);
  dndr[5][1] = -0.125*(1.0+r0)*(1.0+r2);
  dndr[6][1] = -dndr[4][1];
  dndr[7][1] = -dndr[5][1];
  
  dndr[0][2] = -0.125*(1.0-r0)*(1.0-r1);
  dndr[1][2] = -0.125*(1.0+r0)*(1.0-r1);
  dndr[2][2] = -0.125*(1.0-r0)*(1.0+r1);
  dndr[3][2] = -0.125*(1.0+r0)*(1.0+r1);
  dndr[4][2] = -dndr[0][2];
  dndr[5][2] = -dndr[1][2];
  dndr[6][2] = -dndr[2][2];
  dndr[7][2] = -dndr[3][2];
  
  double dxdr[3][3]  = {
    { 0.0, 0.0, 0.0 },
    { 0.0, 0.0, 0.0 },
    { 0.0, 0.0, 0.0 },
  };
  
  for(int inode=0;inode<8;inode++){
    dxdr[0][0] += coords[inode][0]*dndr[inode][0];
    dxdr[0][1] += coords[inode][0]*dndr[inode][1];
    dxdr[0][2] += coords[inode][0]*dndr[inode][2];
    dxdr[1][0] += coords[inode][1]*dndr[inode][0];
    dxdr[1][1] += coords[inode][1]*dndr[inode][1];
    dxdr[1][2] += coords[inode][1]*dndr[inode][2];
    dxdr[2][0] += coords[inode][2]*dndr[inode][0];
    dxdr[2][1] += coords[inode][2]*dndr[inode][1];
    dxdr[2][2] += coords[inode][2]*dndr[inode][2];
  }
  
  detjac = dxdr[0][0]*dxdr[1][1]*dxdr[2][2]
  + dxdr[1][0]*dxdr[2][1]*dxdr[0][2]
  + dxdr[2][0]*dxdr[0][1]*dxdr[1][2]
  - dxdr[0][0]*dxdr[2][1]*dxdr[1][2]
  - dxdr[1][0]*dxdr[0][1]*dxdr[2][2]
  - dxdr[2][0]*dxdr[1][1]*dxdr[0][2];
  
  const double inv_jac = 1.0 / detjac;
  
  double drdx[3][3];
  drdx[0][0] = inv_jac*( dxdr[1][1]*dxdr[2][2]-dxdr[1][2]*dxdr[2][1] );
  drdx[0][1] = inv_jac*( dxdr[0][2]*dxdr[2][1]-dxdr[0][1]*dxdr[2][2] );
  drdx[0][2] = inv_jac*( dxdr[0][1]*dxdr[1][2]-dxdr[0][2]*dxdr[1][1] );
  drdx[1][0] = inv_jac*( dxdr[1][2]*dxdr[2][0]-dxdr[1][0]*dxdr[2][2] );
  drdx[1][1] = inv_jac*( dxdr[0][0]*dxdr[2][2]-dxdr[0][2]*dxdr[2][0] );
  drdx[1][2] = inv_jac*( dxdr[0][2]*dxdr[1][0]-dxdr[0][0]*dxdr[1][2] );
  drdx[2][0] = inv_jac*( dxdr[1][0]*dxdr[2][1]-dxdr[1][1]*dxdr[2][0] );
  drdx[2][1] = inv_jac*( dxdr[0][1]*dxdr[2][0]-dxdr[0][0]*dxdr[2][1] );
  drdx[2][2] = inv_jac*( dxdr[0][0]*dxdr[1][1]-dxdr[0][1]*dxdr[1][0] );
  
  for(int inode=0;inode<8;inode++){
    dndx[inode][0] = dndr[inode][0]*drdx[0][0] + dndr[inode][1]*drdx[1][0] + dndr[inode][2]*drdx[2][0];
    dndx[inode][1] = dndr[inode][0]*drdx[0][1] + dndr[inode][1]*drdx[1][1] + dndr[inode][2]*drdx[2][1];
    dndx[inode][2] = dndr[inode][0]*drdx[0][2] + dndr[inode][1]*drdx[1][2] + dndr[inode][2]*drdx[2][2];
  }
}

// area of a triangle
DFM2_INLINE double TriArea2D(const double p0[], const double p1[], const double p2[]){
  return 0.5*((p1[0]-p0[0])*(p2[1]-p0[1])-(p2[0]-p0[0])*(p1[1]-p0[1]));
}

DFM2_INLINE double TriArea3D(const double v1[3], const double v2[3], const double v3[3]){
  double x, y, z;
  x = ( v2[1] - v1[1] )*( v3[2] - v1[2] ) - ( v3[1] - v1[1] )*( v2[2] - v1[2] );
  y = ( v2[2] - v1[2] )*( v3[0] - v1[0] ) - ( v3[2] - v1[2] )*( v2[0] - v1[0] );
  z = ( v2[0] - v1[0] )*( v3[1] - v1[1] ) - ( v3[0] - v1[0] )*( v2[1] - v1[1] );
  return 0.5*sqrt( x*x + y*y + z*z );
}

DFM2_INLINE double Distance3D(const double p0[3], const double p1[3]){
  return sqrt( (p1[0]-p0[0])*(p1[0]-p0[0]) + (p1[1]-p0[1])*(p1[1]-p0[1]) + (p1[2]-p0[2])*(p1[2]-p0[2]) );
}

DFM2_INLINE double Dot3D(const double a[], const double b[]){
  return a[0]*b[0]+a[1]*b[1]+a[2]*b[2];
}

DFM2_INLINE void Cross3D(double r[3], const double v1[3], const double v2[3]){
  r[0] = v1[1]*v2[2] - v2[1]*v1[2];
  r[1] = v1[2]*v2[0] - v2[2]*v1[0];
  r[2] = v1[0]*v2[1] - v2[0]*v1[1];
}

DFM2_INLINE void  UnitNormalAreaTri3D(double n[3], double& a, const double v1[3], const double v2[3], const double v3[3]){
  n[0] = ( v2[1] - v1[1] )*( v3[2] - v1[2] ) - ( v3[1] - v1[1] )*( v2[2] - v1[2] );
  n[1] = ( v2[2] - v1[2] )*( v3[0] - v1[0] ) - ( v3[2] - v1[2] )*( v2[0] - v1[0] );
  n[2] = ( v2[0] - v1[0] )*( v3[1] - v1[1] ) - ( v3[0] - v1[0] )*( v2[1] - v1[1] );
  a = sqrt(n[0]*n[0]+n[1]*n[1]+n[2]*n[2])*0.5;
  const double invlen = 0.5/a;
  n[0]*=invlen;	n[1]*=invlen;	n[2]*=invlen;
}

DFM2_INLINE void TriDlDx(double dldx[][2],
                         double const_term[],
                         const double p0[],
                         const double p1[],
                         const double p2[])
{
  const double area = femem3::TriArea2D(p0, p1, p2);
  const double tmp1 = 0.5/area;
  
  const_term[0] = tmp1*(p1[0]*p2[1]-p2[0]*p1[1]);
  const_term[1] = tmp1*(p2[0]*p0[1]-p0[0]*p2[1]);
  const_term[2] = tmp1*(p0[0]*p1[1]-p1[0]*p0[1]);
  
  dldx[0][0] = tmp1*(p1[1]-p2[1]);
  dldx[1][0] = tmp1*(p2[1]-p0[1]);
  dldx[2][0] = tmp1*(p0[1]-p1[1]);
  
  dldx[0][1] = tmp1*(p2[0]-p1[0]);
  dldx[1][1] = tmp1*(p0[0]-p2[0]);
  dldx[2][1] = tmp1*(p1[0]-p0[0]);
  /*
   assert( fabs( dldx[0][0]+dldx[1][0]+dldx[2][0] ) < 1.0e-15 );
   assert( fabs( dldx[0][1]+dldx[1][1]+dldx[2][1] ) < 1.0e-15 );
   
   assert( fabs( const_term[0]+dldx[0][0]*p0[0]+dldx[0][1]*p0[1] - 1.0 ) < 1.0e-10 );
   assert( fabs( const_term[0]+dldx[0][0]*p1[0]+dldx[0][1]*p1[1] ) < 1.0e-10 );
   assert( fabs( const_term[0]+dldx[0][0]*p2[0]+dldx[0][1]*p2[1] ) < 1.0e-10 );
   
   assert( fabs( const_term[1]+dldx[1][0]*p0[0]+dldx[1][1]*p0[1] ) < 1.0e-10 );
   assert( fabs( const_term[1]+dldx[1][0]*p1[0]+dldx[1][1]*p1[1] - 1.0 ) < 1.0e-10 );
   assert( fabs( const_term[1]+dldx[1][0]*p2[0]+dldx[1][1]*p2[1] ) < 1.0e-10 );
   
   assert( fabs( const_term[2]+dldx[2][0]*p0[0]+dldx[2][1]*p0[1] ) < 1.0e-10 );
   assert( fabs( const_term[2]+dldx[2][0]*p1[0]+dldx[2][1]*p1[1] ) < 1.0e-10 );
   assert( fabs( const_term[2]+dldx[2][0]*p2[0]+dldx[2][1]*p2[1] - 1.0 ) < 1.0e-10 );
   */
}


DFM2_INLINE void MakeConstMatrix3D
 (double C[6][6],
  double lambda,
  double myu,
  const double Gu[3][3])
{
  const double GuGu2[6] = {
    Dot3D(Gu[0],Gu[0]), // 0 xx
    Dot3D(Gu[1],Gu[1]), // 1 yy
    Dot3D(Gu[2],Gu[2]), // 2 zz
    Dot3D(Gu[0],Gu[1]), // 3 xy
    Dot3D(Gu[1],Gu[2]), // 4 yz
    Dot3D(Gu[2],Gu[0])  // 5 zx
  };
  C[0][0] = lambda*GuGu2[0]*GuGu2[0] + 2*myu*(GuGu2[0]*GuGu2[0]); // 00(0):00(0) 00(0):00(0)
  C[0][1] = lambda*GuGu2[0]*GuGu2[1] + 2*myu*(GuGu2[3]*GuGu2[3]); // 00(0):11(1) 01(3):01(3)
  C[0][2] = lambda*GuGu2[0]*GuGu2[2] + 2*myu*(GuGu2[5]*GuGu2[5]); // 00(0):22(2) 02(5):02(5)
  C[0][3] = lambda*GuGu2[0]*GuGu2[3] + 2*myu*(GuGu2[0]*GuGu2[3]); // 00(0):01(3) 00(0):01(3)
  C[0][4] = lambda*GuGu2[0]*GuGu2[4] + 2*myu*(GuGu2[3]*GuGu2[5]); // 00(0):12(4) 01(3):02(5)
  C[0][5] = lambda*GuGu2[0]*GuGu2[5] + 2*myu*(GuGu2[0]*GuGu2[5]); // 00(0):20(5) 00(0):02(5)
  C[1][0] = lambda*GuGu2[1]*GuGu2[0] + 2*myu*(GuGu2[3]*GuGu2[3]); // 11(1):00(0) 01(3):01(3)
  C[1][1] = lambda*GuGu2[1]*GuGu2[1] + 2*myu*(GuGu2[1]*GuGu2[1]); // 11(1):11(1) 11(1):11(1)
  C[1][2] = lambda*GuGu2[1]*GuGu2[2] + 2*myu*(GuGu2[4]*GuGu2[4]); // 11(1):22(2) 12(4):12(4)
  C[1][3] = lambda*GuGu2[1]*GuGu2[3] + 2*myu*(GuGu2[1]*GuGu2[3]); // 11(1):01(3) 11(1):01(3)
  C[1][4] = lambda*GuGu2[1]*GuGu2[4] + 2*myu*(GuGu2[1]*GuGu2[4]); // 11(1):12(4) 11(1):12(4)
  C[1][5] = lambda*GuGu2[1]*GuGu2[5] + 2*myu*(GuGu2[3]*GuGu2[4]); // 11(1):20(5) 12(4):10(3)
  C[2][0] = lambda*GuGu2[2]*GuGu2[0] + 2*myu*(GuGu2[5]*GuGu2[5]); // 22(2):00(0) 02(5):02(5)
  C[2][1] = lambda*GuGu2[2]*GuGu2[1] + 2*myu*(GuGu2[4]*GuGu2[4]); // 22(2):11(1) 12(4):12(4)
  C[2][2] = lambda*GuGu2[2]*GuGu2[2] + 2*myu*(GuGu2[2]*GuGu2[2]); // 22(2):22(2) 22(2):22(2)
  C[2][3] = lambda*GuGu2[2]*GuGu2[3] + 2*myu*(GuGu2[4]*GuGu2[5]); // 22(2):01(3) 12(4):02(5)
  C[2][4] = lambda*GuGu2[2]*GuGu2[4] + 2*myu*(GuGu2[2]*GuGu2[4]); // 22(2):12(4) 22(2):12(4)
  C[2][5] = lambda*GuGu2[2]*GuGu2[5] + 2*myu*(GuGu2[2]*GuGu2[5]); // 22(2):02(5) 22(2):02(5)
  C[3][0] = lambda*GuGu2[3]*GuGu2[0] + 2*myu*(GuGu2[3]*GuGu2[0]); // 01(3):00(0) 00(0):01(3)
  C[3][1] = lambda*GuGu2[3]*GuGu2[1] + 2*myu*(GuGu2[3]*GuGu2[1]); // 01(3):11(1) 11(1):01(3)
  C[3][2] = lambda*GuGu2[3]*GuGu2[2] + 2*myu*(GuGu2[4]*GuGu2[5]); // 01(3):22(2) 12(4):02(5)
  C[3][3] = lambda*GuGu2[3]*GuGu2[3] + 1*myu*(GuGu2[0]*GuGu2[1] + GuGu2[3]*GuGu2[3]); // 01(3):01(3) 00(0):11(1) 01(3):01(3)
  C[3][4] = lambda*GuGu2[3]*GuGu2[4] + 1*myu*(GuGu2[3]*GuGu2[4] + GuGu2[1]*GuGu2[5]); // 01(3):12(4) 01(3):12(4) 11(1):02(5)
  C[3][5] = lambda*GuGu2[3]*GuGu2[5] + 1*myu*(GuGu2[5]*GuGu2[3] + GuGu2[0]*GuGu2[4]); // 01(3):20(5) 02(5):10(3) 00(0):12(4)
  C[4][0] = lambda*GuGu2[4]*GuGu2[0] + 2*myu*(GuGu2[3]*GuGu2[5]); // 12(4):00(0) 01(3):02(5)
  C[4][1] = lambda*GuGu2[4]*GuGu2[1] + 2*myu*(GuGu2[1]*GuGu2[4]); // 12(4):11(1) 11(1):12(4)
  C[4][2] = lambda*GuGu2[4]*GuGu2[2] + 2*myu*(GuGu2[2]*GuGu2[4]); // 12(4):22(2) 22(2):12(4)
  C[4][3] = lambda*GuGu2[4]*GuGu2[3] + 1*myu*(GuGu2[3]*GuGu2[4] + GuGu2[1]*GuGu2[5]); // 12(4):01(3) 10(3):21(4) 11(1):20(5)
  C[4][4] = lambda*GuGu2[4]*GuGu2[4] + 1*myu*(GuGu2[1]*GuGu2[2] + GuGu2[4]*GuGu2[4]); // 12(4):12(4) 11(1):22(2) 12(4):21(4)
  C[4][5] = lambda*GuGu2[4]*GuGu2[5] + 1*myu*(GuGu2[4]*GuGu2[5] + GuGu2[3]*GuGu2[2]); // 12(4):20(5) 12(4):20(5) 10(3):22(2)
  C[5][0] = lambda*GuGu2[5]*GuGu2[0] + 2*myu*(GuGu2[0]*GuGu2[5]); // 02(5):00(0) 00(0):02(5)
  C[5][1] = lambda*GuGu2[5]*GuGu2[1] + 2*myu*(GuGu2[3]*GuGu2[4]); // 02(5):11(1) 10(3):12(4)
  C[5][2] = lambda*GuGu2[5]*GuGu2[2] + 2*myu*(GuGu2[2]*GuGu2[5]); // 02(5):22(2) 22(2):02(5)
  C[5][3] = lambda*GuGu2[5]*GuGu2[3] + 1*myu*(GuGu2[0]*GuGu2[4] + GuGu2[3]*GuGu2[5]); // 02(5):01(3) 00(0):21(4) 01(3):20(5)
  C[5][4] = lambda*GuGu2[5]*GuGu2[4] + 1*myu*(GuGu2[3]*GuGu2[2] + GuGu2[5]*GuGu2[4]); // 02(5):12(4) 01(3):22(2) 02(5):21(4)
  C[5][5] = lambda*GuGu2[5]*GuGu2[5] + 1*myu*(GuGu2[5]*GuGu2[5] + GuGu2[0]*GuGu2[2]); // 02(5):20(5) 02(5):20(5) 00(0):22(2)
}

DFM2_INLINE void MakePositiveDefinite_Sim22
 (const double s2[3],double s3[3])
{
  const double b = (s2[0]+s2[1])*0.5;
  const double d = (s2[0]-s2[1])*(s2[0]-s2[1])*0.25 + s2[2]*s2[2];
  const double e = sqrt(d);
  if( b-e > 1.0e-20 ){
    s3[0] = s2[0];
    s3[1] = s2[1];
    s3[2] = s2[2];
    return;
  }
  if( b+e < 0 ){
    s3[0] = 0;
    s3[1] = 0;
    s3[2] = 0;
    return;
  }
  const double l = b+e;
  double t0[2] = { s2[0]-l, s2[2]   };
  double t1[2] = { s2[2],   s2[1]-l };
  //  std::cout << t0[0]*t1[1]-t0[1]*t1[0] << std::endl;
  const double sqlen_t0 = t0[0]*t0[0]+t0[1]*t0[1];
  const double sqlen_t1 = t1[0]*t1[0]+t1[1]*t1[1];
  if( sqlen_t0 > sqlen_t1 ){
    if( sqlen_t0 < 1.0e-20 ){
      s3[0] = 0;
      s3[1] = 0;
      s3[2] = 0;
      return;
    }
    const double invlen_t0 = 1.0/sqrt(sqlen_t0);
    t0[0] *= invlen_t0;
    t0[1] *= invlen_t0;
    s3[0] = l*t0[0]*t0[0];
    s3[1] = l*t0[1]*t0[1];
    s3[2] = l*t0[0]*t0[1];
  }
  else{
    if( sqlen_t1 < 1.0e-20 ){
      s3[0] = 0;
      s3[1] = 0;
      s3[2] = 0;
      return;
    }
    const double invlen_t1 = 1.0/sqrt(sqlen_t1);
    t1[0] *= invlen_t1;
    t1[1] *= invlen_t1;
    s3[0] = l*t1[0]*t1[0];
    s3[1] = l*t1[1]*t1[1];
    s3[2] = l*t1[0]*t1[1];
  }
  return;
}

DFM2_INLINE void MakeCurvetureDKT
 (double B1[][3], double B2[][3][2],
  const double coord0[], const double coord1[], const double coord2[],
  const double l1, const double l2 )
{
  const double l0 = 1-l1-l2;
  const double vec0[2] = { coord2[0]-coord1[0], coord2[1]-coord1[1] };
  const double vec1[2] = { coord0[0]-coord2[0], coord0[1]-coord2[1] };
  const double vec2[2] = { coord1[0]-coord0[0], coord1[1]-coord0[1] };
  const double invsqlen0 = 1.0/(vec0[0]*vec0[0]+vec0[1]*vec0[1]);
  const double invsqlen1 = 1.0/(vec1[0]*vec1[0]+vec1[1]*vec1[1]);
  const double invsqlen2 = 1.0/(vec2[0]*vec2[0]+vec2[1]*vec2[1]);
  double p0=-6*vec0[0]*invsqlen0, q0=3*vec0[0]*vec0[1]*invsqlen0, r0=3*vec0[1]*vec0[1]*invsqlen0, t0=-6*vec0[1]*invsqlen0;
  double p1=-6*vec1[0]*invsqlen1, q1=3*vec1[0]*vec1[1]*invsqlen1, r1=3*vec1[1]*vec1[1]*invsqlen1, t1=-6*vec1[1]*invsqlen1;
  double p2=-6*vec2[0]*invsqlen2, q2=3*vec2[0]*vec2[1]*invsqlen2, r2=3*vec2[1]*vec2[1]*invsqlen2, t2=-6*vec2[1]*invsqlen2;
  
  double H1[4][3];
  H1[0][0]=-(l1-l0)*t2+l2*t1;  H1[0][1]=+(l1-l0)*t2+l2*t0;  H1[0][2]=-l2*(t0+t1);
  H1[1][0]= (l2-l0)*t1-l1*t2;  H1[1][1]=+l1*(t0+t2);    H1[1][2]=-l1*t0-(l2-l0)*t1;
  H1[2][0]= (l1-l0)*p2-l2*p1;  H1[2][1]=-(l1-l0)*p2-l2*p0;  H1[2][2]=+l2*(p0+p1);
  H1[3][0]=-(l2-l0)*p1+l1*p2;  H1[3][1]=-l1*(p0+p2);    H1[3][2]=+l1*p0+(l2-l0)*p1;
  
  double H2[4][3][2];
  H2[0][0][0]=-1+(l1-l0)*r2+l2*r1;    H2[0][0][1]=-(l1-l0)*q2-l2*q1;
  H2[0][1][0]= 1+(l1-l0)*r2-l2*r0;    H2[0][1][1]=-(l1-l0)*q2+l2*q0;
  H2[0][2][0]=-l2*(r0-r1);        H2[0][2][1]= l2*(q0-q1);
  
  H2[1][0][0]=-1+l1*r2+(l2-l0)*r1;    H2[1][0][1]=-l1*q2-(l2-l0)*q1;
  H2[1][1][0]=-l1*(r0-r2);        H2[1][1][1]= l1*(q0-q2);
  H2[1][2][0]= 1-l1*r0+(l2-l0)*r1;    H2[1][2][1]= l1*q0-(l2-l0)*q1;
  
  H2[2][0][0]=-(l1-l0)*q2-l2*q1;      H2[2][0][1]= 2-6*l0-(l1-l0)*r2-l2*r1;
  H2[2][1][0]=-(l1-l0)*q2+l2*q0;      H2[2][1][1]=-2+6*l1-(l1-l0)*r2+l2*r0;
  H2[2][2][0]= l2*(q0-q1);        H2[2][2][1]= l2*(r0-r1);
  
  H2[3][0][0]=-l1*q2-(l2-l0)*q1;      H2[3][0][1]= 2-6*l0-l1*r2-(l2-l0)*r1;
  H2[3][1][0]= l1*(q0-q2);        H2[3][1][1]= l1*(r0-r2);
  H2[3][2][0]= l1*q0-(l2-l0)*q1;      H2[3][2][1]=-2+6*l2+l1*r0-(l2-l0)*r1;
  
  double dldx[3][2];
  double const_term[3];
  TriDlDx(dldx,const_term,coord0,coord1,coord2);
  
  for(unsigned int i=0;i<3;i++){
    B1[0][i] =  dldx[1][0]*H1[2][i]+dldx[2][0]*H1[3][i];
    B1[1][i] = -dldx[1][1]*H1[0][i]-dldx[2][1]*H1[1][i];
    B1[2][i] =  dldx[1][1]*H1[2][i]+dldx[2][1]*H1[3][i] - dldx[1][0]*H1[0][i]-dldx[2][0]*H1[1][i];
  }
  for(unsigned int i=0;i<3;i++){
    B2[0][i][0] =  dldx[1][0]*H2[2][i][0]+dldx[2][0]*H2[3][i][0];
    B2[0][i][1] =  dldx[1][0]*H2[2][i][1]+dldx[2][0]*H2[3][i][1];
    B2[1][i][0] = -dldx[1][1]*H2[0][i][0]-dldx[2][1]*H2[1][i][0];
    B2[1][i][1] = -dldx[1][1]*H2[0][i][1]-dldx[2][1]*H2[1][i][1];
    B2[2][i][0] =  dldx[1][1]*H2[2][i][0]+dldx[2][1]*H2[3][i][0] - dldx[1][0]*H2[0][i][0]-dldx[2][0]*H2[1][i][0];
    B2[2][i][1] =  dldx[1][1]*H2[2][i][1]+dldx[2][1]*H2[3][i][1] - dldx[1][0]*H2[0][i][1]-dldx[2][0]*H2[1][i][1];
  }
}

}
}


/*
const static unsigned int NIntTriGauss[3] = { 1, 3, 7 };
const static double TriGauss[3][7][3] =
{
  { // liner
    { 0.3333333333, 0.3333333333, 1.0 },
    { 0.0, 0.0, 0.0 },
    { 0.0, 0.0, 0.0 },
    { 0.0, 0.0, 0.0 },
    { 0.0, 0.0, 0.0 },
    { 0.0, 0.0, 0.0 },
  },
  { // quadratic
    { 0.1666666667, 0.1666666667, 0.3333333333 },
    { 0.6666666667, 0.1666666667, 0.3333333333 },
    { 0.1666666667, 0.6666666667, 0.3333333333 },
    { 0.0, 0.0, 0.0 },
    { 0.0, 0.0, 0.0 },
    { 0.0, 0.0, 0.0 },
  },
  { // cubic
    { 0.1012865073, 0.1012865073, 0.1259391805 },
    { 0.7974269854, 0.1012865073, 0.1259391805 },
    { 0.1012865073, 0.7974269854, 0.1259391805 },
    { 0.4701420641, 0.0597158718, 0.1323941527 },
    { 0.4701420641, 0.4701420641, 0.1323941527 },
    { 0.0597158718, 0.4701420641, 0.1323941527 },
    { 0.3333333333, 0.3333333333, 0.225 },
  }
};
 */



// compute energy and its 1st and 2nd derivative for cloth bending
DFM2_INLINE void delfem2::WdWddW_Bend
(double& W,  // (out) strain energy
 double dW[4][3], // (out) 1st derivative of energy
 double ddW[4][4][3][3], // (out) 2nd derivative of energy
 ////
 const double C[4][3], // (in) undeformed triangle vertex positions
 const double c[4][3], // (in) deformed triangle vertex positions
 double stiff)
{
  const double A0 = femem3::TriArea3D(C[0],C[2],C[3]);
  const double A1 = femem3::TriArea3D(C[1],C[3],C[2]);
  const double L0 = femem3::Distance3D(C[2],C[3]);
  const double H0 = A0*2.0/L0;
  const double H1 = A1*2.0/L0;
  const double e23[3] = { C[3][0]-C[2][0], C[3][1]-C[2][1], C[3][2]-C[2][2] };
  const double e02[3] = { C[2][0]-C[0][0], C[2][1]-C[0][1], C[2][2]-C[0][2] };
  const double e03[3] = { C[3][0]-C[0][0], C[3][1]-C[0][1], C[3][2]-C[0][2] };
  const double e12[3] = { C[2][0]-C[1][0], C[2][1]-C[1][1], C[2][2]-C[1][2] };
  const double e13[3] = { C[3][0]-C[1][0], C[3][1]-C[1][1], C[3][2]-C[1][2] };
  double cot023, cot032;
  {
    const double r2 = -femem3::Dot3D(e02,e23);
    const double r3 = +femem3::Dot3D(e03,e23);
    cot023 = r2/H0;
    cot032 = r3/H0;
  }
  double cot123, cot132;
  {
    const double r2 = -femem3::Dot3D(e12,e23);
    const double r3 = +femem3::Dot3D(e13,e23);
    cot123 = r2/H1;
    cot132 = r3/H1;
  }
  const double tmp0 = stiff/((A0+A1)*L0*L0);
  const double K[4] = { -cot023-cot032, -cot123-cot132, cot032+cot132, cot023+cot123 };
  
  // compute 2nd derivative of energy
  for(int i=0;i<4*4*3*3;i++){ (&ddW[0][0][0][0])[i] = 0; }
  for(int ino=0;ino<4;ino++){
    for(int jno=0;jno<4;jno++){
      const double tmp = K[ino]*K[jno]*tmp0;
      ddW[ino][jno][0][0] = tmp;
      ddW[ino][jno][1][1] = tmp;
      ddW[ino][jno][2][2] = tmp;
    }
  }
  // compute 1st derivative of energy
  W = 0.0;
  for(int ino=0;ino<4;ino++){
    for(int idim=0;idim<3;idim++){
      dW[ino][idim] = 0;
      for(int jno=0;jno<4;jno++){
        for(int jdim=0;jdim<3;jdim++){
          dW[ino][idim] += ddW[ino][jno][idim][jdim]*c[jno][jdim];
        }
      }
      W += dW[ino][idim]*c[ino][idim];
    }
  }
}

DFM2_INLINE void delfem2::WdWddW_CST
(double& W, // (out) energy
 double dW[3][3], // (out) 1st derivative of energy
 double ddW[3][3][3][3], // (out) 2nd derivative of energy
 ////
 const double C[3][3], // (in) undeformed triangle vertex positions
 const double c[3][3], // (in) deformed triangle vertex positions
 const double lambda, // (in) Lame's 1st parameter
 const double myu)     // (in) Lame's 2nd parameter
{
  double Gd[3][3] = { // undeformed edge vector
    { C[1][0]-C[0][0], C[1][1]-C[0][1], C[1][2]-C[0][2] },
    { C[2][0]-C[0][0], C[2][1]-C[0][1], C[2][2]-C[0][2] }, { 0,0,0 } };
  double Area;
  femem3::UnitNormalAreaTri3D(Gd[2], Area, C[0], C[1], C[2]);
  
  double Gu[2][3]; // inverse of Gd
  {
    femem3::Cross3D(Gu[0], Gd[1], Gd[2]);
    const double invtmp1 = 1.0/femem3::Dot3D(Gu[0],Gd[0]);
    Gu[0][0] *= invtmp1;	Gu[0][1] *= invtmp1;	Gu[0][2] *= invtmp1;
    //
    femem3::Cross3D(Gu[1], Gd[2], Gd[0]);
    const double invtmp2 = 1.0/femem3::Dot3D(Gu[1],Gd[1]);
    Gu[1][0] *= invtmp2;	Gu[1][1] *= invtmp2;	Gu[1][2] *= invtmp2;
  }
  
  const double gd[2][3] = { // deformed edge vector
    { c[1][0]-c[0][0], c[1][1]-c[0][1], c[1][2]-c[0][2] },
    { c[2][0]-c[0][0], c[2][1]-c[0][1], c[2][2]-c[0][2] } };
  
  const double E2[3] = {  // green lagrange strain (with engineer's notation)
    0.5*( femem3::Dot3D(gd[0],gd[0]) - femem3::Dot3D(Gd[0],Gd[0]) ),
    0.5*( femem3::Dot3D(gd[1],gd[1]) - femem3::Dot3D(Gd[1],Gd[1]) ),
    1.0*( femem3::Dot3D(gd[0],gd[1]) - femem3::Dot3D(Gd[0],Gd[1]) ) };
  const double GuGu2[3] = { femem3::Dot3D(Gu[0],Gu[0]), femem3::Dot3D(Gu[1],Gu[1]), femem3::Dot3D(Gu[1],Gu[0]) };
  const double Cons2[3][3] = { // constitutive tensor
    { lambda*GuGu2[0]*GuGu2[0] + 2*myu*(GuGu2[0]*GuGu2[0]),
      lambda*GuGu2[0]*GuGu2[1] + 2*myu*(GuGu2[2]*GuGu2[2]),
      lambda*GuGu2[0]*GuGu2[2] + 2*myu*(GuGu2[0]*GuGu2[2]) },
    { lambda*GuGu2[1]*GuGu2[0] + 2*myu*(GuGu2[2]*GuGu2[2]),
      lambda*GuGu2[1]*GuGu2[1] + 2*myu*(GuGu2[1]*GuGu2[1]),
      lambda*GuGu2[1]*GuGu2[2] + 2*myu*(GuGu2[2]*GuGu2[1]) },
    { lambda*GuGu2[2]*GuGu2[0] + 2*myu*(GuGu2[0]*GuGu2[2]),
      lambda*GuGu2[2]*GuGu2[1] + 2*myu*(GuGu2[2]*GuGu2[1]),
      lambda*GuGu2[2]*GuGu2[2] + 1*myu*(GuGu2[0]*GuGu2[1] + GuGu2[2]*GuGu2[2]) } };
  const double S2[3] = {  // 2nd Piola-Kirchhoff stress
    Cons2[0][0]*E2[0] + Cons2[0][1]*E2[1] + Cons2[0][2]*E2[2],
    Cons2[1][0]*E2[0] + Cons2[1][1]*E2[1] + Cons2[1][2]*E2[2],
    Cons2[2][0]*E2[0] + Cons2[2][1]*E2[1] + Cons2[2][2]*E2[2] };
  
  // compute energy
  W = 0.5*Area*(E2[0]*S2[0] + E2[1]*S2[1] + E2[2]*S2[2]);
  
  // compute 1st derivative
  const double dNdr[3][2] = { {-1.0, -1.0}, {+1.0, +0.0}, {+0.0, +1.0} };
  for(int ino=0;ino<3;ino++){
    for(int idim=0;idim<3;idim++){
      dW[ino][idim] = Area*
      (+S2[0]*gd[0][idim]*dNdr[ino][0]
       +S2[2]*gd[0][idim]*dNdr[ino][1]
       +S2[2]*gd[1][idim]*dNdr[ino][0]
       +S2[1]*gd[1][idim]*dNdr[ino][1]);
    }
  }
  
  double S3[3] = { S2[0], S2[1], S2[2] };
  femem3::MakePositiveDefinite_Sim22(S2,S3);
  
  // compute second derivative
  for(int ino=0;ino<3;ino++){
    for(int jno=0;jno<3;jno++){
      for(int idim=0;idim<3;idim++){
        for(int jdim=0;jdim<3;jdim++){
          double dtmp0 = 0;
          dtmp0 += gd[0][idim]*dNdr[ino][0]*Cons2[0][0]*gd[0][jdim]*dNdr[jno][0];
          dtmp0 += gd[0][idim]*dNdr[ino][0]*Cons2[0][1]*gd[1][jdim]*dNdr[jno][1];
          dtmp0 += gd[0][idim]*dNdr[ino][0]*Cons2[0][2]*gd[0][jdim]*dNdr[jno][1];
          dtmp0 += gd[0][idim]*dNdr[ino][0]*Cons2[0][2]*gd[1][jdim]*dNdr[jno][0];
          dtmp0 += gd[1][idim]*dNdr[ino][1]*Cons2[1][0]*gd[0][jdim]*dNdr[jno][0];
          dtmp0 += gd[1][idim]*dNdr[ino][1]*Cons2[1][1]*gd[1][jdim]*dNdr[jno][1];
          dtmp0 += gd[1][idim]*dNdr[ino][1]*Cons2[1][2]*gd[0][jdim]*dNdr[jno][1];
          dtmp0 += gd[1][idim]*dNdr[ino][1]*Cons2[1][2]*gd[1][jdim]*dNdr[jno][0];
          dtmp0 += gd[0][idim]*dNdr[ino][1]*Cons2[2][0]*gd[0][jdim]*dNdr[jno][0];
          dtmp0 += gd[0][idim]*dNdr[ino][1]*Cons2[2][1]*gd[1][jdim]*dNdr[jno][1];
          dtmp0 += gd[0][idim]*dNdr[ino][1]*Cons2[2][2]*gd[0][jdim]*dNdr[jno][1];
          dtmp0 += gd[0][idim]*dNdr[ino][1]*Cons2[2][2]*gd[1][jdim]*dNdr[jno][0];
          dtmp0 += gd[1][idim]*dNdr[ino][0]*Cons2[2][0]*gd[0][jdim]*dNdr[jno][0];
          dtmp0 += gd[1][idim]*dNdr[ino][0]*Cons2[2][1]*gd[1][jdim]*dNdr[jno][1];
          dtmp0 += gd[1][idim]*dNdr[ino][0]*Cons2[2][2]*gd[0][jdim]*dNdr[jno][1];
          dtmp0 += gd[1][idim]*dNdr[ino][0]*Cons2[2][2]*gd[1][jdim]*dNdr[jno][0];
          ddW[ino][jno][idim][jdim] = dtmp0*Area;
        }
      }
      const double dtmp1 = Area*
      (+S3[0]*dNdr[ino][0]*dNdr[jno][0]
       +S3[2]*dNdr[ino][0]*dNdr[jno][1]
       +S3[2]*dNdr[ino][1]*dNdr[jno][0]
       +S3[1]*dNdr[ino][1]*dNdr[jno][1]);
      ddW[ino][jno][0][0] += dtmp1;
      ddW[ino][jno][1][1] += dtmp1;
      ddW[ino][jno][2][2] += dtmp1;
    }
  }
}



// compute energy and its 1st and 2nd derivative for contact against object
DFM2_INLINE void delfem2::WdWddW_Contact
(double& W,  // (out) energy
 double dW[3], // (out) 1st derivative of energy
 double ddW[3][3], // (out) 2nd derivative of energy
 ////
 const double c[3], // (in) deformed triangle vertex positions
 double stiff_contact,
 double contact_clearance,
 const CInput_Contact& input )
{
  double n[3];
  double pd = input.penetrationNormal(n[0],n[1],n[2], c[0],c[1],c[2]);
  pd += contact_clearance;
  if( pd  < 0 ){
    W = 0;
    dW[0] = 0;  dW[1] = 0;  dW[2] = 0;
    ddW[0][0] = 0;  ddW[0][1] = 0;  ddW[0][2] = 0;
    ddW[1][0] = 0;  ddW[1][1] = 0;  ddW[1][2] = 0;
    ddW[2][0] = 0;  ddW[2][1] = 0;  ddW[2][2] = 0;
    return;
  }
  W = 0.5*stiff_contact*pd*pd;
  
  dW[0] = -stiff_contact*pd*n[0];
  dW[1] = -stiff_contact*pd*n[1];
  dW[2] = -stiff_contact*pd*n[2];
  
  ddW[0][0] = stiff_contact*n[0]*n[0];
  ddW[0][1] = stiff_contact*n[0]*n[1];
  ddW[0][2] = stiff_contact*n[0]*n[2];
  ddW[1][0] = stiff_contact*n[1]*n[0];
  ddW[1][1] = stiff_contact*n[1]*n[1];
  ddW[1][2] = stiff_contact*n[1]*n[2];
  ddW[2][0] = stiff_contact*n[2]*n[0];
  ddW[2][1] = stiff_contact*n[2]*n[1];
  ddW[2][2] = stiff_contact*n[2]*n[2];
}


// --------------------------------------------------------------------------------------


const static unsigned int NIntTetGauss[4] = {
  1, 4, 5, 16
};
const static double TetGauss[4][16][4] =
{
  {	// order-1    1point
    { 0.25, 0.25, 0.25, 1.0 },
  },
  {	// order-2    4point
    { 0.585410196624968, 0.138196601125015, 0.138196601125015, 0.25 },
    { 0.138196601125015, 0.585410196624968, 0.138196601125015, 0.25 },
    { 0.138196601125015, 0.138196601125015, 0.585410196624968, 0.25 },
    { 0.138196601125015, 0.138196601125015, 0.138196601125015, 0.25 },
  },
  {	// order-3    5point
    { 0.25, 0.25, 0.25, -0.8 },
    { 0.5, 0.1666666666666667, 0.1666666666666667, 0.45 },
    { 0.1666666666666667, 0.5, 0.1666666666666667, 0.45 },
    { 0.1666666666666667, 0.1666666666666667, 0.5, 0.45 },
    { 0.1666666666666667, 0.1666666666666667, 0.1666666666666667, 0.45 },
  },
  {	// order-4    16point
    { 0.7716429020672371, 0.07611903264425430, 0.07611903264425430, 0.05037379410012282 },
    { 0.07611903264425430, 0.7716429020672371, 0.07611903264425430, 0.05037379410012282 },
    { 0.07611903264425430, 0.07611903264425430, 0.7716429020672371, 0.05037379410012282 },
    { 0.07611903264425430, 0.07611903264425430, 0.07611903264425430, 0.05037379410012282 },

    { 0.1197005277978019, 0.4042339134672644, 0.4042339134672644, 0.06654206863329239 },
    { 0.4042339134672644, 0.1197005277978019, 0.4042339134672644, 0.06654206863329239 },
    { 0.4042339134672644, 0.4042339134672644, 0.1197005277978019, 0.06654206863329239 },

    { 0.07183164526766925, 0.4042339134672644, 0.4042339134672644, 0.06654206863329239 },
    { 0.4042339134672644, 0.07183164526766925, 0.4042339134672644, 0.06654206863329239 },
    { 0.4042339134672644, 0.4042339134672644, 0.07183164526766925, 0.06654206863329239 },

    { 0.1197005277978019, 0.07183164526766925, 0.4042339134672644, 0.06654206863329239 },
    { 0.4042339134672644, 0.1197005277978019, 0.07183164526766925, 0.06654206863329239 },
    { 0.07183164526766925, 0.4042339134672644, 0.1197005277978019, 0.06654206863329239 },

    { 0.07183164526766925, 0.1197005277978019, 0.4042339134672644, 0.06654206863329239 },
    { 0.4042339134672644, 0.07183164526766925, 0.1197005277978019, 0.06654206863329239 },
    { 0.1197005277978019, 0.4042339134672644, 0.07183164526766925, 0.06654206863329239 },
  }
};


DFM2_INLINE void delfem2::ddW_SolidLinear_Tet3D
(double* eKmat,
 double lambda, double myu,
 double vol, double dldx[4][3],
 bool is_add,
 unsigned int nstride)
{
  if( !is_add ){
    for(unsigned int i=0;i<4*4*nstride*nstride;++i){ eKmat[i] = 0.0; }
  }
  for (int ino = 0; ino<4; ino++){
    for (int jno = 0; jno<4; jno++){
      double* pK = eKmat+(nstride*nstride)*(ino*4+jno);
      pK[0*nstride+0] += vol*(lambda*dldx[ino][0]*dldx[jno][0]+myu*dldx[jno][0]*dldx[ino][0]);
      pK[0*nstride+1] += vol*(lambda*dldx[ino][0]*dldx[jno][1]+myu*dldx[jno][0]*dldx[ino][1]);
      pK[0*nstride+2] += vol*(lambda*dldx[ino][0]*dldx[jno][2]+myu*dldx[jno][0]*dldx[ino][2]);
      pK[1*nstride+0] += vol*(lambda*dldx[ino][1]*dldx[jno][0]+myu*dldx[jno][1]*dldx[ino][0]);
      pK[1*nstride+1] += vol*(lambda*dldx[ino][1]*dldx[jno][1]+myu*dldx[jno][1]*dldx[ino][1]);
      pK[1*nstride+2] += vol*(lambda*dldx[ino][1]*dldx[jno][2]+myu*dldx[jno][1]*dldx[ino][2]);
      pK[2*nstride+0] += vol*(lambda*dldx[ino][2]*dldx[jno][0]+myu*dldx[jno][2]*dldx[ino][0]);
      pK[2*nstride+1] += vol*(lambda*dldx[ino][2]*dldx[jno][1]+myu*dldx[jno][2]*dldx[ino][1]);
      pK[2*nstride+2] += vol*(lambda*dldx[ino][2]*dldx[jno][2]+myu*dldx[jno][2]*dldx[ino][2]);
      const double dtmp1 = dldx[ino][0]*dldx[jno][0]+dldx[ino][1]*dldx[jno][1]+dldx[ino][2]*dldx[jno][2];
      pK[0*nstride+0] += vol*myu*dtmp1;
      pK[1*nstride+1] += vol*myu*dtmp1;
      pK[2*nstride+2] += vol*myu*dtmp1;
    }
  }
}

DFM2_INLINE void delfem2::ddW_MassConsistentVal3D_Tet3D
(double* eMmat,
 double rho, double vol,
 bool is_add,
 unsigned int nstride)
{
  if( !is_add ){
    for(unsigned int i=0;i<4*4*nstride*nstride;++i){ eMmat[i] = 0.0; }
  }
  const double dtmp1 = vol*rho*0.05;
  for(int ino=0;ino<4;ino++){
    for(int jno=0;jno<4;jno++){
      double* pM = eMmat+(nstride*nstride)*(ino*4+jno);
      pM[0*nstride+0] += dtmp1;
      pM[1*nstride+1] += dtmp1;
      pM[2*nstride+2] += dtmp1;
    }
    {
      double* pM = eMmat+(nstride*nstride)*(ino*4+ino);
      pM[0*nstride+0] += dtmp1;
      pM[1*nstride+1] += dtmp1;
      pM[2*nstride+2] += dtmp1;
    }
  }
}

DFM2_INLINE void delfem2::EMat_Poisson_Tet3D
(double eres[4],
 double emat[4][4],
 const double alpha, const double source,
 const double coords[4][3],
 const double value[4])
{
  const int nno = 4;
  const int ndim = 3;
  //
  eres[0] = 0;  eres[1] = 0;  eres[2] = 0;  eres[3] = 0;
  for (int i = 0; i<16; ++i){ (&emat[0][0])[i] = 0.0; }
  const double area = femutil::TetVolume3D(coords[0], coords[1], coords[2], coords[3]);
  //
  double dldx[nno][ndim], const_term[nno];
  TetDlDx(dldx, const_term, coords[0], coords[1], coords[2], coords[3]);
  
  for (int ino = 0; ino<nno; ino++){
    for (int jno = 0; jno<nno; jno++){
      emat[ino][jno] = alpha*area*(dldx[ino][0]*dldx[jno][0]+dldx[ino][1]*dldx[jno][1]+dldx[ino][2]*dldx[jno][2]);
    }
  }
  for (int ino = 0; ino<nno; ino++){
    eres[ino] = source*area*0.25;
  }
  for (int ino = 0; ino<nno; ino++){
    for (int jno = 0; jno<nno; jno++){
      eres[ino] -= emat[ino][jno]*value[jno];
    }
  }
}

DFM2_INLINE void delfem2::EMat_Diffusion_Newmark_Tet3D
(double eres[4],
 double emat[4][4],
 const double alpha, const double source,
 const double dt_timestep, const double gamma_newmark, const double rho,
 const double coords[4][3],
 const double value[4], const double velo[4])
{
  const int nno = 4;
  const int ndim = 3;
  
  eres[0] = 0;  eres[1] = 0;  eres[2] = 0;  eres[3] = 0;
  for (int i=0; i<16; ++i){ (&emat[0][0])[i] = 0.0; }
  
  const double vol = femutil::TetVolume3D(coords[0],coords[1],coords[2],coords[3]);
  double dldx[nno][ndim], const_term[nno];
  TetDlDx(dldx,const_term,coords[0],coords[1],coords[2],coords[3]);
  
  // ----------------------
  
  double eCmat[nno][nno];
  for(int ino=0;ino<nno;ino++){
    for(int jno=0;jno<nno;jno++){
      eCmat[ino][jno] = alpha*vol*(dldx[ino][0]*dldx[jno][0]+dldx[ino][1]*dldx[jno][1]+dldx[ino][2]*dldx[jno][2]);
    }
  }
  double eMmat[nno][nno];
  {
    const double dtmp1 = rho*vol*0.05;
    for(int ino=0;ino<nno;ino++){
      for(int jno=0;jno<nno;jno++){
        eMmat[ino][jno] = dtmp1;
      }
      eMmat[ino][ino] += dtmp1;
    }
  }
		
  for(int ino=0;ino<nno;ino++){
    eres[ino] = source*vol*0.25;
  }
  
  {
    const double dtmp1 = gamma_newmark*dt_timestep;
    for(int i=0;i<nno*nno;i++){
      (&emat[0][0])[i] = (&eMmat[0][0])[i]+dtmp1*(&eCmat[0][0])[i];
    }
  }
  for(int ino=0;ino<nno;ino++){
    for(int jno=0;jno<nno;jno++){
      eres[ino]	-= eCmat[ino][jno]*(value[jno]+dt_timestep*velo[jno]) + eMmat[ino][jno]*velo[jno];
    }
  }
}


DFM2_INLINE void delfem2::stress_LinearSolid_TetP2
(double stress[3][3],
const double l0, const double l1, const double l2, const double l3,
const double vol, const double lambda, const double myu,
const double g_x, const double g_y, const double g_z,
const double dldx[4][3],
const double disp[10][3])
{
  /*
  double N[10] = {
    l0*(2*l0-1),
    l1*(2*l1-1),
    l2*(2*l2-1),
    l3*(2*l3-1),
    4*l0*l1,
    4*l1*l2,
    4*l0*l2,
    4*l0*l3,
    4*l1*l3,
    4*l2*l3,
  };
   */
  double dNdx[10][3];
  for (unsigned int i = 0; i<3; i++){
    dNdx[0][i] = (4*l0-1)*dldx[0][i];
    dNdx[1][i] = (4*l1-1)*dldx[1][i];
    dNdx[2][i] = (4*l2-1)*dldx[2][i];
    dNdx[3][i] = (4*l3-1)*dldx[3][i];
    dNdx[4][i] = 4*dldx[0][i]*l1+4*l0*dldx[1][i];
    dNdx[5][i] = 4*dldx[1][i]*l2+4*l1*dldx[2][i];
    dNdx[6][i] = 4*dldx[0][i]*l2+4*l0*dldx[2][i];
    dNdx[7][i] = 4*dldx[0][i]*l3+4*l0*dldx[3][i];
    dNdx[8][i] = 4*dldx[1][i]*l3+4*l1*dldx[3][i];
    dNdx[9][i] = 4*dldx[2][i]*l3+4*l2*dldx[3][i];
  }
  double dudx[3][3] = { { 0, 0, 0 }, { 0, 0, 0 }, { 0, 0, 0 } };
  for (unsigned int ino = 0; ino<10; ino++){
    for (unsigned int i = 0; i<3; i++){
      for (unsigned int j = 0; j<3; j++){
        dudx[i][j] += disp[ino][i]*dNdx[ino][j];
      }
    }
  }
  double strain[3][3];
  for (unsigned int i = 0; i<3; i++){
    for (unsigned int j = 0; j<3; j++){
      strain[i][j] = 0.5*(dudx[i][j]+dudx[j][i]);
    }
  }
  {
    for (unsigned int i = 0; i<3; i++){
      for (unsigned int j = 0; j<3; j++){
        stress[i][j] = 2*myu*strain[i][j];
      }
      stress[i][i] += lambda*(strain[0][0]+strain[1][1]+strain[2][2]);
    }
  }
}


DFM2_INLINE void delfem2::matRes_LinearSolid_TetP2
(double emat[10][10][3][3],
double eres[10][3],
const double vol, const double lambda, const double myu,
const double g_x, const double g_y, const double g_z,
const double rho,
const double dldx[4][3],
const double disp[10][3])
{
  for (unsigned int i = 0; i<10*10*3*3; i++){ (&emat[0][0][0][0])[i] = 0; }
  for (unsigned int i = 0; i<10*3; i++){ (&eres[0][0])[i] = 0; }
  unsigned int nOrder = 2;
  unsigned int nInt = NIntTetGauss[nOrder];
  for (unsigned int iint = 0; iint<nInt; iint++){
    double l0 = TetGauss[nOrder][iint][0];
    double l1 = TetGauss[nOrder][iint][1];
    double l2 = TetGauss[nOrder][iint][2];
    double l3 = (1-l0-l1-l2);
    double w = TetGauss[nOrder][iint][3];
    double N[10] = {
      l0*(2*l0-1),
      l1*(2*l1-1),
      l2*(2*l2-1),
      l3*(2*l3-1),
      4*l0*l1,
      4*l1*l2,
      4*l0*l2,
      4*l0*l3,
      4*l1*l3,
      4*l2*l3,
    };
    double dNdx[10][3];
    for (unsigned int i = 0; i<3; i++){
      dNdx[0][i] = (4*l0-1)*dldx[0][i];
      dNdx[1][i] = (4*l1-1)*dldx[1][i];
      dNdx[2][i] = (4*l2-1)*dldx[2][i];
      dNdx[3][i] = (4*l3-1)*dldx[3][i];
      dNdx[4][i] = 4*dldx[0][i]*l1+4*l0*dldx[1][i];
      dNdx[5][i] = 4*dldx[1][i]*l2+4*l1*dldx[2][i];
      dNdx[6][i] = 4*dldx[0][i]*l2+4*l0*dldx[2][i];
      dNdx[7][i] = 4*dldx[0][i]*l3+4*l0*dldx[3][i];
      dNdx[8][i] = 4*dldx[1][i]*l3+4*l1*dldx[3][i];
      dNdx[9][i] = 4*dldx[2][i]*l3+4*l2*dldx[3][i];
    }
    /*
    {
    double tN[4] = {0,0,0,0};
    for(unsigned int ino=0;ino<10;ino++){
    tN[0] += N[ino];
    tN[1] += dNdx[ino][0];
    tN[2] += dNdx[ino][1];
    tN[3] += dNdx[ino][2];
    }
    std::cout << tN[0] << "   " << tN[1] << " " << tN[2] << " " << tN[3] << std::endl;
    }
    */
    for (unsigned int ino = 0; ino<10; ino++){
      for (unsigned int jno = 0; jno<10; jno++){
        emat[ino][jno][0][0] += w*vol*(lambda*dNdx[ino][0]*dNdx[jno][0]+myu*dNdx[jno][0]*dNdx[ino][0]);
        emat[ino][jno][0][1] += w*vol*(lambda*dNdx[ino][0]*dNdx[jno][1]+myu*dNdx[jno][0]*dNdx[ino][1]);
        emat[ino][jno][0][2] += w*vol*(lambda*dNdx[ino][0]*dNdx[jno][2]+myu*dNdx[jno][0]*dNdx[ino][2]);
        emat[ino][jno][1][0] += w*vol*(lambda*dNdx[ino][1]*dNdx[jno][0]+myu*dNdx[jno][1]*dNdx[ino][0]);
        emat[ino][jno][1][1] += w*vol*(lambda*dNdx[ino][1]*dNdx[jno][1]+myu*dNdx[jno][1]*dNdx[ino][1]);
        emat[ino][jno][1][2] += w*vol*(lambda*dNdx[ino][1]*dNdx[jno][2]+myu*dNdx[jno][1]*dNdx[ino][2]);
        emat[ino][jno][2][0] += w*vol*(lambda*dNdx[ino][2]*dNdx[jno][0]+myu*dNdx[jno][2]*dNdx[ino][0]);
        emat[ino][jno][2][1] += w*vol*(lambda*dNdx[ino][2]*dNdx[jno][1]+myu*dNdx[jno][2]*dNdx[ino][1]);
        emat[ino][jno][2][2] += w*vol*(lambda*dNdx[ino][2]*dNdx[jno][2]+myu*dNdx[jno][2]*dNdx[ino][2]);
        const double dtmp1 = dNdx[ino][0]*dNdx[jno][0]+dNdx[ino][1]*dNdx[jno][1]+dNdx[ino][2]*dNdx[jno][2];
        emat[ino][jno][0][0] += w*vol*myu*dtmp1;
        emat[ino][jno][1][1] += w*vol*myu*dtmp1;
        emat[ino][jno][2][2] += w*vol*myu*dtmp1;
      }
    }
    for (unsigned int ino = 0; ino<10; ino++){
      eres[ino][0] += w*vol*rho*g_x*N[ino];
      eres[ino][1] += w*vol*rho*g_y*N[ino];
      eres[ino][2] += w*vol*rho*g_z*N[ino];
      for (unsigned int jno = 0; jno<10; jno++){
        eres[ino][0] -= emat[ino][jno][0][0]*disp[jno][0]+emat[ino][jno][0][1]*disp[jno][1]+emat[ino][jno][0][2]*disp[jno][2];
        eres[ino][1] -= emat[ino][jno][1][0]*disp[jno][0]+emat[ino][jno][1][1]*disp[jno][1]+emat[ino][jno][1][2]*disp[jno][2];
        eres[ino][2] -= emat[ino][jno][2][0]*disp[jno][0]+emat[ino][jno][2][1]*disp[jno][1]+emat[ino][jno][2][2]*disp[jno][2];
      }
    }
  }
}


DFM2_INLINE void delfem2::EMat_SolidLinear_Static_Tet
(double emat[4][4][3][3],
 double eres[4][3],
 const double myu, const double lambda,
 const double P[4][3],
 const double disp[4][3],
 bool is_add)
{
  const double vol = femutil::TetVolume3D(P[0], P[1], P[2], P[3]);
  double dldx[4][3];
  {
    double const_term[4];    
    TetDlDx(dldx, const_term, P[0], P[1], P[2], P[3]);
  }
  // ----------------------
  ddW_SolidLinear_Tet3D(&emat[0][0][0][0],
                        lambda, myu, vol, dldx, is_add, 3);
  if( !is_add ){
    for(int i=0;i<12;++i){ (&eres[0][0])[i] = 0.0; }
  }
  for (int ino = 0; ino<4; ino++){
    for (int jno = 0; jno<4; jno++){
      eres[ino][0] -= emat[ino][jno][0][0]*disp[jno][0]+emat[ino][jno][0][1]*disp[jno][1]+emat[ino][jno][0][2]*disp[jno][2];
      eres[ino][1] -= emat[ino][jno][1][0]*disp[jno][0]+emat[ino][jno][1][1]*disp[jno][1]+emat[ino][jno][1][2]*disp[jno][2];
      eres[ino][2] -= emat[ino][jno][2][0]*disp[jno][0]+emat[ino][jno][2][1]*disp[jno][1]+emat[ino][jno][2][2]*disp[jno][2];
    }
  }
}

DFM2_INLINE void delfem2::MakeMat_LinearSolid3D_Static_Q1
(const double myu, const double lambda,
 const double rho, const double g_x, const double g_y, const double g_z,
 const double coords[8][3],
 const double disp[8][3],
 //
 double emat[8][8][3][3],
 double eres[8][3])
{
  namespace lcl = ::delfem2::femem3;
  const int nDegInt = 2;
  const int nInt = lcl::NIntLineGauss[nDegInt];
  const double (*Gauss)[2] = lcl::LineGauss[nDegInt];
  
  for(unsigned int i=0;i<8*8*3*3;i++){ *( &emat[0][0][0][0]+i) = 0.0; }
  for(unsigned int i=0;i<    8*3;i++){ *( &eres[0][0]      +i) = 0.0; }
  
  double vol = 0.0;
  for(int ir1=0;ir1<nInt;ir1++){
  for(int ir2=0;ir2<nInt;ir2++){
  for(int ir3=0;ir3<nInt;ir3++){
    const double r1 = Gauss[ir1][0];
    const double r2 = Gauss[ir2][0];
    const double r3 = Gauss[ir3][0];
    double detjac, detwei, dndx[8][3], an[8];
    lcl::ShapeFunc_Hex8(r1,r2,r3,coords,detjac,dndx,an);
    detwei = detjac*Gauss[ir1][1]*Gauss[ir2][1]*Gauss[ir3][1];
    vol += detwei;
    for(int ino=0;ino<8;ino++){
    for(int jno=0;jno<8;jno++){
      double dtmp1 = 0.0;
      for(int idim=0;idim<3;idim++){
        for(int jdim=0;jdim<3;jdim++){
          emat[ino][jno][idim][jdim] += detwei*( lambda*dndx[ino][idim]*dndx[jno][jdim]
                                                +myu*dndx[jno][idim]*dndx[ino][jdim] );
        }
        dtmp1 += dndx[ino][idim]*dndx[jno][idim];
      }
      for(int idim=0;idim<3;idim++){
        emat[ino][jno][idim][idim] += detwei*myu*dtmp1;
      }
    }
    }
    for(int ino=0;ino<8;ino++){
      eres[ino][0] += detwei*rho*g_x*an[ino];
      eres[ino][1] += detwei*rho*g_y*an[ino];
      eres[ino][2] += detwei*rho*g_z*an[ino];
    }
  }
  }
  }
  for (int ino = 0; ino<8; ino++){
  for (int jno = 0; jno<8; jno++){
    eres[ino][0] -= emat[ino][jno][0][0]*disp[jno][0]+emat[ino][jno][0][1]*disp[jno][1]+emat[ino][jno][0][2]*disp[jno][2];
    eres[ino][1] -= emat[ino][jno][1][0]*disp[jno][0]+emat[ino][jno][1][1]*disp[jno][1]+emat[ino][jno][1][2]*disp[jno][2];
    eres[ino][2] -= emat[ino][jno][2][0]*disp[jno][0]+emat[ino][jno][2][1]*disp[jno][1]+emat[ino][jno][2][2]*disp[jno][2];
  }
  }
}



DFM2_INLINE void delfem2::EMat_SolidLinear_NewmarkBeta_MeshTet3D(
    double eres[4][3],
    double emat[4][4][3][3],
    const double myu, const double lambda,
    const double rho, const double g_x, const double g_y, const double g_z,
    const double dt, const double gamma_newmark,  const double beta_newmark,
    const double disp[4][3], const double velo[4][3], const double acc[4][3],
    const double P[4][3],
    bool is_initial_iter)
{
  const int nno = 4;
  const int ndim = 3;
  
  const double vol = femutil::TetVolume3D(P[0],P[1],P[2],P[3]);
  double dldx[nno][ndim];		// spatial derivative of linear shape function
  {
    double zero_order_term[nno];	// const term of shape function
    TetDlDx(dldx, zero_order_term,   P[0],P[1],P[2],P[3]);
  }
  
  double eKmat[nno][nno][ndim][ndim];
  ddW_SolidLinear_Tet3D(&eKmat[0][0][0][0],
                        lambda,myu,
                        vol, dldx, false, 3);
  
  double eMmat[nno][nno][ndim][ndim];
  ddW_MassConsistentVal3D_Tet3D(&eMmat[0][0][0][0],
                                rho,vol,false,3);
  
  // calc external force
  for(int ino=0;ino<nno;ino++){
    eres[ino][0] = vol*rho*g_x*0.25;
    eres[ino][1] = vol*rho*g_y*0.25;
    eres[ino][2] = vol*rho*g_z*0.25;
  }
  
  ////////////////////////////////////////////////////////////////////////////////////
  
  {	// calc coeff matrix for newmark-beta
    double dtmp1 = beta_newmark*dt*dt;
    for(int i=0;i<nno*nno*ndim*ndim;i++){
      (&emat[0][0][0][0])[i] = (&eMmat[0][0][0][0])[i]+dtmp1*(&eKmat[0][0][0][0])[i];
    }
  }
  
  // calc element redisual vector
  for(int ino=0;ino<nno;ino++){
    for(int jno=0;jno<nno;jno++){
      eres[ino][0] -= eKmat[ino][jno][0][0]*disp[jno][0]+eKmat[ino][jno][0][1]*disp[jno][1]+eKmat[ino][jno][0][2]*disp[jno][2];
      eres[ino][1] -= eKmat[ino][jno][1][0]*disp[jno][0]+eKmat[ino][jno][1][1]*disp[jno][1]+eKmat[ino][jno][1][2]*disp[jno][2];
      eres[ino][2] -= eKmat[ino][jno][2][0]*disp[jno][0]+eKmat[ino][jno][2][1]*disp[jno][1]+eKmat[ino][jno][2][2]*disp[jno][2];
    }
    for(int jno=0;jno<nno;jno++){
      eres[ino][0] -= eMmat[ino][jno][0][0]*acc[jno][0]+eMmat[ino][jno][0][1]*acc[jno][1]+eMmat[ino][jno][0][2]*acc[jno][2];
      eres[ino][1] -= eMmat[ino][jno][1][0]*acc[jno][0]+eMmat[ino][jno][1][1]*acc[jno][1]+eMmat[ino][jno][1][2]*acc[jno][2];
      eres[ino][2] -= eMmat[ino][jno][2][0]*acc[jno][0]+eMmat[ino][jno][2][1]*acc[jno][1]+eMmat[ino][jno][2][2]*acc[jno][2];
    }
  }
  if( is_initial_iter ){
    for(int ino=0;ino<nno;ino++){
      for(int jno=0;jno<nno;jno++){
        eres[ino][0] -= dt*(eKmat[ino][jno][0][0]*velo[jno][0]+eKmat[ino][jno][0][1]*velo[jno][1]+eKmat[ino][jno][0][2]*velo[jno][2]);
        eres[ino][1] -= dt*(eKmat[ino][jno][1][0]*velo[jno][0]+eKmat[ino][jno][1][1]*velo[jno][1]+eKmat[ino][jno][1][2]*velo[jno][2]);
        eres[ino][2] -= dt*(eKmat[ino][jno][2][0]*velo[jno][0]+eKmat[ino][jno][2][1]*velo[jno][1]+eKmat[ino][jno][2][2]*velo[jno][2]);
      }
      for(int jno=0;jno<nno;jno++){
        eres[ino][0] -= 0.5*dt*dt*(eKmat[ino][jno][0][0]*acc[jno][0]+eKmat[ino][jno][0][1]*acc[jno][1]+eKmat[ino][jno][0][2]*acc[jno][2]);
        eres[ino][1] -= 0.5*dt*dt*(eKmat[ino][jno][1][0]*acc[jno][0]+eKmat[ino][jno][1][1]*acc[jno][1]+eKmat[ino][jno][1][2]*acc[jno][2]);
        eres[ino][2] -= 0.5*dt*dt*(eKmat[ino][jno][2][0]*acc[jno][0]+eKmat[ino][jno][2][1]*acc[jno][1]+eKmat[ino][jno][2][2]*acc[jno][2]);
      }
    }
  }
}

DFM2_INLINE void delfem2::MakeMat_Stokes3D_Static_P1P1
(double alpha, double g_x, double g_y, double g_z,
 const double coords[4][3],
 const double velo[4][3], const double press[4],
 double emat_uu[4][4][3][3], double emat_up[4][4][3], double emat_pu[4][4][3], double emat_pp[4][4],
 double eres_u[4][3], double eres_p[4])
{
  const unsigned int nno = 4;
  const unsigned int ndim = 3;
  
  const double vol = femutil::TetVolume3D(coords[0],coords[1],coords[2],coords[3]);
  
  double dldx[nno][ndim], const_term[nno];
  TetDlDx(dldx, const_term,   coords[0], coords[1], coords[2], coords[3]);
  
  // --------------------------------------------------
  
  for(unsigned int i=0;i<nno*nno*ndim*ndim;i++){ *(&emat_uu[0][0][0][0]+i) = 0.0; }
  for(unsigned int ino=0;ino<nno;ino++){
    for(unsigned int jno=0;jno<nno;jno++){
      const double dtmp1 = vol*alpha*(dldx[jno][0]*dldx[ino][0]+dldx[jno][1]*dldx[ino][1]+dldx[jno][2]*dldx[ino][2]);
      emat_uu[ino][jno][0][0] = dtmp1;
      emat_uu[ino][jno][1][1] = dtmp1;
      emat_uu[ino][jno][2][2] = dtmp1;
    }
  }
  for(unsigned int i=0;i<nno*nno*ndim;i++){ *(&emat_up[0][0][0]+i) = 0.0; }
  for(unsigned int ino=0;ino<nno;ino++){
    for(unsigned int jno=0;jno<nno;jno++){
      emat_up[ino][jno][0] += vol*dldx[ino][0]*0.25;
      emat_up[ino][jno][1] += vol*dldx[ino][1]*0.25;
      emat_up[ino][jno][2] += vol*dldx[ino][2]*0.25;
    }
  }
  for(unsigned int i=0;i<nno*nno*ndim;i++){ *(&emat_pu[0][0][0]+i) = 0.0; }
  for(unsigned int ino=0;ino<nno;ino++){
    for(unsigned int jno=0;jno<nno;jno++){
      emat_pu[ino][jno][0] += vol*dldx[jno][0]*0.25;
      emat_pu[ino][jno][1] += vol*dldx[jno][1]*0.25;
      emat_pu[ino][jno][2] += vol*dldx[jno][2]*0.25;
    }
  }
  for(unsigned int i=0;i<nno*nno;i++){ *(&emat_pp[0][0]+i) = 0.0; }
  double tau; // relaxation parameter
  {
    const double h = pow(vol/3.14, 0.3333333333)*2;
    tau = -h*h/alpha*0.1;
  }
  for(unsigned int ino=0;ino<nno;ino++){
    for(unsigned int jno=0;jno<nno;jno++){
      emat_pp[ino][jno] = vol*tau*(dldx[jno][0]*dldx[ino][0]+dldx[jno][1]*dldx[ino][1]+dldx[jno][2]*dldx[ino][2]);
    }
  }
  
  for(unsigned int ino=0;ino<nno;ino++){
    eres_u[ino][0] = vol*g_x*0.25;
    eres_u[ino][1] = vol*g_y*0.25;
    eres_u[ino][2] = vol*g_z*0.25;
  }
  eres_p[0] = 0;
  eres_p[1] = 0;
  eres_p[2] = 0;
  eres_p[3] = 0;
  
  ////////////////////////////////////////////////////////////////////////////////////
  
  for(unsigned int ino=0;ino<nno;ino++){
    for(unsigned int jno=0;jno<nno;jno++){
      eres_u[ino][0] -= emat_uu[ino][jno][0][0]*velo[jno][0]+emat_uu[ino][jno][0][1]*velo[jno][1]+emat_uu[ino][jno][0][2]*velo[jno][2];
      eres_u[ino][1] -= emat_uu[ino][jno][1][0]*velo[jno][0]+emat_uu[ino][jno][1][1]*velo[jno][1]+emat_uu[ino][jno][1][2]*velo[jno][2];
      eres_u[ino][2] -= emat_uu[ino][jno][2][0]*velo[jno][0]+emat_uu[ino][jno][2][1]*velo[jno][1]+emat_uu[ino][jno][2][2]*velo[jno][2];
    }
    for(unsigned int jno=0;jno<nno;jno++){
      eres_u[ino][0] -= emat_up[ino][jno][0]*press[jno];
      eres_u[ino][1] -= emat_up[ino][jno][1]*press[jno];
      eres_u[ino][2] -= emat_up[ino][jno][2]*press[jno];
    }
  }
  for(unsigned int ino=0;ino<nno;ino++){
    eres_p[ino] = 0.0;
    for(unsigned int jno=0;jno<nno;jno++){
      eres_p[ino] -= emat_pu[ino][jno][0]*velo[jno][0]+emat_pu[ino][jno][1]*velo[jno][1]+emat_pu[ino][jno][2]*velo[jno][2];
    }
    for(unsigned int jno=0;jno<nno;jno++){
      eres_p[ino] -= emat_pp[ino][jno]*press[jno];
    }
  }
}



DFM2_INLINE void delfem2::MakeMat_Stokes3D_Static_P1
(double alpha, double g_x, double g_y, double g_z,
 const double coords[4][3],
 const double velo_press[4][4],
 double emat[4][4][4][4],
 double eres[4][4])
{
  const int nno = 4;
  
  const double velo[4][3] = {
    {velo_press[0][0],velo_press[0][1],velo_press[0][2]},
    {velo_press[1][0],velo_press[1][1],velo_press[1][2]},
    {velo_press[2][0],velo_press[2][1],velo_press[2][2]},
    {velo_press[3][0],velo_press[3][1],velo_press[3][2]} };
  const double press[4] = { velo_press[0][3], velo_press[1][3], velo_press[2][3], velo_press[3][3] };
  ////
  double emat_uu[4][4][3][3], emat_up[4][4][3], emat_pu[4][4][3], emat_pp[4][4];
  double eres_u[4][3], eres_p[4];
  MakeMat_Stokes3D_Static_P1P1(alpha, g_x, g_y, g_z,
                               coords, velo, press,
                               emat_uu, emat_up, emat_pu, emat_pp,
                               eres_u, eres_p);
  ////
  for(int ino=0;ino<nno;ino++){
    for(int jno=0;jno<nno;jno++){
      emat[ino][jno][0][0] = emat_uu[ino][jno][0][0];
      emat[ino][jno][0][1] = emat_uu[ino][jno][0][1];
      emat[ino][jno][0][2] = emat_uu[ino][jno][0][2];
      ////
      emat[ino][jno][1][0] = emat_uu[ino][jno][1][0];
      emat[ino][jno][1][1] = emat_uu[ino][jno][1][1];
      emat[ino][jno][1][2] = emat_uu[ino][jno][1][2];
      ////
      emat[ino][jno][2][0] = emat_uu[ino][jno][2][0];
      emat[ino][jno][2][1] = emat_uu[ino][jno][2][1];
      emat[ino][jno][2][2] = emat_uu[ino][jno][2][2];
      ////
      emat[ino][jno][0][3] = emat_up[ino][jno][0];
      emat[ino][jno][1][3] = emat_up[ino][jno][1];
      emat[ino][jno][2][3] = emat_up[ino][jno][2];
      ////
      emat[ino][jno][3][0] = emat_pu[ino][jno][0];
      emat[ino][jno][3][1] = emat_pu[ino][jno][1];
      emat[ino][jno][3][2] = emat_pu[ino][jno][2];
      ////
      emat[ino][jno][3][3] = emat_pp[ino][jno];
    }
  }
  for(int ino=0;ino<nno;ino++){
    eres[ino][0] = eres_u[ino][0];
    eres[ino][1] = eres_u[ino][1];
    eres[ino][2] = eres_u[ino][2];
    eres[ino][3] = eres_p[ino];
  }
}



DFM2_INLINE void delfem2::MakeMat_Stokes3D_Dynamic_Newmark_P1P1(
    double alpha, double rho, double g_x, double g_y, double g_z,
    const double dt_timestep, const double gamma_newmark,
    const double coords[4][3],
    const double velo[4][3], const double press[4], const double acc[4][3], const double apress[4],
    double emat_uu[4][4][3][3], double emat_up[4][4][3], double emat_pu[4][4][3], double emat_pp[4][4],
    double eres_u[4][3], double eres_p[4])
{
  //	std::cout << "AddMat_Stokes2D_NonStatic_Newmark_P1P1" << std::endl;
  
  const int nno = 4;
  const int ndim = 3;
  ////
  const double vol = femutil::TetVolume3D(coords[0],coords[1],coords[2],coords[3]);
  double dldx[nno][ndim], const_term[nno];
  TetDlDx(dldx, const_term,   coords[0],coords[1],coords[2],coords[3]);
  
  // ------------------------------
  double eCmat_uu[4][4][3][3];
  for(int ino=0;ino<nno;ino++){
    for(int jno=0;jno<nno;jno++){
      const double dtmp1 = vol*alpha*(dldx[jno][0]*dldx[ino][0]+dldx[jno][1]*dldx[ino][1]+dldx[jno][2]*dldx[ino][2]);
      eCmat_uu[ino][jno][0][0] = dtmp1;
      eCmat_uu[ino][jno][0][1] = 0.0;
      eCmat_uu[ino][jno][0][2] = 0.0;
      ////
      eCmat_uu[ino][jno][1][0] = 0.0;
      eCmat_uu[ino][jno][1][1] = dtmp1;
      eCmat_uu[ino][jno][1][2] = 0.0;
      ////
      eCmat_uu[ino][jno][2][0] = 0.0;
      eCmat_uu[ino][jno][2][1] = 0.0;
      eCmat_uu[ino][jno][2][2] = dtmp1;
    }
  }
  
  double eCmat_up[4][4][3];
  for(int ino=0;ino<nno;ino++){
    for(int jno=0;jno<nno;jno++){
      eCmat_up[ino][jno][0] = vol*dldx[ino][0]*0.25;
      eCmat_up[ino][jno][1] = vol*dldx[ino][1]*0.25;
      eCmat_up[ino][jno][2] = vol*dldx[ino][2]*0.25;
    }
  }
  
  double eCmat_pu[4][4][3];
  for(int ino=0;ino<nno;ino++){
    for(int jno=0;jno<nno;jno++){
      eCmat_pu[ino][jno][0] = vol*dldx[jno][0]*0.25;
      eCmat_pu[ino][jno][1] = vol*dldx[jno][1]*0.25;
      eCmat_pu[ino][jno][2] = vol*dldx[jno][2]*0.25;
    }
  }
  
  double tau;
  {
    const double h = pow( vol / 3.14, 0.3333333 )*2;
    tau = -h*h/alpha*0.1;
  }
  
  double eCmat_pp[4][4];
  for(int ino=0;ino<nno;ino++){
    for(int jno=0;jno<nno;jno++){
      eCmat_pp[ino][jno] = vol*tau*(dldx[jno][0]*dldx[ino][0]+dldx[jno][1]*dldx[ino][1]+dldx[jno][2]*dldx[ino][2]);
    }
  }
  
  // -------------------------
  double eMmat_uu[4][4][3][3];
  {
    const double dtmp1 = vol*rho*0.05;
    for(int ino=0;ino<nno;ino++){
      for(int jno=0;jno<nno;jno++){
        eMmat_uu[ino][jno][0][0] = dtmp1;
        eMmat_uu[ino][jno][0][1] = 0.0;
        eMmat_uu[ino][jno][0][2] = 0.0;
        ////
        eMmat_uu[ino][jno][1][0] = 0.0;
        eMmat_uu[ino][jno][1][1] = dtmp1;
        eMmat_uu[ino][jno][1][2] = 0.0;
        ////
        eMmat_uu[ino][jno][2][0] = 0.0;
        eMmat_uu[ino][jno][2][1] = 0.0;
        eMmat_uu[ino][jno][2][2] = dtmp1;
      }
      eMmat_uu[ino][ino][0][0] += dtmp1;
      eMmat_uu[ino][ino][1][1] += dtmp1;
      eMmat_uu[ino][ino][2][2] += dtmp1;
    }
  }
  
  for(int ino=0;ino<nno;ino++){
    eres_u[ino][0] = vol*g_x*0.25;
    eres_u[ino][1] = vol*g_y*0.25;
    eres_u[ino][2] = vol*g_z*0.25;
  }
  eres_p[0] = 0;
  eres_p[1] = 0;
  eres_p[2] = 0;
  eres_p[3] = 0;
  
  // --------------------------------
  {
    double dtmp1 = gamma_newmark*dt_timestep;
    for(int ino=0;ino<nno;ino++){
      for(int jno=0;jno<nno;jno++){
        for(int idim=0;idim<ndim;idim++){
          for(int jdim=0;jdim<ndim;jdim++){
            emat_uu[ino][jno][idim][jdim] = eMmat_uu[ino][jno][idim][jdim]+dtmp1*eCmat_uu[ino][jno][idim][jdim];
          }
        }
      }
    }
    
    for(int ino=0;ino<nno;ino++){
      for(int jno=0;jno<nno;jno++){
        for(int idim=0;idim<ndim;idim++){
          emat_up[ino][jno][idim] = dtmp1*eCmat_up[ino][jno][idim];
          emat_pu[ino][jno][idim] = dtmp1*eCmat_pu[ino][jno][idim];
        }
      }
    }
    
    for(int ino=0;ino<nno;ino++){
      for(int jno=0;jno<nno;jno++){
        emat_pp[ino][jno] = dtmp1*eCmat_pp[ino][jno];
      }
    }
  }
  
  for(int ino=0;ino<nno;ino++){
    for(int idim=0;idim<ndim;idim++){
      for(int jno=0;jno<nno;jno++){
        for(int jdim=0;jdim<ndim;jdim++){
          eres_u[ino][idim] -= eCmat_uu[ino][jno][idim][jdim]*(velo[jno][jdim]+dt_timestep*acc[jno][jdim])
          + eMmat_uu[ino][jno][idim][jdim]*acc[jno][jdim];
        }
      }
      for(int jno=0;jno<nno;jno++){
        eres_u[ino][idim] -= eCmat_up[ino][jno][idim]*(press[jno]+dt_timestep*apress[jno]);
      }
    }
  }
  
  for(int ino=0;ino<nno;ino++){
    for(int jno=0;jno<nno;jno++){
      for(int jdim=0;jdim<ndim;jdim++){
        eres_p[ino] -= eCmat_pu[ino][jno][jdim]*(velo[jno][jdim]+dt_timestep*acc[jno][jdim]);
      }
    }
    for(int jno=0;jno<nno;jno++){
      eres_p[ino] -= eCmat_pp[ino][jno]*(press[jno]+dt_timestep*apress[jno]);
    }
  }
}


DFM2_INLINE void delfem2::MakeMat_Stokes3D_Dynamic_P1
(double alpha, double rho, double g_x, double g_y, double g_z,
 const double dt_timestep, const double gamma_newmark,
 const double coords[4][3],
 const double velo_press[4][4], const double acc_apress[4][4],
 double emat[4][4][4][4],
 double eres[4][4])
{
  const int nno = 4;
  const int ndim = 4;
  //
  const double vol = femutil::TetVolume3D(coords[0],coords[1],coords[2],coords[3]);
  double dldx[nno][3], const_term[nno];
  TetDlDx(dldx, const_term,   coords[0], coords[1], coords[2], coords[3]);
  //
  double tau;
  {
    const double h = pow( vol / 3.14, 0.333333333 )*2;
    tau = -h*h/alpha*0.1;
  }
  // -------------------------
  double eCmat[4][4][4][4];
  for(int ino=0;ino<nno;ino++){
    for(int jno=0;jno<nno;jno++){
      const double dtmp1 = vol*alpha*(dldx[jno][0]*dldx[ino][0]+dldx[jno][1]*dldx[ino][1]+dldx[jno][2]*dldx[ino][2]);
      eCmat[ino][jno][0][0] = dtmp1;
      eCmat[ino][jno][1][1] = dtmp1;
      eCmat[ino][jno][2][2] = dtmp1;
      eCmat[ino][jno][0][1] = 0.0;
      eCmat[ino][jno][0][2] = 0.0;
      eCmat[ino][jno][1][0] = 0.0;
      eCmat[ino][jno][1][2] = 0.0;
      eCmat[ino][jno][2][0] = 0.0;
      eCmat[ino][jno][2][1] = 0.0;
      eCmat[ino][jno][3][0] = vol*dldx[jno][0]*0.25;
      eCmat[ino][jno][3][1] = vol*dldx[jno][1]*0.25;
      eCmat[ino][jno][3][2] = vol*dldx[jno][2]*0.25;
      eCmat[ino][jno][0][3] = vol*dldx[ino][0]*0.25;
      eCmat[ino][jno][1][3] = vol*dldx[ino][1]*0.25;
      eCmat[ino][jno][2][3] = vol*dldx[ino][2]*0.25;
      eCmat[ino][jno][3][3] = vol*tau*(dldx[jno][0]*dldx[ino][0]+dldx[jno][1]*dldx[ino][1]+dldx[jno][2]*dldx[ino][2]);
    }
  }
  
  ////////////////
  double eMmat[4][4][4][4];
  ddW_MassConsistentVal3D_Tet3D(&eMmat[0][0][0][0], rho, vol, false, 4);
  
  for(int ino=0;ino<nno;ino++){
    eres[ino][0] = vol*g_x*0.25;
    eres[ino][1] = vol*g_y*0.25;
    eres[ino][2] = vol*g_z*0.25;
    eres[ino][3] = 0.0;
  }
  
  ////////////////////////////////
  
  {
    double dtmp1 = gamma_newmark*dt_timestep;
    for(int i=0;i<nno*nno*ndim*ndim;++i){
      (&emat[0][0][0][0])[i] = (&eMmat[0][0][0][0])[i]+dtmp1*(&eCmat[0][0][0][0])[i];
    }
  }
  for(int ino=0;ino<nno;ino++){
    for(int idim=0;idim<ndim;idim++){
      for(int jno=0;jno<nno;jno++){
        for(int jdim=0;jdim<ndim;jdim++){
          eres[ino][idim]
          -= eCmat[ino][jno][idim][jdim]*(velo_press[jno][jdim]+dt_timestep*acc_apress[jno][jdim])
          + eMmat[ino][jno][idim][jdim]*acc_apress[jno][jdim];
        }
      }
    }
  }
}

// ------------------------------


DFM2_INLINE void delfem2::MakeMat_PlateBendingDKT
(double emat_ww[3][3],
 double emat_wr[3][3][2],
 double emat_rw[3][3][2],
 double emat_rr[3][3][2][2],
 double eres_w[3],
 double eres_r[3][2],
 const double young, const double poisson, const double thickness,
 const double coord[][2], const double w[], const double rot[][2])
{
  namespace lcl = ::delfem2::femem3;
  const unsigned int ndim = 2;
  const unsigned int nno = 3;
  
  for(unsigned int i=0;i<nno*nno;    i++){ *(&emat_ww[0][0]      +i) = 0.0; }
  for(unsigned int i=0;i<nno*nno*2;  i++){ *(&emat_wr[0][0][0]   +i) = 0.0; }
  for(unsigned int i=0;i<nno*nno*2;  i++){ *(&emat_rw[0][0][0]   +i) = 0.0; }
  for(unsigned int i=0;i<nno*nno*2*2;i++){ *(&emat_rr[0][0][0][0]+i) = 0.0; }
  
  double dmat[3][3];
  {
    const double dtmp1 = young*thickness*thickness*thickness/(12.0*(1.0-poisson*poisson));
    dmat[0][0] = dtmp1;      dmat[0][1] = dtmp1*poisson;  dmat[0][2] = 0;
    dmat[1][0] = dtmp1*poisson;  dmat[1][1] = dtmp1;      dmat[1][2] = 0;
    dmat[2][0] = 0.0;      dmat[2][1] = 0.0;      dmat[2][2] = dtmp1*(1-poisson)*0.5;
  }
  double B1[3][nno];
  double B2[3][nno][ndim];
  const double area = lcl::TriArea2D(coord[0],coord[1],coord[2]);
  const double dtmp1 = area/3.0;
  for(unsigned int iw=0;iw<3;iw++){
    if(      iw == 0 ){ femem3::MakeCurvetureDKT(B1,B2,coord[0],coord[1],coord[2],0.5,0.5); }
    else if( iw == 1 ){ femem3::MakeCurvetureDKT(B1,B2,coord[0],coord[1],coord[2],0.0,0.5); }
    else if( iw == 2 ){ femem3::MakeCurvetureDKT(B1,B2,coord[0],coord[1],coord[2],0.5,0.0); }
    for(unsigned int ino=0;ino<nno;ino++){
      for(unsigned int jno=0;jno<nno;jno++){
        for(unsigned int k=0;k<3;k++){
          for(unsigned int l=0;l<3;l++){
            emat_ww[ino][jno] += dtmp1*B1[k][ino]*dmat[k][l]*B1[l][jno];
            for(unsigned int idim=0;idim<ndim;idim++){
              for(unsigned int jdim=0;jdim<ndim;jdim++){
                emat_rr[ino][jno][idim][jdim] += dtmp1*B2[k][ino][idim]*dmat[k][l]*B2[l][jno][jdim];
              }
            }
            for(unsigned int idim=0;idim<ndim;idim++){
              emat_rw[ino][jno][idim] += dtmp1*B2[k][ino][idim]*dmat[k][l]*B1[l][jno];
              emat_wr[ino][jno][idim] += dtmp1*B1[k][ino]*dmat[k][l]*B2[l][jno][idim];
            }
          }
        }
      }
    }
  }
  
  for(unsigned int ino=0;ino<nno;ino++){
    for(unsigned int idim=0;idim<2;idim++){
      eres_r[ino][idim] = 0.0;
      for(unsigned int jno=0;jno<nno;jno++){
        eres_r[ino][idim] -= emat_rr[ino][jno][idim][0]*rot[jno][0]
        + emat_rr[ino][jno][idim][1]*rot[jno][1]
        + emat_rw[ino][jno][idim]*w[jno];
      }
    }
  }
  for(unsigned int ino=0;ino<nno;ino++){
    eres_w[ino] = 0.0;
    for(unsigned int jno=0;jno<nno;jno++){
      eres_w[ino] -= emat_ww[ino][jno]*w[jno]
      + emat_wr[ino][jno][0]*rot[jno][0]
      + emat_wr[ino][jno][1]*rot[jno][1];
    }
  }
}