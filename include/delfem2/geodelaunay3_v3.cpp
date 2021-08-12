/*
 * Copyright (c) 2019 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */

#include <cstdlib>
#include <cmath>
#include <stack>

#include "delfem2/geodelaunay3_v3.h"

#ifndef M_PI
#  define M_PI 3.14159265358979323846
#endif

namespace delfem2{
namespace delaunay3{

template <typename T>
T Volume_Tet3(
    const T v1[3],
    const T v2[3],
    const T v3[3],
    const T v4[3])
{
  return
      ((v2[0]-v1[0])*((v3[1]-v1[1])*(v4[2]-v1[2])-(v4[1]-v1[1])*(v3[2]-v1[2]))
       -(v2[1]-v1[1])*((v3[0]-v1[0])*(v4[2]-v1[2])-(v4[0]-v1[0])*(v3[2]-v1[2]))
       +(v2[2]-v1[2])*((v3[0]-v1[0])*(v4[1]-v1[1])-(v4[0]-v1[0])*(v3[1]-v1[1]))
      ) * 0.16666666666666666666666666666667;
}

}
}

template <typename T>
double delfem2::SquareCircumradius(
 const CVec3<T>& ipo0,
 const CVec3<T>& ipo1,
 const CVec3<T>& ipo2,
 const CVec3<T>& ipo3)
{
  double base[3][3] = {
    { ipo1.p[0]-ipo0.p[0], ipo1.p[1]-ipo0.p[1], ipo1.p[2]-ipo0.p[2] },
    { ipo2.p[0]-ipo0.p[0], ipo2.p[1]-ipo0.p[1], ipo2.p[2]-ipo0.p[2] },
    { ipo3.p[0]-ipo0.p[0], ipo3.p[1]-ipo0.p[1], ipo3.p[2]-ipo0.p[2] }
  };
  double s[6] = {
    base[0][0]*base[0][0]+base[0][1]*base[0][1]+base[0][2]*base[0][2],
    base[1][0]*base[1][0]+base[1][1]*base[1][1]+base[1][2]*base[1][2],
    base[2][0]*base[2][0]+base[2][1]*base[2][1]+base[2][2]*base[2][2],
    base[1][0]*base[2][0]+base[1][1]*base[2][1]+base[1][2]*base[2][2],
    base[2][0]*base[0][0]+base[2][1]*base[0][1]+base[2][2]*base[0][2],
    base[0][0]*base[1][0]+base[0][1]*base[1][1]+base[0][2]*base[1][2],
  };
  const double vol = Volume_Tet(ipo0,ipo1,ipo2,ipo3)*6.0;
  if( vol < 1.0e-20 ){ assert(0); }
  const double inv_det = 1.0 / (vol*vol);
  double t[6] = {
    (s[1]*s[2]-s[3]*s[3])*0.5*inv_det,
    (s[2]*s[0]-s[4]*s[4])*0.5*inv_det,
    (s[0]*s[1]-s[5]*s[5])*0.5*inv_det,
    (s[4]*s[5]-s[0]*s[3])*0.5*inv_det,
    (s[5]*s[3]-s[1]*s[4])*0.5*inv_det,
    (s[3]*s[4]-s[2]*s[5])*0.5*inv_det,
  };
  double u[3] = {
    t[0]*s[0]+t[5]*s[1]+t[4]*s[2],
    t[5]*s[0]+t[1]*s[1]+t[3]*s[2],
    t[4]*s[0]+t[3]*s[1]+t[2]*s[2],
  };
  return  0.5*(u[0]*s[0]+u[1]*s[1]+u[2]*s[2]);
  /*
   const double square_radius = 0.5*(u[0]*s[0]+u[1]*s[1]+u[2]*s[2]);
   CVector3 vec1;
   vec1.p[0] = base[0][0]*u[0]+base[1][0]*u[1]+base[2][0]*u[2] + ipo0.p[0];
   vec1.p[1] = base[0][1]*u[0]+base[1][1]*u[1]+base[2][1]*u[2] + ipo0.p[1];
   vec1.p[2] = base[0][2]*u[0]+base[1][2]*u[1]+base[2][2]*u[2] + ipo0.p[2];
   std::cout << square_radius << " ";
   std::cout << SquareLength(vec1,ipo0) << " ";
   std::cout << SquareLength(vec1,ipo1) << " ";
   std::cout << SquareLength(vec1,ipo2) << " ";
   std::cout << SquareLength(vec1,ipo3) << std::endl;;
   return square_radius;
   */
}

template <typename T>
delfem2::CVec3<T> delfem2::CircumCenter(
    const CVec3<T>& ipo0,
    const CVec3<T>& ipo1,
    const CVec3<T>& ipo2,
    const CVec3<T>& ipo3)
{
  namespace lcl = delfem2::delaunay3;
  double base[3][3] = {
    { ipo1.p[0]-ipo0.p[0], ipo1.p[1]-ipo0.p[1], ipo1.p[2]-ipo0.p[2] },
    { ipo2.p[0]-ipo0.p[0], ipo2.p[1]-ipo0.p[1], ipo2.p[2]-ipo0.p[2] },
    { ipo3.p[0]-ipo0.p[0], ipo3.p[1]-ipo0.p[1], ipo3.p[2]-ipo0.p[2] }
  };
  double s[6] = {
    base[0][0]*base[0][0]+base[0][1]*base[0][1]+base[0][2]*base[0][2],
    base[1][0]*base[1][0]+base[1][1]*base[1][1]+base[1][2]*base[1][2],
    base[2][0]*base[2][0]+base[2][1]*base[2][1]+base[2][2]*base[2][2],
    base[1][0]*base[2][0]+base[1][1]*base[2][1]+base[1][2]*base[2][2],
    base[2][0]*base[0][0]+base[2][1]*base[0][1]+base[2][2]*base[0][2],
    base[0][0]*base[1][0]+base[0][1]*base[1][1]+base[0][2]*base[1][2],
  };
  const double vol = lcl::Volume_Tet3(ipo0.p,ipo1.p,ipo2.p,ipo3.p)*6.0;
  if( vol < 1.0e-20 ){ assert(0); }
  const double inv_det = 1.0 / (vol*vol);
  double t[6] = {
    (s[1]*s[2]-s[3]*s[3])*0.5*inv_det,
    (s[2]*s[0]-s[4]*s[4])*0.5*inv_det,
    (s[0]*s[1]-s[5]*s[5])*0.5*inv_det,
    (s[4]*s[5]-s[0]*s[3])*0.5*inv_det,
    (s[5]*s[3]-s[1]*s[4])*0.5*inv_det,
    (s[3]*s[4]-s[2]*s[5])*0.5*inv_det,
  };
  double u[3] = {
    t[0]*s[0]+t[5]*s[1]+t[4]*s[2],
    t[5]*s[0]+t[1]*s[1]+t[3]*s[2],
    t[4]*s[0]+t[3]*s[1]+t[2]*s[2],
  };
  //    const double square_radius = 0.5*(u[0]*s[0]+u[1]*s[1]+u[2]*s[2]);
  CVec3<T> vec1;
  vec1.p[0] = base[0][0]*u[0]+base[1][0]*u[1]+base[2][0]*u[2] + ipo0.p[0];
  vec1.p[1] = base[0][1]*u[0]+base[1][1]*u[1]+base[2][1]*u[2] + ipo0.p[1];
  vec1.p[2] = base[0][2]*u[0]+base[1][2]*u[1]+base[2][2]*u[2] + ipo0.p[2];
  return vec1;
}
#ifdef DFM2_STATIC_LIBRARY
template delfem2::CVec3d delfem2::CircumCenter(const CVec3d& ipo0,
                                         const CVec3d& ipo1,
                                         const CVec3d& ipo2,
                                         const CVec3d& ipo3);
#endif
  
  
// ------------------------------------------------

/**
 * @brief check if Delaunay condition satisfied
 * @details
 * 0 : p3 is inside circum circle on the p0,p1,p2
 * 1 :       on
 * 2 :       outsdie
 */
template <typename T>
int delfem2::DetDelaunay
(const CVec3<T>& p0,
const CVec3<T>& p1,
const CVec3<T>& p2,
const CVec3<T>& p3)
{
  const double area = Area_Tri(p0, p1, p2);
  if (fabs(area) < 1.0e-10){
    return 3;
  }
  const double tmp_val = 1.0/(area*area*16.0);

  const double dtmp0 = SquareDistance(p1, p2);
  const double dtmp1 = SquareDistance(p0, p2);
  const double dtmp2 = SquareDistance(p0, p1);

  const double etmp0 = tmp_val*dtmp0*(dtmp1+dtmp2-dtmp0);
  const double etmp1 = tmp_val*dtmp1*(dtmp0+dtmp2-dtmp1);
  const double etmp2 = tmp_val*dtmp2*(dtmp0+dtmp1-dtmp2);

  const CVec3<T> out_center = etmp0*p0+etmp1*p1+etmp2*p2;

  const double qradius = SquareDistance(out_center, p0);
  const double qdistance = SquareDistance(out_center, p3);

  const double tol = 1.0e-20;
  if (qdistance > qradius*(1.0+tol)){ return 2; }	// outside the circumcircle
  else{
    if (qdistance < qradius*(1.0-tol)){ return 0; }	// inside the circumcircle
    else{ return 1; }	// on the circumcircle
  }
  return 0;
}
#ifdef DFM2_STATIC_LIBRARY
template int delfem2::DetDelaunay(const CVec3d& p0,
                               const CVec3d& p1,
                               const CVec3d& p2,
                               const CVec3d& p3);
#endif
  
// -------------------------------------------

/**
 * @brief curcumradius of a tetrahedra
 */
template <typename T>
double delfem2::Circumradius
(const CVec3<T>& ipo0,
 const CVec3<T>& ipo1,
 const CVec3<T>& ipo2,
 const CVec3<T>& ipo3)
{
  return sqrt( SquareCircumradius(ipo0,ipo1,ipo2,ipo3) );
}

