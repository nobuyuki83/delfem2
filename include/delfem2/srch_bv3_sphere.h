/*
 * Copyright (c) 2019 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */


#ifndef DFM2_SRHBV3SPHERE_H
#define DFM2_SRHBV3SPHERE_H

#include <math.h>
#include <assert.h>
#include <vector>
#include <iostream>
#include <sstream>

namespace delfem2 {

/**
 * @class bounding volume with sphere
 */
template <typename REAL>
class CBV3_Sphere
{
public:
  CBV3_Sphere(){
    c[0]=c[1]=c[2]=0;
    r = -1; // if r is negative this is not active yet
  }
  void AddPoint(const REAL p[3], REAL R){
    if( R < 0 ){ return; }
    if( r < 0 ){ // empty
      c[0]=p[0]; c[1]=p[1]; c[2]=p[2]; r=R;
      return;
    }
    const REAL L = this->Distance3(p,c);
    if( r>L+R ){ return; } // including
    if( R>L+r){ // included
      c[0]=p[0]; c[1]=p[1]; c[2]=p[2]; r=R;
      return;
    }
    if( fabs(L) <= 1.0e-5*fabs(r+R) ){ // almost co-centric
      r = L+R;
      return;
    }
    const REAL r0 = 0.5*(L+r-R)/L;
    const REAL r1 = 0.5*(L+R-r)/L;
    assert( r0 >= 0 && r1 >= 0 );
    c[0] = r0*c[0] + r1*p[0];
    c[1] = r0*c[1] + r1*p[1];
    c[2] = r0*c[2] + r1*p[2];
    r = 0.5*(L+r+R);
    return;
  }
  bool IsIntersect(const CBV3_Sphere& bb) const
  {
    const double L = this->Distance3(bb.c,c);
    if( L > bb.r + r ) return false;
    return true;
  }
  void Set_Inactive() {
    r = -1.0;
  }
  bool IsActive() const {
    return r >= 0;
  }
  template <typename REAL1>
  bool IsIntersectLine(const REAL1 src[3], const REAL1 dir[3]) const {
    REAL ratio = dir[0]*(c[0]-src[0]) + dir[1]*(c[1]-src[1]) + dir[2]*(c[2]-src[2]);
    ratio = ratio/(dir[0]*dir[0] + dir[1]*dir[1] + dir[2]*dir[2]);
    const REAL px = src[0] + ratio*dir[0];
    const REAL py = src[1] + ratio*dir[1];
    const REAL pz = src[2] + ratio*dir[2];
    const REAL L = sqrt((px-c[0])*(px-c[0]) + (py-c[1])*(py-c[1]) + (pz-c[2])*(pz-c[2]));
    assert( fabs(dir[0]*(px-c[0]) + dir[1]*(py-c[1]) + dir[2]*(pz-c[2])) < 1.0e-4 );
    if( L <= r ){ return true; }
    return false;
  }
  bool IsIntersectRay(const double src[3], const double dir[3]) const {
    const double L0 = sqrt((src[0]-c[0])*(src[0]-c[0]) + (src[1]-c[1])*(src[1]-c[1]) + (src[2]-c[2])*(src[2]-c[2]));
    if( L0 <= r ){ return true; } // source included
    double ratio = dir[0]*(c[0]-src[0]) + dir[1]*(c[1]-src[1]) + dir[2]*(c[2]-src[2]);
    ratio = ratio/(dir[0]*dir[0] + dir[1]*dir[1] + dir[2]*dir[2]);
    if( ratio < 0 ){ return false; }
    const double px = src[0] + ratio*dir[0];
    const double py = src[1] + ratio*dir[1];
    const double pz = src[2] + ratio*dir[2];
    const double L = sqrt((px-c[0])*(px-c[0]) + (py-c[1])*(py-c[1]) + (pz-c[2])*(pz-c[2]));
    assert( fabs(dir[0]*(px-c[0]) + dir[1]*(py-c[1]) + dir[2]*(pz-c[2])) < 1.0e-10 );
    if( L <= r ){ return true; }
    return false;
  }
  CBV3_Sphere& operator+=(const CBV3_Sphere& bb)
  {
    this->AddPoint(bb.c, bb.r);
    return *this;
  }
  bool isInclude_Point(double x, double y, double z) const {
    const double L = (x-c[0])*(x-c[0]) + (y-c[1])*(y-c[1]) + (z-c[2])*(z-c[2]);
    if( L < r*r ){ return true; }
    return false;
  }
  bool IsInclude(const CBV3_Sphere& bv, double margin) const {
    const double LL = (bv.c[0]-c[0])*(bv.c[0]-c[0]) + (bv.c[1]-c[1])*(bv.c[1]-c[1]) + (bv.c[2]-c[2])*(bv.c[2]-c[2]);
    const double L = sqrt(LL);
    if( L+bv.r <= r+margin ){ return true; }
    return false;
  }
  /**
   * @brief minimum and maximum distance of this bounding box from a point (x,y,z)
   * do nothing when this bounding box is inactive
   */
  void Range_DistToPoint(REAL& min0, REAL& max0,
                         REAL x, REAL y, REAL z) const {
    if( r < 0 ){ return; }
    const REAL L = sqrt((x-c[0])*(x-c[0]) + (y-c[1])*(y-c[1]) + (z-c[2])*(z-c[2]));
    if( L < r ){
      min0 = 0;
      max0 = r+L;
      return;
    }
    min0 = L-r;
    max0 = L+r;
  }
  void SetPoints4(
      const REAL p0[3],
      const REAL p1[3],
      const REAL p2[3],
      const REAL p3[4],
      REAL cc) // clearance
  {
    assert(cc>=0);
    // the center of the gravity
    c[0] = (p0[0]+p1[0]+p2[0]+p3[0])*0.25;
    c[1] = (p0[1]+p1[1]+p2[1]+p3[1])*0.25;
    c[2] = (p0[2]+p1[2]+p2[2]+p3[2])*0.25;
    // distance to input points
    const REAL r0 = CBV3_Sphere<REAL>::Distance3(c,p0);
    const REAL r1 = CBV3_Sphere<REAL>::Distance3(c,p1);
    const REAL r2 = CBV3_Sphere<REAL>::Distance3(c,p2);
    const REAL r3 = CBV3_Sphere<REAL>::Distance3(c,p3);
    // pick the maximum distance
    r = (r1>r0)?r1:r0;
    r = (r2>r)?r2:r;
    r = (r3>r)?r3:r;
    r += cc;
  }
public:
  // the order of this declarations should not be changed since it is used by cuda BVH.
  REAL r;
  REAL c[3];
private:
  static REAL Distance3(const REAL p[3], const REAL q[3]) {
    return sqrt((p[0]-q[0])*(p[0]-q[0]) + (p[1]-q[1])*(p[1]-q[1]) + (p[2]-q[2])*(p[2]-q[2]) );
  }
  /*
  static float Distance3(const float p[3], const float q[3]) {
    return sqrtf((p[0]-q[0])*(p[0]-q[0]) + (p[1]-q[1])*(p[1]-q[1]) + (p[2]-q[2])*(p[2]-q[2]) );
  }
   */
};

//! @brief 3D bounding volume of sphere with "float" precision
using CBV3f_Sphere = CBV3_Sphere<float>;
//! @brief 3D bounding volume of sphere with "double" precision
using CBV3d_Sphere = CBV3_Sphere<double>;

  
} // namespace delfem2

#endif
