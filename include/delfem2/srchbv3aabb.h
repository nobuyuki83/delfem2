/*
 * Copyright (c) 2019 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */


#ifndef DFM2_SRHBV3AABB_H
#define DFM2_SRHBV3AABB_H

#include <math.h>
#include <assert.h>
#include <vector>
#include <iostream>
#include <sstream>

namespace delfem2 {

/**
 * @class 3D axis alligned bounding box
 */
template <typename REAL>
class CBV3_AABB
{
public:
	CBV3_AABB(): bbmin{+1,0,0},bbmax{-1,0,0} {}
	CBV3_AABB(const CBV3_AABB<REAL>& bb ):
    bbmin{bb.bbmin[0], bb.bbmin[1], bb.bbmin[2]},
    bbmax{bb.bbmax[0], bb.bbmax[1], bb.bbmax[2]} {}
  CBV3_AABB(const std::vector<REAL>& bbmin, const std::vector<REAL>& bbmax):
    bbmin{bbmin[0], bbmin[1], bbmin[2]},
    bbmax{bbmax[0], bbmax[1], bbmax[2]} {}
  CBV3_AABB(const REAL bbmin[3], const REAL bbmax[3]):
    bbmin{bbmin[0], bbmin[1], bbmin[2]},
    bbmax{bbmax[0], bbmax[1], bbmax[2]} {}
  CBV3_AABB(const std::vector<REAL>& aabbvec3):
    bbmin{aabbvec3[0], aabbvec3[2], aabbvec3[4]},
    bbmax{aabbvec3[1], aabbvec3[3], aabbvec3[5]} {}
  // ------------
  bool IsActive() const {
    if( bbmin[0] > bbmax[0] ){ return false; }
    return true;
  }
  void Set_Inactive() {
	  bbmin[0] = +1;
	  bbmax[0] = -1;
	}
  REAL DiagonalLength() const{
    REAL x0 = bbmax[0] - bbmin[0];
    REAL y0 = bbmax[1] - bbmin[1];
    REAL z0 = bbmax[2] - bbmin[2];
    return sqrt(x0*x0+y0*y0+z0*z0);
  }
  REAL MaxLength() const{
    REAL x0 = bbmax[0] - bbmin[0];
    REAL y0 = bbmax[1] - bbmin[1];
    REAL z0 = bbmax[2] - bbmin[2];
    if( x0 > y0 && x0 > z0 ){ return x0; }
    else if( y0 > x0 && y0 > z0 ){ return y0; }
    return z0;
  }
  void SetCenterWidth(REAL cx, REAL cy, REAL cz,
                      REAL wx, REAL wy, REAL wz)
  {
    bbmin[0] = cx-wx*0.5; bbmax[0] = cx+wx*0.5;
    bbmin[1] = cy-wy*0.5; bbmax[1] = cy+wy*0.5;
    bbmin[2] = cz-wz*0.5; bbmax[2] = cz+wz*0.5;
  }
  void GetCenterWidth(REAL& cx, REAL& cy, REAL& cz,
                      REAL& wx, REAL& wy, REAL& wz)
  {
    cx = (bbmax[0]+bbmin[0])*0.5;
    cy = (bbmax[1]+bbmin[1])*0.5;
    cz = (bbmax[2]+bbmin[2])*0.5;
    wx = (bbmax[0]-bbmin[0]);
    wy = (bbmax[1]-bbmin[1]);
    wz = (bbmax[2]-bbmin[2]);
  }
	CBV3_AABB<REAL>& operator+=(const CBV3_AABB<REAL>& bb)
	{
		if( !bb.IsActive() ) return *this;
		if( !this->IsActive() ){
			bbmax[0] = bb.bbmax[0];	bbmin[0] = bb.bbmin[0];
			bbmax[1] = bb.bbmax[1];	bbmin[1] = bb.bbmin[1];
			bbmax[2] = bb.bbmax[2];	bbmin[2] = bb.bbmin[2];
			return *this;
		}
    bbmin[0] = ( bbmin[0] < bb.bbmin[0] ) ? bbmin[0] : bb.bbmin[0];
    bbmin[1] = ( bbmin[1] < bb.bbmin[1] ) ? bbmin[1] : bb.bbmin[1];
    bbmin[2] = ( bbmin[2] < bb.bbmin[2] ) ? bbmin[2] : bb.bbmin[2];
		bbmax[0] = ( bbmax[0] > bb.bbmax[0] ) ? bbmax[0] : bb.bbmax[0];
		bbmax[1] = ( bbmax[1] > bb.bbmax[1] ) ? bbmax[1] : bb.bbmax[1];
		bbmax[2] = ( bbmax[2] > bb.bbmax[2] ) ? bbmax[2] : bb.bbmax[2];
		return *this;
	}
  void Add_AABBVec3(const std::vector<REAL>& bbvec)
  {
    assert(bbvec.size()==6);
    CBV3_AABB aabb0(bbvec);
    (*this) += aabb0;
  }
  void Set_AABBVec3(const std::vector<REAL>& bbvec)
  {
    assert(bbvec.size()==6);
    this->bbmin[0] = bbvec[0];
    this->bbmin[1] = bbvec[2];
    this->bbmin[2] = bbvec[4];
    this->bbmax[0] = bbvec[1];
    this->bbmax[1] = bbvec[3];
    this->bbmax[2] = bbvec[5];
  }
  void Add_MinMax(const std::vector<REAL>& bbmin_,
                  const std::vector<REAL>& bbmax_)
  {
    assert(bbmin_.size()==3);
    assert(bbmax_.size()==3);
    CBV3_AABB aabb0(bbmin_.data(), bbmax_.data());
    (*this) += aabb0;
  }
  void Set_MinMax(const std::vector<REAL>& bbmin_,
                  const std::vector<REAL>& bbmax_)
  {
    assert(bbmin_.size()==3);
    assert(bbmax_.size()==3);
    this->bbmin[0] = bbmin_[0];
    this->bbmin[1] = bbmin_[1];
    this->bbmin[2] = bbmin_[2];
    this->bbmax[0] = bbmax_[0];
    this->bbmax[1] = bbmax_[1];
    this->bbmax[2] = bbmax_[2];
  }
  bool IsIntersect(const CBV3_AABB<REAL>& bb) const
  {
    if( !IsActive() ) return false;
    if( !bb.IsActive() ) return false;
    if( bbmin[0] > bb.bbmax[0] ) return false;
    if( bbmin[1] > bb.bbmax[1] ) return false;
    if( bbmin[2] > bb.bbmax[2] ) return false;
    //
    if( bbmax[0] < bb.bbmin[0] ) return false;
    if( bbmax[1] < bb.bbmin[1] ) return false;
    if( bbmax[2] < bb.bbmin[2] ) return false;
    return true;
  }
  CBV3_AABB<REAL>& operator+=(const REAL v[3])
	{
		if( !IsActive() ){
			bbmax[0] = bbmin[0] = v[0];
			bbmax[1] = bbmin[1] = v[1];
			bbmax[2] = bbmin[2] = v[2];
			return *this;
		}
    bbmin[0] = ( bbmin[0] < v[0] ) ? bbmin[0] : v[0];
    bbmin[1] = ( bbmin[1] < v[1] ) ? bbmin[1] : v[1];
    bbmin[2] = ( bbmin[2] < v[2] ) ? bbmin[2] : v[2];
    //
    bbmax[0] = ( bbmax[0] > v[0] ) ? bbmax[0] : v[0];
    bbmax[1] = ( bbmax[1] > v[1] ) ? bbmax[1] : v[1];
		bbmax[2] = ( bbmax[2] > v[2] ) ? bbmax[2] : v[2];
		return *this;
	}
	/**
	 * expand the aabb such that a point is included with a margin
	 * @param x x-coordinte
	 * @param y y-coorinate
	 * @param z z-coordinate
	 * @param eps margin. if negative, do nothing
	 */
  void AddPoint(const REAL p[3], REAL eps){
    if( eps < 0 ){ return; }
    if( !this->IsActive() ){ // something inside
      bbmin[0] = p[0]-eps;  bbmax[0] = p[0]+eps;
      bbmin[1] = p[1]-eps;  bbmax[1] = p[1]+eps;
      bbmin[2] = p[2]-eps;  bbmax[2] = p[2]+eps;
      return;
    }
    bbmin[0] = ( bbmin[0] < p[0]-eps ) ? bbmin[0] : p[0]-eps;
    bbmin[1] = ( bbmin[1] < p[1]-eps ) ? bbmin[1] : p[1]-eps;
    bbmin[2] = ( bbmin[2] < p[2]-eps ) ? bbmin[2] : p[2]-eps;
    bbmax[0] = ( bbmax[0] > p[0]+eps ) ? bbmax[0] : p[0]+eps;
    bbmax[1] = ( bbmax[1] > p[1]+eps ) ? bbmax[1] : p[1]+eps;
    bbmax[2] = ( bbmax[2] > p[2]+eps ) ? bbmax[2] : p[2]+eps;
  }
  REAL MinimumDistance(REAL x, REAL y, REAL z) const
  {
    REAL x0, y0, z0;
    if(      x < bbmin[0] ){ x0 = bbmin[0]; }
    else if( x < bbmax[0] ){ x0 = x;     }
    else{                    x0 = bbmax[0]; }
    if(      y < bbmin[1] ){ y0 = bbmin[1]; }
    else if( y < bbmax[1] ){ y0 = y;     }
    else{                    y0 = bbmax[1]; }
    if(      z < bbmin[2] ){ z0 = bbmin[2]; }
    else if( z < bbmax[2] ){ z0 = z;     }
    else{                    z0 = bbmax[2]; }
    return sqrt( (x0-x)*(x0-x) + (y0-y)*(y0-y) + (z0-z)*(z0-z) );
  }
  bool isInclude_Point(REAL x, REAL y, REAL z) const
  {
   if( !IsActive() ) return false;
   if(   x >= bbmin[0] && x <= bbmax[0]
      && y >= bbmin[1] && y <= bbmax[1]
      && z >= bbmin[2] && z <= bbmax[2] ) return true;
   return false;
  }
  bool IsInclude_AABB3(const CBV3_AABB<REAL>& bb) const
  {
    if( !this->IsActive() ) return false;
    if( !bb.IsActive() ) return true;
    if( this->bbmin[0] <= bb.bbmin[0] &&
        this->bbmin[1] <= bb.bbmin[1] &&
        this->bbmin[2] <= bb.bbmin[2] &&
        this->bbmax[0] >= bb.bbmax[0] &&
        this->bbmax[1] >= bb.bbmax[1] &&
        this->bbmax[2] >= bb.bbmax[2] ) return true;
    return false;
  }
  std::vector<REAL> AABBVec3(){
    return {
      bbmin[0], bbmax[0],
      bbmin[1], bbmax[1],
      bbmin[2], bbmax[2] };
  }
  std::vector<REAL> Center(){
    return {
      (bbmin[0]+bbmax[0])*0.5,
      (bbmin[1]+bbmax[1])*0.5,
      (bbmin[2]+bbmax[2])*0.5 };
  }
  std::vector<REAL> Point3D_Vox(){
    std::vector<REAL> aP(3*8);
    aP[0*3+0]=bbmin[0]; aP[0*3+1]=bbmin[1]; aP[0*3+2]=bbmin[2];
    aP[1*3+0]=bbmax[0]; aP[1*3+1]=bbmin[1]; aP[1*3+2]=bbmin[2];
    aP[2*3+0]=bbmin[0]; aP[2*3+1]=bbmax[1]; aP[2*3+2]=bbmin[2];
    aP[3*3+0]=bbmax[0]; aP[3*3+1]=bbmax[1]; aP[3*3+2]=bbmin[2];
    aP[4*3+0]=bbmin[0]; aP[4*3+1]=bbmin[1]; aP[4*3+2]=bbmax[2];
    aP[5*3+0]=bbmax[0]; aP[5*3+1]=bbmin[1]; aP[5*3+2]=bbmax[2];
    aP[6*3+0]=bbmin[0]; aP[6*3+1]=bbmax[1]; aP[6*3+2]=bbmax[2];
    aP[7*3+0]=bbmax[0]; aP[7*3+1]=bbmax[1]; aP[7*3+2]=bbmax[2];
    return aP;
  }
  std::string str(){
    std::stringstream ss;
    ss << bbmin[0] << " " << bbmax[0] << " " << bbmin[1] << " " << bbmax[1] << " " << bbmin[2] << " " << bbmax[2];
    return ss.str();
  }
public:
  // the order of this declarations should not be changed since it is used by cuda BVH.
  REAL bbmin[3], bbmax[3];
};
using CBV3d_AABB = CBV3_AABB<double>;
using CBV3f_AABB = CBV3_AABB<float>;
  
} // namespace delfem2

#endif
