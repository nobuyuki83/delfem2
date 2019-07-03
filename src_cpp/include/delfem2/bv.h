/*
 * Copyright (c) 2019 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */


#ifndef BV_H
#define BV_H

#include <math.h>
#include <assert.h>
#include <vector>
#include <iostream>

//! 3D bounding box class
class CBV3D_AABB
{
public:
	CBV3D_AABB(){
		x_min=0;	x_max=0;
		y_min=0;	y_max=0;
		z_min=0;	z_max=0;
		is_active = false;
	}
	CBV3D_AABB(double x_min0,double x_max0,
          double y_min0,double y_max0,
          double z_min0,double z_max0)
  : x_min(x_min0),x_max(x_max0),
  y_min(y_min0),y_max(y_max0),
  z_min(z_min0),z_max(z_max0)
	{
		assert( x_min <= x_max );
		assert( y_min <= y_max );
		assert( z_min <= z_max );
		is_active = true;
	}
	CBV3D_AABB(const CBV3D_AABB& bb )
  : x_min(bb.x_min),x_max(bb.x_max),
  y_min(bb.y_min),y_max(bb.y_max),
  z_min(bb.z_min),z_max(bb.z_max),
  is_active(bb.is_active){
  }
  CBV3D_AABB(const std::vector<double>& minmaxXYZ)
  {
    if( minmaxXYZ[0] > minmaxXYZ[1] ){
      x_min=0;  x_max=0;
      y_min=0;  y_max=0;
      z_min=0;  z_max=0;
      is_active = false;
      return;
    }
    x_min = minmaxXYZ[0];
    x_max = minmaxXYZ[1];
    y_min = minmaxXYZ[2];
    y_max = minmaxXYZ[3];
    z_min = minmaxXYZ[4];
    z_max = minmaxXYZ[5];
    is_active = true;
  }
  double DiagonalLength() const{
    double x0 = x_max - x_min;
    double y0 = y_max - y_min;
    double z0 = z_max - z_min;
    return sqrt(x0*x0+y0*y0+z0*z0);
  }
  double MaxLength() const{
    double x0 = x_max - x_min;
    double y0 = y_max - y_min;
    double z0 = z_max - z_min;
    if( x0 > y0 && x0 > z0 ){ return x0; }
    else if( y0 > x0 && y0 > z0 ){ return y0; }
    return z0;
  }
  void SetMinMaxXYZ(double x_min, double x_max,
                    double y_min, double y_max,
                    double z_min, double z_max)
  {
    this->x_min = x_min;  this->x_max = x_max;
    this->y_min = y_min;  this->y_max = y_max;
    this->z_min = z_min;  this->z_max = z_max;
  }
  void SetCenterWidth(double cx, double cy, double cz,
                      double wx, double wy, double wz)
  {
    x_min = cx-wx*0.5; x_max = cx+wx*0.5;
    y_min = cy-wy*0.5; y_max = cy+wy*0.5;
    z_min = cz-wz*0.5; z_max = cz+wz*0.5;
  }
  void GetCenterWidth(double& cx, double& cy, double& cz,
                      double& wx, double& wy, double& wz)
  {
    cx = (x_max+x_min)*0.5;
    cy = (y_max+y_min)*0.5;
    cz = (z_max+z_min)*0.5;
    wx = (x_max-x_min);
    wy = (y_max-y_min);
    wz = (z_max-z_min);
  }
	CBV3D_AABB& operator+=(const CBV3D_AABB& bb)
	{
		if( !bb.is_active ) return *this;
		if( !is_active ){
			x_max = bb.x_max;	x_min = bb.x_min;
			y_max = bb.y_max;	y_min = bb.y_min;
			z_max = bb.z_max;	z_min = bb.z_min;
      this->is_active = bb.is_active;
			return *this;
		}
    x_min = ( x_min < bb.x_min ) ? x_min : bb.x_min;
		x_max = ( x_max > bb.x_max ) ? x_max : bb.x_max;
    y_min = ( y_min < bb.y_min ) ? y_min : bb.y_min;
		y_max = ( y_max > bb.y_max ) ? y_max : bb.y_max;
    z_min = ( z_min < bb.z_min ) ? z_min : bb.z_min;
		z_max = ( z_max > bb.z_max ) ? z_max : bb.z_max;
		return *this;
	}
  void Add_AABBMinMax(const std::vector<double>& aabb){
    assert(aabb.size()==6);
    CBV3D_AABB aabb0(aabb);
    (*this) += aabb0;
  }
  bool IsIntersect(const CBV3D_AABB& bb) const
  {
    if( !is_active ) return false;
    if( !bb.is_active ) return false;
    if( x_min > bb.x_max ) return false;
    if( x_max < bb.x_min ) return false;
    if( y_min > bb.y_max ) return false;
    if( y_max < bb.y_min ) return false;
    if( z_min > bb.z_max ) return false;
    if( z_max < bb.z_min ) return false;
    return true;
  }
  CBV3D_AABB& operator+=(const double v[3])
	{
		if( !is_active ){
			x_max = v[0];	x_min = v[0];
			y_max = v[1];	y_min = v[1];
			z_max = v[2];	z_min = v[2];
      this->is_active = true;
			return *this;
		}
    x_min = ( x_min < v[0] ) ? x_min : v[0];
		x_max = ( x_max > v[0] ) ? x_max : v[0];
    y_min = ( y_min < v[1] ) ? y_min : v[1];
		y_max = ( y_max > v[1] ) ? y_max : v[1];
    z_min = ( z_min < v[2] ) ? z_min : v[2];
		z_max = ( z_max > v[2] ) ? z_max : v[2];
		return *this;
	}
  void AddPoint(double x, double y, double z, double eps){
    if( eps <= 0 ){ return; }
    if( is_active ){ // something inside
      x_min = ( x_min < x-eps ) ? x_min : x-eps;
      x_max = ( x_max > x+eps ) ? x_max : x+eps;
      y_min = ( y_min < y-eps ) ? y_min : y-eps;
      y_max = ( y_max > y+eps ) ? y_max : y+eps;
      z_min = ( z_min < z-eps ) ? z_min : z-eps;
      z_max = ( z_max > z+eps ) ? z_max : z+eps;
    }
    else{ // empty
      is_active = true;
      x_min = x-eps;  x_max = x+eps;
      y_min = y-eps;  y_max = y+eps;
      z_min = z-eps;  z_max = z+eps;
    }
    return;
  }
  double MinimumDistance(double x, double y, double z) const
  {
    double x0, y0, z0;
    if(      x < x_min ){ x0 = x_min; }
    else if( x < x_max ){ x0 = x;     }
    else{                 x0 = x_max; }
    if(      y < y_min ){ y0 = y_min; }
    else if( y < y_max ){ y0 = y;     }
    else{                 y0 = y_max; }
    if(      z < z_min ){ z0 = z_min; }
    else if( z < z_max ){ z0 = z;     }
    else{                 z0 = z_max; }
    return sqrt( (x0-x)*(x0-x) + (y0-y)*(y0-y) + (z0-z)*(z0-z) );
  }
  bool IsInside(double x, double y, double z) const
  {
   if( !is_active ) return false;
   if(   x >= x_min && x <= x_max
      && y >= y_min && y <= y_max
      && z >= z_min && z <= z_max ) return true;
   return false;
  }
  std::vector<double> MinMaxXYZ(){
    std::vector<double> mm(6);
    if( !this->is_active ){
      mm[0] = +1;
      mm[1] = -1;
      return mm;
    }
    mm[0] = x_min;  mm[1] = x_max;
    mm[2] = y_min;  mm[3] = y_max;
    mm[4] = z_min;  mm[5] = z_max;
    return mm;
  }
  std::vector<double> Center(){
    std::vector<double> mm(3);
    mm[0] = (x_min+x_max)*0.5;
    mm[1] = (y_min+y_max)*0.5;
    mm[2] = (z_min+z_max)*0.5;
    return mm;
  }
  std::vector<double> Point3D_Vox(){
    std::vector<double> aP(3*8);
    aP[0*3+0]=x_min; aP[0*3+1]=y_min; aP[0*3+2]=z_min;
    aP[1*3+0]=x_max; aP[1*3+1]=y_min; aP[1*3+2]=z_min;
    aP[2*3+0]=x_min; aP[2*3+1]=y_max; aP[2*3+2]=z_min;
    aP[3*3+0]=x_max; aP[3*3+1]=y_max; aP[3*3+2]=z_min;
    aP[4*3+0]=x_min; aP[4*3+1]=y_min; aP[4*3+2]=z_max;
    aP[5*3+0]=x_max; aP[5*3+1]=y_min; aP[5*3+2]=z_max;
    aP[6*3+0]=x_min; aP[6*3+1]=y_max; aP[6*3+2]=z_max;
    aP[7*3+0]=x_max; aP[7*3+1]=y_max; aP[7*3+2]=z_max;
    return aP;
  }
  std::string str(){
    return std::string(std::to_string(x_min)+" "+std::to_string(x_max)+" "+
                       std::to_string(y_min)+" "+std::to_string(y_max)+" "+
                       std::to_string(z_min)+" "+std::to_string(z_max));
  }
public:
	double x_min,x_max,  y_min,y_max,  z_min,z_max;
	bool is_active;	//!< false if there is nothing inside
};

//////////////////////////////////////////////////////////////////////

class CBV3D_Sphere
{
public:
  bool is_active;
  double cx,cy,cz,r;
public:
  CBV3D_Sphere(){
    is_active = false;
    cx=cy=cz=r=0;
  }
  void AddPoint(double x,double y,double z, double R){
    if( R <= 0 ){ return; }
    if( !is_active ){ // empty
      is_active = true;
      cx=x; cy=y; cz=z; r=R;
      return;
    }
    const double L = sqrt((x-cx)*(x-cx) + (y-cy)*(y-cy) + (z-cz)*(z-cz));
    if( L<r+R ){ // overlapping
      if( r>L+R ){ return; } // included
      if( R>L+r){
        cx=x; cy=y; cz=z; r=R;
        return;
      }
    }
    if( fabs(L) < 1.0e-5*fabs(r+R) ){
      r = L+R;
      return;
    }
    const double r0 = 0.5*(1+(r-R)/L);
    const double r1 = 0.5*(1+(R-r)/L);
    cx = r0*cx + r1*x;
    cy = r0*cy + r1*y;
    cz = r0*cz + r1*z;
    r = 0.5*(L+r+R);
    return;
  }
  bool IsIntersect(const CBV3D_Sphere& bb) const
  {
    if( !is_active ) return false;
    if( !bb.is_active ) return false;
    const double L = sqrt((bb.cx-cx)*(bb.cx-cx) + (bb.cy-cy)*(bb.cy-cy) + (bb.cz-cz)*(bb.cz-cz));
    if( L > bb.r + r ) return false;
    return true;
  }
  CBV3D_Sphere& operator+=(const CBV3D_Sphere& bb)
  {
    if( !bb.is_active ) return *this;
    this->AddPoint(bb.cx,bb.cy,bb.cz, bb.r);
    return *this;
  }
  bool isInclude_Point(double x, double y, double z) const {
    if( !is_active ){ return false; }
    const double L = (x-cx)*(x-cx) + (y-cy)*(y-cy) + (z-cz)*(z-cz);
    if( L < r*r ){ return true; }
    return false;
  }
  void getRange_Point(double& min0, double& max0,
                      double x, double y, double z) const {
    const double L = sqrt((x-cx)*(x-cx) + (y-cy)*(y-cy) + (z-cz)*(z-cz));
    if( L < r ){
      min0 = 0;
      max0 = r-L;
      return;
    }
    min0 = L-r;
    max0 = L+r;
  }
};

#endif
