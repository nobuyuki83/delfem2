#ifndef aabb_h
#define aabb_h

#include <math.h>
#include <assert.h>

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
	CBV3D_AABB( const CBV3D_AABB& bb )
  : x_min(bb.x_min),x_max(bb.x_max),
  y_min(bb.y_min),y_max(bb.y_max),
  z_min(bb.z_min),z_max(bb.z_max),
  is_active(bb.is_active){}
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
		x_max = ( x_max > bb.x_max ) ? x_max : bb.x_max;
		x_min = ( x_min < bb.x_min ) ? x_min : bb.x_min;
		y_max = ( y_max > bb.y_max ) ? y_max : bb.y_max;
		y_min = ( y_min < bb.y_min ) ? y_min : bb.y_min;
		z_max = ( z_max > bb.z_max ) ? z_max : bb.z_max;
		z_min = ( z_min < bb.z_min ) ? z_min : bb.z_min;
		return *this;
	}
  bool IsIntersect(const CBV3D_AABB& bb) const
  {
    if( !is_active ) return false;
    if( !bb.is_active ) return false;
    if( x_max < bb.x_min ) return false;
    if( x_min > bb.x_max ) return false;
    if( y_max < bb.y_min ) return false;
    if( y_min > bb.y_max ) return false;
    if( z_max < bb.z_min ) return false;
    if( z_min > bb.z_max ) return false;
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
		x_max = ( x_max > v[0] ) ? x_max : v[0];
		x_min = ( x_min < v[0] ) ? x_min : v[0];
		y_max = ( y_max > v[1] ) ? y_max : v[1];
		y_min = ( y_min < v[1] ) ? y_min : v[1];
		z_max = ( z_max > v[2] ) ? z_max : v[2];
		z_min = ( z_min < v[2] ) ? z_min : v[2];
		return *this;
	}
  void AddPoint(double x, double y, double z, double eps){
    if( eps <= 0 ){ return; }
    if( is_active ){ // something inside
      x_min = ( x_min < x-eps ) ? x_min : x-eps;
      y_min = ( y_min < y-eps ) ? y_min : y-eps;
      z_min = ( z_min < z-eps ) ? z_min : z-eps;
      x_max = ( x_max > x+eps ) ? x_max : x+eps;
      y_max = ( y_max > y+eps ) ? y_max : y+eps;
      z_max = ( z_max > z+eps ) ? z_max : z+eps;
    }
    else{ // empty
      is_active = true;
      x_min = x-eps;
      y_min = y-eps;
      z_min = z-eps;
      x_max = x+eps;
      y_max = y+eps;
      z_max = z+eps;
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
