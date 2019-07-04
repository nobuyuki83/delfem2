/*
 * Copyright (c) 2019 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */


#include <vector>
#include <math.h>
#include <iostream>
#include <assert.h>

#include <fstream>
#include <time.h>

#include "delfem2/sdf.h"
#include "delfem2/hash3d.h"

static void Cross3D(double r[3], const double v1[3], const double v2[3]){
	r[0] = v1[1]*v2[2] - v2[1]*v1[2];
	r[1] = v1[2]*v2[0] - v2[2]*v1[0];
	r[2] = v1[0]*v2[1] - v2[0]*v1[1];
}


static inline double Dot3D(const double p0[3], const double p1[3]){
	return p0[0]*p1[0]+p0[1]*p1[1]+p0[2]*p1[2];
}


static inline double Length3D(const double v[3]){
	return sqrt( v[0]*v[0] + v[1]*v[1] + v[2]*v[2] );
}

static inline double SqLength3D(const double v[3]){
	return v[0]*v[0] + v[1]*v[1] + v[2]*v[2];
}

static double Distance3D(const double p0[3], const double p1[3]){
  return sqrt((p1[0]-p0[0])*(p1[0]-p0[0])+(p1[1]-p0[1])*(p1[1]-p0[1])+(p1[2]-p0[2])*(p1[2]-p0[2]));
}

/*
static inline void NormalTri3D(double n[3], const double v1[3], const double v2[3], const double v3[3]){
	n[0] = ( v2[1] - v1[1] )*( v3[2] - v1[2] ) - ( v3[1] - v1[1] )*( v2[2] - v1[2] );
	n[1] = ( v2[2] - v1[2] )*( v3[0] - v1[0] ) - ( v3[2] - v1[2] )*( v2[0] - v1[0] );
	n[2] = ( v2[0] - v1[0] )*( v3[1] - v1[1] ) - ( v3[0] - v1[0] )*( v2[1] - v1[1] );
}
 */

static double TetVolume3D
(const double v1[3],
 const double v2[3],
 const double v3[3],
 const double v4[3] )
{
	return	
	(   ( v2[0] - v1[0] )*( ( v3[1] - v1[1] )*( v4[2] - v1[2] ) - ( v4[1] - v1[1] )*( v3[2] - v1[2] ) )		
	 -	( v2[1] - v1[1] )*( ( v3[0] - v1[0] )*( v4[2] - v1[2] ) - ( v4[0] - v1[0] )*( v3[2] - v1[2] ) )		
	 +	( v2[2] - v1[2] )*( ( v3[0] - v1[0] )*( v4[1] - v1[1] ) - ( v4[0] - v1[0] )*( v3[1] - v1[1] ) )
	 ) * 0.16666666666666666666666666666667;		
}

/*
static double TriArea3D(const double v1[3], const double v2[3], const double v3[3]){
	double x, y, z;
	x = ( v2[1] - v1[1] )*( v3[2] - v1[2] ) - ( v3[1] - v1[1] )*( v2[2] - v1[2] );
	y = ( v2[2] - v1[2] )*( v3[0] - v1[0] ) - ( v3[2] - v1[2] )*( v2[0] - v1[0] );
	z = ( v2[0] - v1[0] )*( v3[1] - v1[1] ) - ( v3[0] - v1[0] )*( v2[1] - v1[1] );
	return 0.5*sqrt( x*x + y*y + z*z );
}
 */
static void  UnitNormalAreaTri3D(double n[3], double& a, const double v1[3], const double v2[3], const double v3[3]){
	n[0] = ( v2[1] - v1[1] )*( v3[2] - v1[2] ) - ( v3[1] - v1[1] )*( v2[2] - v1[2] );
	n[1] = ( v2[2] - v1[2] )*( v3[0] - v1[0] ) - ( v3[2] - v1[2] )*( v2[0] - v1[0] );
	n[2] = ( v2[0] - v1[0] )*( v3[1] - v1[1] ) - ( v3[0] - v1[0] )*( v2[1] - v1[1] );
	a = sqrt(n[0]*n[0]+n[1]*n[1]+n[2]*n[2])*0.5;
	const double invlen = 0.5/a;
	n[0]*=invlen;	n[1]*=invlen;	n[2]*=invlen;
}

// (v1*v2).v3
inline double ScalarTripleProduct3D(const double v1[3], const double v2[3], const double v3[3]){
	const double x = ( v1[1]*v2[2] - v2[1]*v1[2] )*v3[0];
	const double y = ( v1[2]*v2[0] - v2[2]*v1[0] )*v3[1];
	const double z = ( v1[0]*v2[1] - v2[0]*v1[1] )*v3[2];
	return x+y+z;
}

static bool GetVertical2Vector3D(const double vec0[3], double vec1[3], double vec2[3])
{
  double vec_a[3] = { vec0[0], vec0[1], vec0[2] };
  const double len = Length3D(vec0);
  if( len < 1.0e-20 ){ 
    assert(0);
    return false; 
  }
  vec_a[0] *= 1.0/len;
  vec_a[1] *= 1.0/len;
	vec_a[2] *= 1.0/len;
  double vec_t[3] = {0,1,0};
  Cross3D(vec1,vec_t,vec_a);
  if( Length3D(vec1) < 1.0e-10 ){
    double vec_t[3] = {1,0,0};
    Cross3D(vec1,vec_t,vec_a);  // z????
    Cross3D(vec2,vec_a,vec1);  // x????
  }
  else{
	  const double len = Length3D(vec1);
    vec1[0] *= 1.0/len;
    vec1[1] *= 1.0/len;
		vec1[2] *= 1.0/len;
    Cross3D(vec2,vec_a,vec1);
  }
  return true;
}

static void GetNearest_LineSegPoint3D
(double pn[3],
const double p[3], // point
const double s[3], // source
const double e[3]) // end
{
  const double d[3] = { e[0]-s[0], e[1]-s[1], e[2]-s[2] };
  double t = 0.5;
  if (Dot3D(d, d)>1.0e-20){
    const double ps[3] = { s[0]-p[0], s[1]-p[1], s[2]-p[2] };
    double a = Dot3D(d, d);
    double b = Dot3D(d, ps);
    t = -b/a;
    if (t<0) t = 0;
    if (t>1) t = 1;
  }
  pn[0] = s[0]+t*d[0];
  pn[1] = s[1]+t*d[1];
  pn[2] = s[2]+t*d[2];
  return;
}

static void GetNearest_TrianglePoint3D
(double pn[3],
 double& r0, double& r1,
const double ps[3], // origin point
const double q0[3],
const double q1[3],
const double q2[3])
{
  double area, n012[3]; UnitNormalAreaTri3D(n012, area, q0, q1, q2);
  const double pe[3] = { ps[0]+n012[0], ps[1]+n012[1], ps[2]+n012[2] };
  const double v012 = TetVolume3D(ps, q0, q1, q2);
  if (fabs(v012) > 1.0e-10){
    const double sign = (v012 > 0) ? +1 : -1;
    const double v0 = TetVolume3D(ps, q1, q2, pe)*sign;
    const double v1 = TetVolume3D(ps, q2, q0, pe)*sign;
    const double v2 = TetVolume3D(ps, q0, q1, pe)*sign;
    assert(fabs(v0+v1+v2) > 1.0e-10);
    double inv_v012 = 1.0/(v0+v1+v2);
    r0 = v0*inv_v012;
    r1 = v1*inv_v012;
    const double r2 = (1.0-r0-r1);
    const double tol = 1.0e-4;
    if (r0 > -tol && r1 > -tol && r2 > -tol){
      pn[0] = q0[0]*r0+q1[0]*r1+q2[0]*r2;
      pn[1] = q0[1]*r0+q1[1]*r1+q2[1]*r2;
      pn[2] = q0[2]*r0+q1[2]*r1+q2[2]*r2;
      return;
    }
  }
  double r12[3]; GetNearest_LineSegPoint3D(r12, ps, q1, q2);
  double r20[3]; GetNearest_LineSegPoint3D(r20, ps, q2, q0);
  double r01[3]; GetNearest_LineSegPoint3D(r01, ps, q0, q1);
  const double d12 = Distance3D(r12, ps);
  const double d20 = Distance3D(r20, ps);
  const double d01 = Distance3D(r01, ps);
  if (d12 < d20){
    if (d12 < d01){ // 12 is the smallest
      pn[0] = r12[0];
      pn[1] = r12[1];
      pn[2] = r12[2];
      r0 = 0;
      r1 = Distance3D(pn,q2)/Distance3D(q1,q2);
      return;
    }
  }
  else{
    if (d20 < d01){ // d20 is the smallest
      pn[0] = r20[0];
      pn[1] = r20[1];
      pn[2] = r20[2];
      r0 = Distance3D(pn,q2)/Distance3D(q0,q2);
      r1 = 0;
      return;
    }
  }
  pn[0] = r01[0];
  pn[1] = r01[1];
  pn[2] = r01[2];
  r0 = Distance3D(pn,q1)/Distance3D(q0,q1);
  r1 = 1-r0;
  return;
}

////////////////////////////////////////////////////////////////////////////////////////////

CSignedDistanceField3D_Plane::CSignedDistanceField3D_Plane(double n[3], double o[3])
{
	////
	normal_[0] = n[0];
	normal_[1] = n[1];
	normal_[2] = n[2];
	////
	origin_[0] = o[0];
  origin_[1] = o[1];
	origin_[2] = o[2];
}



double CSignedDistanceField3D_Plane::Projection(double px, double py, double pz,
							  double n[3]) const // normal
{
	n[0] = normal_[0];		
	n[1] = normal_[1];		
	n[2] = normal_[2];		
	return -( normal_[0]*(px-origin_[0]) + normal_[1]*(py-origin_[1]) + normal_[2]*(pz-origin_[2]) );
}


////////////////////////////////////////////////////////////////


CSignedDistanceField3D_Sphere::CSignedDistanceField3D_Sphere(double r, const std::vector<double>& c, bool is_out){
  cent_.resize(3);
	cent_[0] = c[0];
	cent_[1] = c[1];
	cent_[2] = c[2];
	radius_ = r;
	this->is_out_ = is_out;
}

void CSignedDistanceField3D_Sphere::GetMesh
(std::vector<unsigned int>& aTri, std::vector<double>& aXYZ, double elen) const
{
  double pi = 3.1415;
  const unsigned int nlg = 32;
	const unsigned int nlt = 18;	
	const double rlg = 2*pi/nlg;	// longtitude
	const double rlt = 1*pi/nlt;	// latitude
  aXYZ.resize(nlg*(nlt-1)*3+2*3);
	for(unsigned int ilg=0;ilg<nlg;ilg++){    
    for(unsigned int ilt=0;ilt<(nlt-1);ilt++){
      double alt = ilt*rlt+rlt-pi*0.5;
      aXYZ[(ilg*(nlt-1)+ilt)*3+0] = cent_[0]+radius_*cos(alt)*sin(ilg*rlg);
      aXYZ[(ilg*(nlt-1)+ilt)*3+1] = cent_[1]+radius_*cos(alt)*cos(ilg*rlg);
      aXYZ[(ilg*(nlt-1)+ilt)*3+2] = cent_[2]+radius_*sin(alt);
    }
  }
  const unsigned int iu1 = nlg*(nlt-1);
  {
    aXYZ[iu1*3+0] = cent_[0];
    aXYZ[iu1*3+1] = cent_[1];
    aXYZ[iu1*3+2] = cent_[2]+radius_;
  }
  const unsigned int iu2 = nlg*(nlt-1)+1;  
  {    
    aXYZ[iu2*3+0] = cent_[0];
    aXYZ[iu2*3+1] = cent_[1];
    aXYZ[iu2*3+2] = cent_[2]-radius_;   
  }
  ////
  aTri.resize(nlg*(nlt-2)*2*3+nlg*2*3);
	for(unsigned int ilg=0;ilg<nlg;ilg++){        
    for(unsigned int ilt=0;ilt<nlt-2;ilt++){
      unsigned int iug = ( ilg == nlg-1 ) ? 0 : ilg+1;
      aTri[(ilg*(nlt-2)+ilt)*6+0] = ilg*(nlt-1)+ilt;
      aTri[(ilg*(nlt-2)+ilt)*6+1] = iug*(nlt-1)+ilt+1;
      aTri[(ilg*(nlt-2)+ilt)*6+2] = iug*(nlt-1)+ilt;
      ////
      aTri[(ilg*(nlt-2)+ilt)*6+3] = ilg*(nlt-1)+ilt;
      aTri[(ilg*(nlt-2)+ilt)*6+4] = ilg*(nlt-1)+ilt+1;
      aTri[(ilg*(nlt-2)+ilt)*6+5] = iug*(nlt-1)+ilt+1;
    }
  }
  const unsigned int itri1 = nlg*(nlt-2)*2;
	for(unsigned int ilg=0;ilg<nlg;ilg++){        
    unsigned int iug = ( ilg == nlg-1 ) ? 0 : ilg+1;    
    aTri[(itri1+ilg)*3+0] = nlt-2+(nlt-1)*ilg;
    aTri[(itri1+ilg)*3+1] = iu1;
    aTri[(itri1+ilg)*3+2] = nlt-2+(nlt-1)*iug;
  }
  const unsigned int itri2 = nlg*(nlt-2)*2+nlg;
	for(unsigned int ilg=0;ilg<nlg;ilg++){        
    unsigned int iug = ( ilg == nlg-1 ) ? 0 : ilg+1;    
    aTri[(itri2+ilg)*3+0] = (nlt-1)*iug;
    aTri[(itri2+ilg)*3+1] = iu2;
    aTri[(itri2+ilg)*3+2] = (nlt-1)*ilg;    
  }    
}

	
// return penetration depth (inside is positive)
double CSignedDistanceField3D_Sphere::Projection
	(double px, double py, double pz,
	 double n[3]) const // normal outward
{
	double dir[3] = { px-cent_[0], py-cent_[1], pz-cent_[2] };
	const double len = sqrt( dir[0]*dir[0]+dir[1]*dir[1]+dir[2]*dir[2] );		
	const double invlen = 1.0/len;
	if( !is_out_ ){
		n[0] = -dir[0]*invlen;
		n[1] = -dir[1]*invlen;		
		n[2] = -dir[2]*invlen;
		return +len-radius_;
	}
	n[0] = dir[0]*invlen;
	n[1] = dir[1]*invlen;		
	n[2] = dir[2]*invlen;	
	return radius_-len;
}

unsigned int CSignedDistanceField3D_Sphere::FindInOut(double px, double py, double pz) const
{
	double n[3];
	double pd = this->Projection(px, py, pz, n);	
	if( !is_out_ ) pd *= -1.0;
	if( pd > 0 ){ return 0; }
	return 1;
}

bool CSignedDistanceField3D_Sphere::IntersectionPoint
(double p[3], 
 const double o[3], const double d[3]) const 
{
  const double q[3] = { o[0]-cent_[0], o[1]-cent_[1], o[2]-cent_[2] };
  const double a = d[0]*d[0] + d[1]*d[1] + d[2]*d[2];
  const double b = q[0]*d[0] + q[1]*d[1] + q[2]*d[2];
  const double c = q[0]*q[0] + q[1]*q[1] + q[2]*q[2] - radius_*radius_;
  const double det = b*b-a*c;
  if( det < 0 )  return false;
  const double t = (-b+sqrt(det))/a;
  p[0] = o[0] + t*d[0];
  p[1] = o[1] + t*d[1];
  p[2] = o[2] + t*d[2];  
  return true;
}


////////////////////////////////////////////////////////////////

CSignedDistanceField3D_Cylinder::CSignedDistanceField3D_Cylinder
(double r, double cnt[3], double dir[3], bool is_out){
	cent_[0] = cnt[0];
	cent_[1] = cnt[1];
	cent_[2] = cnt[2];
  dir_[0] = dir[0];
  dir_[1] = dir[1];
  dir_[2] = dir[2];
	radius_ = r;
	this->is_out_ = is_out;
}

void CSignedDistanceField3D_Cylinder::GetMesh
(std::vector<unsigned int>& aTri,
 std::vector<double>& aXYZ,
 double elen) const
{
  double pi = 3.1415;
  const unsigned int nlg = 32;
	const unsigned int nlt = 18;	
	const double rlg = 2*pi/nlg;	// longtitude
	const double rlt = 1*pi/nlt;	// latitude
  aXYZ.resize(nlg*(nlt-1)*3+2*3);
	for(unsigned int ilg=0;ilg<nlg;ilg++){    
    for(unsigned int ilt=0;ilt<(nlt-1);ilt++){
      double alt = ilt*rlt+rlt-pi*0.5;
      aXYZ[(ilg*(nlt-1)+ilt)*3+0] = cent_[0]+radius_*cos(alt)*sin(ilg*rlg);
      aXYZ[(ilg*(nlt-1)+ilt)*3+1] = cent_[1]+radius_*cos(alt)*cos(ilg*rlg);
      aXYZ[(ilg*(nlt-1)+ilt)*3+2] = cent_[2]+radius_*sin(alt);
    }
  }
  const unsigned int iu1 = nlg*(nlt-1);
  {
    aXYZ[iu1*3+0] = cent_[0];
    aXYZ[iu1*3+1] = cent_[1];
    aXYZ[iu1*3+2] = cent_[2]+radius_;
  }
  const unsigned int iu2 = nlg*(nlt-1)+1;  
  {    
    aXYZ[iu2*3+0] = cent_[0];
    aXYZ[iu2*3+1] = cent_[1];
    aXYZ[iu2*3+2] = cent_[2]-radius_;   
  }
  ////
  aTri.resize(nlg*(nlt-2)*2*3+nlg*2*3);
	for(unsigned int ilg=0;ilg<nlg;ilg++){        
    for(unsigned int ilt=0;ilt<nlt-2;ilt++){
      unsigned int iug = ( ilg == nlg-1 ) ? 0 : ilg+1;
      aTri[(ilg*(nlt-2)+ilt)*6+0] = ilg*(nlt-1)+ilt;
      aTri[(ilg*(nlt-2)+ilt)*6+1] = iug*(nlt-1)+ilt+1;
      aTri[(ilg*(nlt-2)+ilt)*6+2] = iug*(nlt-1)+ilt;
      ////
      aTri[(ilg*(nlt-2)+ilt)*6+3] = ilg*(nlt-1)+ilt;
      aTri[(ilg*(nlt-2)+ilt)*6+4] = ilg*(nlt-1)+ilt+1;
      aTri[(ilg*(nlt-2)+ilt)*6+5] = iug*(nlt-1)+ilt+1;
    }
  }
  const unsigned int itri1 = nlg*(nlt-2)*2;
	for(unsigned int ilg=0;ilg<nlg;ilg++){        
    unsigned int iug = ( ilg == nlg-1 ) ? 0 : ilg+1;    
    aTri[(itri1+ilg)*3+0] = nlt-2+(nlt-1)*ilg;
    aTri[(itri1+ilg)*3+1] = iu1;
    aTri[(itri1+ilg)*3+2] = nlt-2+(nlt-1)*iug;
  }
  const unsigned int itri2 = nlg*(nlt-2)*2+nlg;
	for(unsigned int ilg=0;ilg<nlg;ilg++){        
    unsigned int iug = ( ilg == nlg-1 ) ? 0 : ilg+1;    
    aTri[(itri2+ilg)*3+0] = (nlt-1)*iug;
    aTri[(itri2+ilg)*3+1] = iu2;
    aTri[(itri2+ilg)*3+2] = (nlt-1)*ilg;    
  }    
}


// return penetration depth (inside is positive)
double CSignedDistanceField3D_Cylinder::Projection
(double px, double py, double pz,
 double n[3]) const // normal outward
{
  double dd = dir_[0]*dir_[0] + dir_[1]*dir_[1] + dir_[2]*dir_[2];
  double pod = (px-cent_[0])*dir_[0]+(py-cent_[1])*dir_[1]+(pz-cent_[2])*dir_[2];
  double a = pod/dd;
  double hp[3] = { cent_[0]+a*dir_[0], cent_[1]+a*dir_[1], cent_[2]+a*dir_[2] };
	double d[3] = { px-hp[0], py-hp[1], pz-hp[2] };
	const double len = sqrt( d[0]*d[0]+d[1]*d[1]+d[2]*d[2] );		
	const double invlen = 1.0/len;
	if( !is_out_ ){
		n[0] = -d[0]*invlen;
		n[1] = -d[1]*invlen;		
		n[2] = -d[2]*invlen;
		return +len-radius_;
	}
	n[0] = d[0]*invlen;
	n[1] = d[1]*invlen;		
	n[2] = d[2]*invlen;	
	return radius_-len;
}

unsigned int CSignedDistanceField3D_Cylinder::FindInOut(double px, double py, double pz) const
{
	double n[3];
	double pd = this->Projection(px, py, pz, n);	
	if( !is_out_ ) pd *= -1.0;
	if( pd > 0 ){ return 0; }
	return 1;
}

bool CSignedDistanceField3D_Cylinder::IntersectionPoint
(double p[3], 
 const double o[3], const double d[3]) const 
{
  const double q[3] = { o[0]-cent_[0], o[1]-cent_[1], o[2]-cent_[2] };
  const double a = d[0]*d[0] + d[1]*d[1] + d[2]*d[2];
  const double b = q[0]*d[0] + q[1]*d[1] + q[2]*d[2];
  const double c = q[0]*q[0] + q[1]*q[1] + q[2]*q[2] - radius_*radius_;
  const double det = b*b-a*c;
  if( det < 0 )  return false;
  const double t = (-b+sqrt(det))/a;
  p[0] = o[0] + t*d[0];
  p[1] = o[1] + t*d[1];
  p[2] = o[2] + t*d[2];
  return true;
}

////////////////////////////////////////////////////////////////


CSignedDistanceField3D_Torus::CSignedDistanceField3D_Torus(){
	cent_[0] = 0;	cent_[1] = 0;	cent_[2] = 0;
	radius_ = 0.5;
	radius_tube_ = 0.2;
}

void CSignedDistanceField3D_Torus::GetMesh
(std::vector<unsigned int>& aTri,
 std::vector<double>& aXYZ,
 double elen) const
{
	const unsigned int nlg = 32;
	const unsigned int nlt = 18;	
	const double rlg = 6.28/nlg;	// longtitude
	const double rlt = 6.28/nlt;	// latitude
  aXYZ.resize(nlg*nlt*3);
	for(unsigned int ilg=0;ilg<nlg;ilg++){    
  for(unsigned int ilt=0;ilt<nlt;ilt++){
    aXYZ[(ilg*nlt+ilt)*3+0] = ( radius_ + radius_tube_*cos(ilt*rlt) )*sin(ilg*rlg);
    aXYZ[(ilg*nlt+ilt)*3+1] = ( radius_ + radius_tube_*cos(ilt*rlt) )*cos(ilg*rlg);
    aXYZ[(ilg*nlt+ilt)*3+2] = radius_tube_*sin(ilt*rlt);
	}
  }
  aTri.resize(nlg*nlt*2*3);
	for(unsigned int ilg=0;ilg<nlg;ilg++){        
  for(unsigned int ilt=0;ilt<nlt;ilt++){
    unsigned int iug = ( ilg == nlg-1 ) ? 0 : ilg+1;
    unsigned int iut = ( ilt == nlt-1 ) ? 0 : ilt+1; 
    aTri[(ilg*nlt+ilt)*6+0] = ilg*nlt+ilt;
    aTri[(ilg*nlt+ilt)*6+2] = iug*nlt+ilt;
    aTri[(ilg*nlt+ilt)*6+1] = iug*nlt+iut;
    ////
    aTri[(ilg*nlt+ilt)*6+3] = ilg*nlt+ilt;
    aTri[(ilg*nlt+ilt)*6+5] = iug*nlt+iut;
    aTri[(ilg*nlt+ilt)*6+4] = ilg*nlt+iut;
  }
  }
}

	
// return penetration depth (inside is positive)
double CSignedDistanceField3D_Torus::Projection
(double px, double py, double pz,
 double n[3]) const // normal outward
{
	double dir[3] = { px-cent_[0], py-cent_[1], pz-cent_[2] };
	const double t = dir[2];
//	dir[0] -= t*norm_[0];
//	dir[1] -= t*norm_[1];
	dir[2] -= t;		
	const double len = sqrt( dir[0]*dir[0]+dir[1]*dir[1]+dir[2]*dir[2] );
	if( len < 1.0e-20 ){	// vertical to the center of torus
		return radius_tube_ - radius_;
	}
	const double invlen = 1.0/len;
	dir[0] *= invlen;
	dir[1] *= invlen;		
	dir[2] *= invlen;
	double p[3];
	p[0] = cent_[0]+radius_*dir[0];		
	p[1] = cent_[1]+radius_*dir[1];		
	p[2] = cent_[2]+radius_*dir[2];		
	double dir2[3] = { px-p[0], py-p[1], pz-p[2] };
	const double len2 = sqrt( dir2[0]*dir2[0]+dir2[1]*dir2[1]+dir2[2]*dir2[2] );
	const double invlen2 = 1.0/len2;
	n[0] = dir2[0]*invlen2;
	n[1] = dir2[1]*invlen2;		
	n[2] = dir2[2]*invlen2;
	//		std::cout << len << " " << len2 << std::endl;
	return radius_tube_-len2;
}
	
unsigned int CSignedDistanceField3D_Torus::FindInOut(double px, double py, double pz) const
{
	double n[3];
	const double pd = this->Projection(px, py, pz, n);
	if( pd > 0 ){ return 0; }
	return 1;
}

/////


double CSignedDistanceField3D_Combine::Projection
(double px, double py, double pz,
 double n[3]) const // normal
{
  if( apCT.size() == 0 ) return -1;  
  double max_dist;
  max_dist = apCT[0]->Projection(px,py,pz, n);
  for(unsigned int ipct=1;ipct<apCT.size();ipct++){
    double dist0,n0[3];
    dist0 = apCT[ipct]->Projection(px,py,pz, n0);
    if( dist0 < max_dist ) continue;
    max_dist = dist0;
    n[0] = n0[0];
    n[1] = n0[1];
    n[2] = n0[2];
  }  
  return max_dist;
}

void CSignedDistanceField3D_Combine::GetMesh
(std::vector<unsigned int>& aTri,
 std::vector<double>& aXYZ,
 double elen) const
{
  aTri.size();
  aXYZ.size();
  for(unsigned int ipct=0;ipct<apCT.size();ipct++){  
    std::vector<unsigned int> aTri0;
    std::vector<double> aXYZ0;    
    apCT[ipct]->GetMesh(aTri0,aXYZ0,elen);
    const unsigned int i0 = aXYZ.size();
    for(unsigned int itri0=0;itri0<aTri0.size();itri0++){
      aTri.push_back(aTri0[itri0*3+0]+i0);
      aTri.push_back(aTri0[itri0*3+1]+i0);
      aTri.push_back(aTri0[itri0*3+2]+i0); 
    }
    for(unsigned int i=0;i<aXYZ0.size();i++){
      aXYZ.push_back(aXYZ0[i]);
    }
  }
}

/////



double CSignedDistanceField3D_Transform::Projection
(double px, double py, double pz,
 double n[3]) const // normal
{
	const double mat[3][3] = {
		{ cos(psi)*cos(theta),	cos(psi)*sin(theta)*sin(phi)-sin(psi)*cos(phi), cos(psi)*sin(theta)*cos(phi)+sin(psi)*sin(phi) },
		{ sin(psi)*cos(theta),	sin(psi)*sin(theta)*sin(phi)+cos(psi)*cos(phi), sin(psi)*sin(theta)*cos(phi)-cos(psi)*sin(phi) },
		{ -sin(theta),			cos(theta)*sin(phi),							cos(theta)*cos(phi) } };
	const double px0 = px-trans[0];
	const double py0 = py-trans[1];
	const double pz0 = pz-trans[2];	
	const double px1 = mat[0][0]*px0 + mat[1][0]*py0 + mat[2][0]*pz0;
	const double py1 = mat[0][1]*px0 + mat[1][1]*py0 + mat[2][1]*pz0;
	const double pz1 = mat[0][2]*px0 + mat[1][2]*py0 + mat[2][2]*pz0;
	double n0[3];
	const double d = this->pCT->Projection(px1,py1,pz1, n0);
	n[0] = mat[0][0]*n0[0] + mat[0][1]*n0[1] + mat[0][2]*n0[2];
	n[1] = mat[1][0]*n0[0] + mat[1][1]*n0[1] + mat[1][2]*n0[2];
	n[2] = mat[2][0]*n0[0] + mat[2][1]*n0[1] + mat[2][2]*n0[2];
	return d;
}


void CSignedDistanceField3D_Transform::GetMesh
(std::vector<unsigned int>& aTri,
 std::vector<double>& aXYZ,
 double elen) const
{
	const double mat[3][3] = {
		{ cos(psi)*cos(theta),	cos(psi)*sin(theta)*sin(phi)-sin(psi)*cos(phi), cos(psi)*sin(theta)*cos(phi)+sin(psi)*sin(phi) },
		{ sin(psi)*cos(theta),	sin(psi)*sin(theta)*sin(phi)+cos(psi)*cos(phi), sin(psi)*sin(theta)*cos(phi)-cos(psi)*sin(phi) },
		{ -sin(theta),			cos(theta)*sin(phi),							cos(theta)*cos(phi) } };
  this->pCT->GetMesh(aTri,aXYZ,elen);
  const unsigned int nnode = aXYZ.size()/3;
  for(unsigned int ino=0;ino<nnode;ino++){
    double px0 = aXYZ[ino*3+0];                     
    double py0 = aXYZ[ino*3+1];
    double pz0 = aXYZ[ino*3+2]; 
    const double px1 = mat[0][0]*px0 + mat[0][1]*py0 + mat[0][2]*pz0;
    const double py1 = mat[1][0]*px0 + mat[1][1]*py0 + mat[1][2]*pz0;
    const double pz1 = mat[2][0]*px0 + mat[2][1]*py0 + mat[2][2]*pz0;
    const double px2 = px1+trans[0];
    const double py2 = py1+trans[1];
    const double pz2 = pz1+trans[2];	    
    aXYZ[ino*3+0] = px2;
    aXYZ[ino*3+1] = py2;
    aXYZ[ino*3+2] = pz2;    
  }
}





CSignedDistanceField3D_Mesh::CSignedDistanceField3D_Mesh(){
  nnode_ = 0;	pXYZs_ = 0;
  ntri_ = 0;	aTri_ = 0;
  pBoxel_ = 0;
  is_hole = false;
}

CSignedDistanceField3D_Mesh::~CSignedDistanceField3D_Mesh(){
  if( pXYZs_  != 0 ){ delete pXYZs_; }
  if( aTri_   != 0 ){ delete aTri_; }
  if( pBoxel_ != 0 ){ delete pBoxel_; }
}

void CSignedDistanceField3D_Mesh::GetCenterWidth(double& cx, double& cy, double& cz,
                                                 double& wx, double& wy, double& wz)
{
  double x_min = pXYZs_[0], x_max = pXYZs_[0];
  double y_min = pXYZs_[1], y_max = pXYZs_[1];
  double z_min = pXYZs_[2], z_max = pXYZs_[2];
  for(unsigned int ino=1;ino<nnode_;ino++){
    x_min = ( x_min < pXYZs_[ino*3+0] ) ? x_min : pXYZs_[ino*3+0];
    x_max = ( x_max > pXYZs_[ino*3+0] ) ? x_max : pXYZs_[ino*3+0];
    y_min = ( y_min < pXYZs_[ino*3+1] ) ? y_min : pXYZs_[ino*3+1];
    y_max = ( y_max > pXYZs_[ino*3+1] ) ? y_max : pXYZs_[ino*3+1];
    z_min = ( z_min < pXYZs_[ino*3+2] ) ? z_min : pXYZs_[ino*3+2];
    z_max = ( z_max > pXYZs_[ino*3+2] ) ? z_max : pXYZs_[ino*3+2];
  }
  cx = (x_min+x_max)*0.5;
  cy = (y_min+y_max)*0.5;
  cz = (z_min+z_max)*0.5;
  wx = x_max-x_min;
  wy = y_max-y_min;
  wz = z_max-z_min;
}

/*
 void CSignedDistanceField3D_Mesh::Load_Off(const std::string& fname)
 {
	std::ifstream fin;
	fin.open(fname.c_str());
	if( fin.fail() ){
 std::cout << "Fail Read Fail" << std::endl;
 return;
	}
	std::string str;
	fin >> str;
	fin >> nnode_ >> ntri_ >> str;
	std::cout << "Load Off Nnode: ntri :" << nnode_ << " " << ntri_ << std::endl;
	pXYZs_ = new double [nnode_*3];
	aTri_ = new unsigned int [ntri_*3];
	for(unsigned int ino=0;ino<nnode_;ino++){
 double x,y,z;
 fin >> x >> y >> z;
 //		std::cout << ino << " " << x << " " << y << " " << z << std::endl;
 pXYZs_[ino*3+0] = x;
 pXYZs_[ino*3+1] = y;
 pXYZs_[ino*3+2] = z;
	}
	for(unsigned int itri=0;itri<ntri_;itri++){
 int itmp, i1, i2, i3;
 fin >> itmp >> i1 >> i2 >> i3;
 aTri_[itri*3+0] = i1;
 aTri_[itri*3+1] = i2;
 aTri_[itri*3+2] = i3;
 //		std::cout << itri << " " << itmp << " " << i1 << " " << i2 << " " << i3 << std::endl;
	}
 }
 
 void CSignedDistanceField3D_Mesh::Load_Gmv(const std::string& fname)
 {
	std::cout << "File load " << fname << std::endl;
	std::ifstream fin;
	fin.open(fname.c_str());
	if( fin.fail() ){
 std::cout << "Fail Read Fail" << std::endl;
 return;
	}
	std::string str;
	fin >> str;
	fin >> str;
	fin >> str;
	fin >> str;
	fin >> str;
	fin >> str;
	fin >> str;
	fin >> str;
	fin >> str;
	fin >> nnode_;
	std::cout << "Nnode " << nnode_ << std::endl;
	pXYZs_ = new double [nnode_*3];
 
	for(unsigned int ino=0;ino<nnode_;ino++){
 double x;
 fin >> x;
 pXYZs_[ino*3+0] = x;
	}
	for(unsigned int ino=0;ino<nnode_;ino++){
 double y;
 fin >> y;
 pXYZs_[ino*3+1] = y;
	}
	for(unsigned int ino=0;ino<nnode_;ino++){
 double z;
 fin >> z;
 pXYZs_[ino*3+2] = z;
	}
	
	fin >> str;
	fin >> ntri_;
	std::cout << "Ntri " << ntri_ << std::endl;
	aTri_ = new unsigned int [ntri_*3];
	for(unsigned int itri=0;itri<ntri_;itri++){
 int itmp, i1, i2, i3;
 fin >> str >> itmp >> i1 >> i2 >> i3;
 aTri_[itri*3+0] = i1-1;
 aTri_[itri*3+1] = i2-1;
 aTri_[itri*3+2] = i3-1;
 //			std::cout << itri << " " << itmp << " " << i1 << " " << i2 << " " << i3 << std::endl;
	}
 }
 
 void CSignedDistanceField3D_Mesh::Load_Ply(const std::string& fname)
 {
	std::cout << "File load " << fname << std::endl;
	std::ifstream fin;
	fin.open(fname.c_str());
	if( fin.fail() ){
 std::cout << "Fail Read Fail" << std::endl;
 return;
	}
 const unsigned int nbuff = 256;
 char buff[nbuff], buff1[nbuff], buff2[nbuff];
	std::string str1,str2;
 fin.getline(buff,nbuff);  // ply
 fin.getline(buff,nbuff);  // format asi 1.0
 for(;;){
 fin.getline(buff,nbuff);
 if( strncmp(buff, "comment ", 8) != 0 ){ break; }
 }
 /////
 sscanf(buff,"%s %s %d",buff1,buff2,&nnode_);
	std::cout << "Nnode " << nnode_ << std::endl;
 ////
 for(;;){
 fin.getline(buff,nbuff);
 if( strncmp(buff, "property ", 9) != 0 ){ break; }
 }
 sscanf(buff,"%s %s %d",buff1,buff2,&ntri_);
 std::cout << "NTri " << ntri_ << std::endl;
 /////
 fin.getline(buff,nbuff);  // property list int int vertex_indices
 fin.getline(buff,nbuff);  // end header
 ////
	pXYZs_ = new double [nnode_*3];
	for(unsigned int ino=0;ino<nnode_;ino++){
 double x,y,z;
 fin >> x >> y >> z;
 //		std::cout << ino << " " << x << " " << y << " " << z << std::endl;
 pXYZs_[ino*3+0] = x;
 pXYZs_[ino*3+1] = y;
 pXYZs_[ino*3+2] = z;
	}
	aTri_ = new unsigned int [ntri_*3];
	for(unsigned int itri=0;itri<ntri_;itri++){
 int itmp, i1, i2, i3;
 fin >> itmp >> i1 >> i2 >> i3;
 aTri_[itri*3+0] = i1;
 aTri_[itri*3+1] = i2;
 aTri_[itri*3+2] = i3;
 //		std::cout << itri << " " << itmp << " " << i1 << " " << i2 << " " << i3 << std::endl;
	}
 }
 */

void CSignedDistanceField3D_Mesh::SetMesh
(const std::vector<unsigned int>& aTri,
 const std::vector<double>& aXYZ)
{
  ntri_ = (int)aTri.size()/3;
  aTri_ = new unsigned int [ntri_*3];
  for(unsigned int i=0;i<ntri_*3;i++){ this->aTri_[i] = aTri[i]; }
  nnode_ = (int)aXYZ.size()/3;
  pXYZs_ = new double [nnode_*3];
  for(unsigned int i=0;i<nnode_*3;i++){ this->pXYZs_[i] = aXYZ[i]; }
}


// return penetration depth (inside is positive)
double CSignedDistanceField3D_Mesh::Projection
(double px, double py, double pz,
 double n[3]) const // normal outward
{
  unsigned int inout;
  if( pBoxel_ != 0 ){
    inout = this->FindInOut_Boxel(px,py,pz);
  }
  else {
    inout = this->FindInOut(px,py,pz);
  }
  ////
  double dist;
  {
    int itri;
    double r0,r1;
    dist = FindNearest(itri, r0, r1, px,py,pz);
  }
  if(      inout == 0 ){ return  dist; }
  else if( inout == 1 ){ return -dist; }
  return -dist;	// if not sure assume out
}

void CSignedDistanceField3D_Mesh::GetMesh
(std::vector<unsigned int>& aTri,
 std::vector<double>& aXYZ,
 double elen) const
{
  aTri.resize(ntri_*3);
  for(unsigned int i=0;i<ntri_*3;i++){ aTri[i] = this->aTri_[i]; }
  aXYZ.resize(nnode_*3);
  for(unsigned int i=0;i<nnode_*3;i++){ aXYZ[i] = this->pXYZs_[i]; }
}


// 0:in 1:out 2:not sure
double CSignedDistanceField3D_Mesh::FindNearest
(int& itri, double& r0, double& r1,
 double px, double py, double pz) const
{
  double dist = -1;
  const double ps[3] = {px,py,pz};
  if( pBoxel_ == 0 ){
    for (unsigned int jtri = 0; jtri<ntri_; ++jtri){
      const int i0 = aTri_[jtri*3+0];
      const int i1 = aTri_[jtri*3+1];
      const int i2 = aTri_[jtri*3+2];
      double pn[3], s0, s1;
      GetNearest_TrianglePoint3D(pn, s0, s1, ps, pXYZs_+i0*3, pXYZs_+i1*3, pXYZs_+i2*3);
      const double d1 = Distance3D(pn, ps);
      if (d1>=dist&&dist>0) continue;
      dist = d1;
      itri = jtri;
      r0 = s0;
      r1 = s1;
    }
  }
  else{
    pBoxel_->Find_NearestTriCand(ps,aIndTriCand);
    for(unsigned int jjtri=0;jjtri<aIndTriCand.size();jjtri++){
      const int jtri = aIndTriCand[jjtri];
      const int i1 = aTri_[jtri*3+0];
      const int i2 = aTri_[jtri*3+1];
      const int i3 = aTri_[jtri*3+2];
      double pn[3], s0, s1;
      GetNearest_TrianglePoint3D(pn, s0, s1, ps, pXYZs_+i1*3, pXYZs_+i2*3, pXYZs_+i3*3);
      const double d1 = Distance3D(pn, ps);
      if (d1>=dist&&dist>0) continue;
      dist = d1;
      itri = jtri;
      r0 = s0;
      r1 = s1;
    }
  }
//  std::cout << r0 << " " <<r1 << " " << dist << std::endl;
  return dist;
}

/*
// 0:in 1:out 2:not sure
double CSignedDistanceField3D_Mesh::Distance_Mesh_Boxel
(double px, double py, double pz,
 double n[3]) const
{
  assert( this->pBoxel_ != 0 );
  double p0[3] = {px,py,pz};
  double dist = pBoxel_->GetWidth()*3;
  
  pBoxel_->Find_NearestTriCand(p0,aIndTriCand);
  //	std::cout << aIndTriCand.size() << std::endl;
  
  for(unsigned int iitri=0;iitri<aIndTriCand.size();iitri++){
    unsigned int itri = aIndTriCand[iitri];
    assert( itri < ntri_ );
    const unsigned int i1 = aTri_[itri*3+0];
    const unsigned int i2 = aTri_[itri*3+1];
    const unsigned int i3 = aTri_[itri*3+2];
    double pn[3];
    GetNearest_TrianglePoint3D(pn, p0, pXYZs_+i1*3, pXYZs_+i2*3, pXYZs_+i3*3);
    double d1 = Distance3D(pn, p0);
    if (d1>=dist) continue;
    dist = d1;
    n[0] = pn[0]-p0[0];
    n[1] = pn[1]-p0[1];
    n[2] = pn[2]-p0[2];
  }
  {
    const double invlen = 1.0/Length3D(n);
    n[0] *= invlen;
    n[1] *= invlen;
    n[2] *= invlen;
  }
  return dist;
}
 */


// 0:in 1:out 2:not sure
unsigned int CSignedDistanceField3D_Mesh::FindInOut_IntersectionRay
(double px, double py, double pz,
 const double dir[3]) const
{
  double p0[3] = {px,py,pz};
  double p1[3] = {px+dir[0],py+dir[1],pz+dir[2]};
  unsigned int icnt = 0;
  for(unsigned int itri=0;itri<ntri_;itri++){
    unsigned int i1 = aTri_[itri*3+0];
    unsigned int i2 = aTri_[itri*3+1];
    unsigned int i3 = aTri_[itri*3+2];
    const double v0 = TetVolume3D(p0, pXYZs_+i1*3, pXYZs_+i2*3, pXYZs_+i3*3);
    const double sign = ( v0 > 0 ) ? 1 : -1;
    const double v1 = TetVolume3D(p0, pXYZs_+i2*3, pXYZs_+i3*3, p1)*sign;
    const double v2 = TetVolume3D(p0, pXYZs_+i3*3, pXYZs_+i1*3, p1)*sign;
    const double v3 = TetVolume3D(p0, pXYZs_+i1*3, pXYZs_+i2*3, p1)*sign;
    if( fabs(v1+v2+v3) < 1.0e-10 ) return 2;	// p0 and p1 is on the triangle
    double inv_v4 = 1.0/fabs(v1+v2+v3);
    const double r1 = v1*inv_v4;
    const double r2 = v2*inv_v4;
    const double r3 = v3*inv_v4;
    const double tol = 1.0e-2;
    if( r1 < -tol || r2 < -tol || r3 < -tol ) continue;	// need tol  ( compare with fabs(v1+v2+v3)? )
    if( r1 < tol || r2 < tol || r3 < tol ) return 2;	// on the edge
    const double dir2[3] = {
      pXYZs_[i1*3+0]*r1 + pXYZs_[i2*3+0]*r2 + pXYZs_[i3*3+0]*r3 - px,
      pXYZs_[i1*3+1]*r1 + pXYZs_[i2*3+1]*r2 + pXYZs_[i3*3+1]*r3 - py,
      pXYZs_[i1*3+2]*r1 + pXYZs_[i2*3+2]*r2 + pXYZs_[i3*3+2]*r3 - pz};
    double dotdir = Dot3D(dir,dir2);
    if( dotdir > 0 ) icnt++;
  }
  if( icnt % 2 == 0 ) return 1;
  return 0;
}


// 0:in 1:out 2:not sure
unsigned int CSignedDistanceField3D_Mesh::FindInOut_IntersectionRay_Boxel
(double px, double py, double pz,
 const double dir[3]) const
{
  assert( pBoxel_ != 0 );
  double p0[3] = {px,py,pz};
  pBoxel_->Find_IntersecTriCand(p0,dir,aIndTriCand);
  //	std::cout << aIndTriCand_.size() << std::endl;
  double p1[3] = {px+dir[0],py+dir[1],pz+dir[2]};
  unsigned int icnt = 0;
  for(unsigned int iitri=0;iitri<aIndTriCand.size();iitri++){
    unsigned int itri = aIndTriCand[iitri];
    if( aFlgTriUsed[itri] == 1 ) continue;
    aFlgTriUsed[itri] = 1;
    unsigned int i1 = aTri_[itri*3+0];
    unsigned int i2 = aTri_[itri*3+1];
    unsigned int i3 = aTri_[itri*3+2];
    const double v0 = TetVolume3D(p0, pXYZs_+i1*3, pXYZs_+i2*3, pXYZs_+i3*3);
    const double sign = ( v0 > 0 ) ? 1 : -1;
    const double v1 = TetVolume3D(p0, pXYZs_+i2*3, pXYZs_+i3*3, p1)*sign;
    const double v2 = TetVolume3D(p0, pXYZs_+i3*3, pXYZs_+i1*3, p1)*sign;
    const double v3 = TetVolume3D(p0, pXYZs_+i1*3, pXYZs_+i2*3, p1)*sign;
    if( fabs(v1+v2+v3) < 1.0e-10 ) goto AMBIGUOUS;	// p0 and p1 is on the triangle
    double inv_v4 = 1.0/fabs(v1+v2+v3);
    const double r1 = v1*inv_v4;
    const double r2 = v2*inv_v4;
    const double r3 = v3*inv_v4;
    const double tol = 1.0e-2;
    if( r1 < -tol || r2 < -tol || r3 < -tol ) continue;	// need tol  ( compare with fabs(v1+v2+v3)? )
    if( r1 < tol || r2 < tol || r3 < tol ) goto AMBIGUOUS;	// on the edge
    double dir2[3] = {
      pXYZs_[i1*3+0]*r1 + pXYZs_[i2*3+0]*r2 + pXYZs_[i3*3+0]*r3 - px,
      pXYZs_[i1*3+1]*r1 + pXYZs_[i2*3+1]*r2 + pXYZs_[i3*3+1]*r3 - py,
      pXYZs_[i1*3+2]*r1 + pXYZs_[i2*3+2]*r2 + pXYZs_[i3*3+2]*r3 - pz};
    double dotdir = Dot3D(dir,dir2);
    if( dotdir > 0 ) icnt++;
  }
  for(unsigned int iitri=0;iitri<aIndTriCand.size();iitri++){
    unsigned int itri = aIndTriCand[iitri];
    aFlgTriUsed[itri] = 0;
  }
  //	std::cout << "Cunt" << icnt << std::endl;
  if( icnt % 2 == 0 ) return 1;
  return 0;
AMBIGUOUS:
  for(unsigned int iitri=0;iitri<aIndTriCand.size();iitri++){
    unsigned int itri = aIndTriCand[iitri];
    aFlgTriUsed[itri] = 0;
  }
  return 2;
}


void CSignedDistanceField3D_Mesh::BuildBoxel()
{
  //	double c[3] = {0.0542,-0.04374532,0.06234};
  double c[3] = {0,0,0};
  double w[3] = {0,0,0};
  this->GetCenterWidth(c[0],c[1],c[2], w[0],w[1],w[2]);
  double width = w[0];
  width = ( w[1] > width ) ? w[1] : width;
  width = ( w[2] > width ) ? w[2] : width;
  if( pBoxel_ != 0 ){ delete pBoxel_; }
  pBoxel_ = new CSpatialHash_Grid3D(32,c,width*0.5*1.13454325);
  for(unsigned int itri=0;itri<ntri_;itri++){
    unsigned int i0 = aTri_[itri*3+0];
    unsigned int i1 = aTri_[itri*3+1];
    unsigned int i2 = aTri_[itri*3+2];
    pBoxel_->AddTri(itri, pXYZs_+i0*3, pXYZs_+i1*3, pXYZs_+i2*3);
  }
  if( !is_hole ){ pBoxel_->BuildOutFlg(); }
  aFlgTriUsed.clear();
  aFlgTriUsed.resize(ntri_,0);
  aIndTriCand.reserve(2048);
}

unsigned int CSignedDistanceField3D_Mesh::FindInOut(double px, double py, double pz) const
{
  unsigned int icnt_in  = 0;
  unsigned int icnt_out = 0;
  for(unsigned int i=0;i<10;i++){
    		const double theta = i*6.28/10+0.1431432154;
    		const double dir[3] = { sin(theta)*cos(theta*2), sin(theta)*sin(theta*2), cos(theta) };
//    const double dir[3] = { 1, 0, 0 };
    unsigned int ires = FindInOut_IntersectionRay(px,py,pz, dir);
    if( ires != 2 && !is_hole ){ return ires; }
    if( ires == 0 ) icnt_in++;
    if( ires == 1 ) icnt_out++;
  }
  if( icnt_in > 5 )  return 0;
  if( icnt_out > 5 ) return 1;
  /*	for(unsigned int i=1;i<20;i++){
   const double theta = i*6.28*543260;
   const double dir[3] = { sin(theta)*cos(theta*3), cos(theta), sin(theta)*sin(theta*3) };
   unsigned int ires = FindInOut_IntersectionRay(px,py,pz, dir);
   if( ires == 0 ) icnt_in++;
   if( ires == 1 ) icnt_out++;
   if( icnt_in  - icnt_out > 5 ) return 0;
   if( icnt_out - icnt_in  > 5 ) return 1;
   }*/
  return 2;
}

unsigned int CSignedDistanceField3D_Mesh::FindInOut_Boxel
(double px, double py, double pz) const
{
  assert( pBoxel_ != 0 );
  double p[3] = { px, py, pz };
  if( pBoxel_->IsOut(p) && !is_hole ){ return 1; }
  unsigned int icnt_in  = 0;
  unsigned int icnt_out = 0;
  for(unsigned int i=0;i<10;i++){
    const double theta = i*6.28/10+0.15432452356673;
    const double dir[3] = { sin(theta)*cos(theta*2), sin(theta)*sin(theta*2), cos(theta) };
    //		const double dir[3] = { -0.35, 0.1342, 0.3 };
    unsigned int ires  = FindInOut_IntersectionRay_Boxel(px,py,pz, dir);
    //		unsigned int ires1 = FindInOut_IntersectionRay(px,py,pz, dir);
    /*		if( ires != ires1 ){
     std::cout << "hoge " << px << " " << py << " " << pz << std::endl;
     }*/
    if( ires != 2 && !is_hole ){ return ires; }
    if( ires == 0 ) icnt_in++;
    if( ires == 1 ) icnt_out++;
  }
  if( icnt_in > 5 )  return 0;
  if( icnt_out > 5 ) return 1;
  return 2;
}

void CSignedDistanceField3D_Mesh::Translate(double x, double y, double z)
{
  for(unsigned int ino=0;ino<nnode_;ino++){
    pXYZs_[ino*3+0] += x;
    pXYZs_[ino*3+1] += y;
    pXYZs_[ino*3+2] += z;
  }
}


bool IsIntersecTri3D
(double& r0, double &r1, double psec[3],
 const double ps[3], const double dir[3],
 const double q0[3], const double q1[3], const double q2[3])
{
  double pe[3] = {ps[0]+dir[0],ps[1]+dir[1],ps[2]+dir[2]};
  const double v012 = TetVolume3D(ps, q0,q1,q2);
  const double sign = ( v012 > 0 ) ? 1 : -1;
  const double v0 = TetVolume3D(ps, q1, q2, pe)*sign;
  const double v1 = TetVolume3D(ps, q2, q0, pe)*sign;
  const double v2 = TetVolume3D(ps, q0, q1, pe)*sign;
  double inv_v4 = 1.0/fabs(v0+v1+v2);
  r0 = v0*inv_v4;
  r1 = v1*inv_v4;
  const double r2 = (1-r0-r1);
  const double tol = 1.0e-2;
  if( r0 < -tol || r1 < -tol || r2 < -tol ) return false;	// need tol  ( compare with fabs(v1+v2+v3)? )
  psec[0] = q0[0]*r0 + q1[0]*r1 + q2[0]*r2;
  psec[1] = q0[1]*r0 + q1[1]*r1 + q2[1]*r2;
  psec[2] = q0[2]*r0 + q1[2]*r1 + q2[2]*r2;
  double dir2[3] = { psec[0]-ps[0], psec[1]-ps[1], psec[2]-ps[2] };
  double dotdir = Dot3D(dir,dir2);
  if( dotdir < 0 ) return false;
  return true;
}

bool CSignedDistanceField3D_Mesh::FindIntersectionTri
(double psec[3], int& itri, double& r0, double& r1,
 const double org[3], const double dir[3]) const
{
  if( pBoxel_ != 0 ){
    std::vector<unsigned int> aIndTriCand;
    pBoxel_->Find_IntersecTriCand(org,dir, aIndTriCand);
    if( aIndTriCand.empty() ) return false;
    double min_dist = -1;
    for(unsigned int ijtri=0;ijtri<aIndTriCand.size();ijtri++){
      const int jtri = aIndTriCand[ijtri];
      const int i0 = aTri_[jtri*3+0];
      const int i1 = aTri_[jtri*3+1];
      const int i2 = aTri_[jtri*3+2];
      double s0, s1, qsec[3];
      bool res = IsIntersecTri3D(s0, s1, qsec, org, dir, pXYZs_+i0*3, pXYZs_+i1*3, pXYZs_+i2*3);
      if( !res ) continue;
      double dist = Distance3D(org, qsec);
      if( min_dist < 0 || dist < min_dist ){
        min_dist = dist;
        itri = jtri;
        r0 = s0;
        r1 = s1;
        psec[0] = qsec[0];
        psec[1] = qsec[1];
        psec[2] = qsec[2];
      }
    }
    if( min_dist >0 ) return true;
  }
  //////////////////////////////////////
  else {
    double min_dist = -1;
    for(unsigned int jtri=0;jtri<ntri_;jtri++){
      unsigned int i0 = aTri_[jtri*3+0];
      unsigned int i1 = aTri_[jtri*3+1];
      unsigned int i2 = aTri_[jtri*3+2];
      double s0, s1, qsec[3];
      bool res = IsIntersecTri3D(s0, s1, qsec, org, dir, pXYZs_+i0*3, pXYZs_+i1*3, pXYZs_+i2*3);
      if( !res ) continue;
      double dist = Distance3D(org, qsec);
      if( min_dist < 0 || dist < min_dist ){
        min_dist = dist;
        itri = jtri;
        r0 = s0;
        r1 = s1;
        psec[0] = qsec[0];
        psec[1] = qsec[1];
        psec[2] = qsec[2];
      }
    }
    if( min_dist >0 ) return true;
  }
  
  return false;
}


bool CSignedDistanceField3D_Mesh::IntersectionPoint
(double p[3],
 const double org[3], const double dir[3]) const
{
  double r0, r1;
  int itri;
  bool res = this->FindIntersectionTri(p,itri,r0,r1,org,dir);
  return res;
}


////////////////


