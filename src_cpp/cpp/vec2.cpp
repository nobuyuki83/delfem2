/*
 * Copyright (c) 2019 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */

#include <cstdlib>

#include "delfem2/vec2.h"

double TriArea2D(const double v1[2], const double v2[2], const double v3[2]){
  double z = ( v2[0] - v1[0] )*( v3[1] - v1[1] ) - ( v3[0] - v1[0] )*( v2[1] - v1[1] );
  return z*0.5;
}

double SqDistance2D(const double v1[2], const double v2[2]){
  return (v1[0]-v2[0])*(v1[0]-v2[0]) + (v1[1]-v2[1])*(v1[1]-v2[1]);
}

double Distance2D(const double v1[2], const double v2[2]){
  return sqrt( (v1[0]-v2[0])*(v1[0]-v2[0]) + (v1[1]-v2[1])*(v1[1]-v2[1]) );
}

static bool IsCrossLines(const double po_s0[], const double po_e0[],
                         const double po_s1[], const double po_e1[] )
{
  const double area1 = TriArea2D(po_s0,po_e0,po_s1);
  const double area2 = TriArea2D(po_s0,po_e0,po_e1);
  if( area1 * area2 > 0.0 ) return false;
  const double area3 = TriArea2D(po_s1,po_e1,po_s0);
  const double area4 = TriArea2D(po_s1,po_e1,po_e0);
  if( area3 * area4 > 0.0 ) return false;
  return true;
}

void noise2D(double noise[2])
{
  double a0 = rand()/(RAND_MAX+1.0);
  double a1 = rand()/(RAND_MAX+1.0);
  double x = sqrt(-2.0*log(a0))*cos(3.1415*2*a1);
  double y = sqrt(-2.0*log(a0))*sin(3.1415*2*a1);
  noise[0] = x;
  noise[1] = y;
}

bool InverseMat2(double invB[4], const double B[4])
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

void matMat2(double AB[4], const double A[4], const double B[4])
{
  AB[0*2+0] = A[0*2+0]*B[0*2+0]+A[0*2+1]*B[1*2+0];
  AB[0*2+1] = A[0*2+0]*B[0*2+1]+A[0*2+1]*B[1*2+1];
  AB[1*2+0] = A[1*2+0]*B[0*2+0]+A[1*2+1]*B[1*2+0];
  AB[1*2+1] = A[1*2+0]*B[0*2+1]+A[1*2+1]*B[1*2+1];
}

void MatVec2(double w[2], const double A[4], const double v[2])
{
  w[0] = A[0*2+0]*v[0]+A[0*2+1]*v[1];
  w[1] = A[1*2+0]*v[0]+A[1*2+1]*v[1];
}

void setNormalized2(double w[2])
{
  double l = sqrt(w[0]*w[0]+w[1]*w[1]);
  double invl = 1.0/l;
  w[0] *= invl;
  w[1] *= invl;
}

double dot2(const double w[2], const double v[2])
{
  return w[0]*v[0]+w[1]*v[1];
}

double len2(const double v[2]){
  return sqrt(v[0]*v[0]+v[1]*v[1]);
}

double sqLen2(const double v[2]){
  return v[0]*v[0]+v[1]*v[1];
}

void gramian2(double AtA[3], const double A[4])
{
  AtA[0] = A[0*2+0]*A[0*2+0]+A[1*2+0]*A[1*2+0];
  AtA[1] = A[0*2+0]*A[0*2+1]+A[1*2+0]*A[1*2+1];
  AtA[2] = AtA[1];
  AtA[3] = A[0*2+1]*A[0*2+1]+A[1*2+1]*A[1*2+1];
}

void VLVt2(double A[4], double l0, double l1, const double V[4])
{
  A[0] = l0*V[0]*V[0]+l1*V[1]*V[1];
  A[1] = l0*V[2]*V[0]+l1*V[3]*V[1];
  A[2] = l0*V[0]*V[2]+l1*V[1]*V[3];
  A[3] = l0*V[2]*V[2]+l1*V[3]*V[3];
}


void RotationalComponentOfMatrix2(double R[4], const double M[4])
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
      if (sqLen2(v0)>eps){ setNormalized2(v0); }
      v1[0] = G[0]-l0;
      v1[1] = G[2];
      if (sqLen2(v1)>eps){ setNormalized2(v1); }
    }
  }
  double V[4] = { v0[0], v1[0], v0[1], v1[1] };
  if (l0<eps){ l0 = 1; }
  if (l1<eps){ l1 = 1; }
  double il0 = 1.0/sqrt(l0);
  double il1 = 1.0/sqrt(l1);
  double invS[4]; VLVt2(invS, il0, il1, V);
  matMat2(R, A, invS);
}


void makeSplineLoop
(const std::vector<double>& aCV,
 std::vector<double>& aVecCurve)
{
  aVecCurve.resize(0);
  const int nCV = (int)aCV.size()/2;
  unsigned int ndiv = 5;
  for(int icv=0;icv<nCV;icv++){
    int icv0=icv;   if( icv0 >= nCV ){ icv0-=nCV; }
    int icv1=icv+1; if( icv1 >= nCV ){ icv1-=nCV; }
    int icv2=icv+2; if( icv2 >= nCV ){ icv2-=nCV; }
    const double p0[2] = { aCV[icv0*2+0], aCV[icv0*2+1] };
    const double p1[2] = { aCV[icv1*2+0], aCV[icv1*2+1] };
    const double p2[2] = { aCV[icv2*2+0], aCV[icv2*2+1] };
    for(unsigned int idiv=0;idiv<ndiv;idiv++){
      const double t = 1.0-(double)idiv/ndiv;
      const double w[3] = {0.5*t*t, -t*t + t + 0.5, 0.5*(1-t)*(1-t) };
      const double px = w[0]*p0[0] + w[1]*p1[0] + w[2]*p2[0];
      const double py = w[0]*p0[1] + w[1]*p1[1] + w[2]*p2[1];
      aVecCurve.push_back(px);
      aVecCurve.push_back(py);
    }
  }
}


void MeanValueCoordinate2D
(double* aW,
 double px, double py,
 const double* aXY, int nv)
{
  for(int iv=0;iv<nv;++iv){ aW[iv] = 0.0; }
  for(int iv=0;iv<nv;++iv){
    CVector2 v0(aXY[iv*2+0]-px,aXY[iv*2+1]-py);
    if( v0.Length() > 1.0e-10 ){ continue; }
    aW[iv] = 1.0;
    return;
  }
  for(int ie=0;ie<nv;++ie){
    int iv0 = (ie+0)%nv;
    int iv1 = (ie+1)%nv;
    CVector2 v0(aXY[iv0*2+0]-px,aXY[iv0*2+1]-py);
    CVector2 v1(aXY[iv1*2+0]-px,aXY[iv1*2+1]-py);
    const double l0 = v0.Length();
    const double l1 = v1.Length();
    if( abs((v0*v1)/(l0*l1)+1) > 1.0e-10 ){ continue; }
    aW[iv0] = l1/(l0+l1);
    aW[iv1] = l0/(l0+l1);
    return;
  }
  double sum = 0;
  for(int ie=0;ie<nv;++ie){
    int iv0 = (ie+0)%nv;
    int iv1 = (ie+1)%nv;
    int iv2 = (ie+2)%nv;
    CVector2 v0(aXY[iv0*2+0]-px,aXY[iv0*2+1]-py);
    CVector2 v1(aXY[iv1*2+0]-px,aXY[iv1*2+1]-py);
    CVector2 v2(aXY[iv2*2+0]-px,aXY[iv2*2+1]-py);
    double c01 = (v0*v1)/(v0.Length()*v1.Length());
    double s01 = (Cross(v0,v1)>0)?1:-1;
    double c12 = (v1*v2)/(v1.Length()*v2.Length());
    double s12 = (Cross(v1,v2)>0)?1:-1;
    double t01 = s01*sqrt((1-c01)/(1+c01));
    double t12 = s12*sqrt((1-c12)/(1+c12));
    double w1 =  (t01+t12)/v1.Length();
    aW[iv1] = w1;
    sum += w1;
  }
  for(int iv=0;iv<nv;++iv){
    aW[iv] /= sum;
  }
}

////////////////////////////////////////////////////////////////////////////////////////////////

std::ostream &operator<<(std::ostream &output, const CVector2& v)
{
  output.setf(std::ios::scientific);
  output<<v.x<<" "<<v.y;
  return output;
}

std::istream &operator>>(std::istream &input, CVector2& v)
{
  input>>v.x>>v.y;
  return input;
}

CVector2 operator*(double c, const CVector2& v0)
{
  return CVector2(v0.x*c,v0.y*c);
}

CVector2 operator*(const CVector2& v0, double c)
{
  return CVector2(v0.x*c,v0.y*c);
}

double operator * (const CVector2& lhs, const CVector2& rhs)
{
  return lhs.x*rhs.x + lhs.y*rhs.y;
}

double operator ^ (const CVector2& lhs, const CVector2& rhs)
{
  return lhs.x*rhs.y - lhs.y*rhs.x;
}

//! divide by real number
CVector2 operator/ (const CVector2& vec, double d)
{
  CVector2 temp = vec;
  temp /= d;
  return temp;
}

CVector2 rotate(const CVector2& p0, double theta)
{
  CVector2 p1;
  double c = cos(theta), s = sin(theta);
  p1.x = c*p0.x-s*p0.y;
  p1.y = s*p0.x+c*p0.y;
  return p1;
}

CVector2 rotate90(const CVector2& p0)
{
  return CVector2(-p0.y,p0.x);
}

///////////////////////////////////////////////////////////////////////////////////////////////


CVector2 Mat2Vec(const double A[4], const CVector2& v)
{
  CVector2 w;
  w.x = A[0]*v.x+A[1]*v.y;
  w.y = A[2]*v.x+A[3]*v.y;
  return w;
}

//! Area of the Triangle
double TriArea
(const CVector2& v1,
 const CVector2& v2,
 const CVector2& v3)
{
  return 0.5*( (v2.x-v1.x)*(v3.y-v1.y) - (v3.x-v1.x)*(v2.y-v1.y) );
}

inline double Cross(const CVector2& v1, const CVector2& v2){
  return v1.x*v2.y - v2.x*v1.y;
}

inline double SquareLength
(const CVector2& ipo0, const CVector2& ipo1)
{
  return  ( ipo1.x - ipo0.x )*( ipo1.x - ipo0.x ) + ( ipo1.y - ipo0.y )*( ipo1.y - ipo0.y );
}

double SquareLength(const CVector2& point){
  return  point.x*point.x + point.y*point.y;
}

double Length(const CVector2& point){
  return  point.Length();
}

// Distance between two points
double Distance
(const CVector2& ipo0,
 const CVector2& ipo1)
{
  return  sqrt( ( ipo1.x - ipo0.x )*( ipo1.x - ipo0.x ) + ( ipo1.y - ipo0.y )*( ipo1.y - ipo0.y ) );
}

// Distance between two points
double SquareDistance
(const CVector2& ipo0,
 const CVector2& ipo1)
{
  return  ( ipo1.x - ipo0.x )*( ipo1.x - ipo0.x ) + ( ipo1.y - ipo0.y )*( ipo1.y - ipo0.y );
}

// Hight of a triangle : between v1 and line of v2-v3
double TriHeight(const CVector2& v1, const CVector2& v2, const CVector2& v3){
  const double area = TriArea(v1,v2,v3);
  const double len = sqrt( SquareLength(v2,v3) );
  return area*2.0/len;
}

// compute dot product
double Dot(const CVector2& ipo0, const CVector2& ipo1)
{
  return  ipo0.x*ipo1.x + ipo0.y*ipo1.y;
}

// get parameter 't' of the line against point. t=0 is po_s, t=1 is po_e
// this one has implementation in header because GetDist_LineSeg_Point below refers this
double FindNearestPointParameter_Line_Point
(const CVector2& po_c,
 const CVector2& po_s, const CVector2& po_e)
{
  const CVector2& es = po_e-po_s;
  const CVector2& sc = po_s-po_c;
  const double a = SquareLength(es);
  const double b = Dot(es, sc);
  return -b/a;
}

CVector2 GetNearest_LineSeg_Point
(const CVector2& po_c,
 const CVector2& po_s, const CVector2& po_e)
{
  double t = FindNearestPointParameter_Line_Point(po_c, po_s, po_e);
  if (t < 0){ return po_s; }
  if (t > 1){ return po_e; }
  return po_s+t*(po_e-po_s);
}

double GetDist_LineSeg_Point
(const CVector2& po_c,
 const CVector2& po_s, const CVector2& po_e)
{
  CVector2 p = GetNearest_LineSeg_Point(po_c, po_s,po_e);
  return Distance(p, po_c);
}

bool IsCross_LineSeg_LineSeg
(const CVector2& po_s0, const CVector2& po_e0,
 const CVector2& po_s1, const CVector2& po_e1 )
{
  {
    const double min0x = ( po_s0.x < po_e0.x ) ? po_s0.x : po_e0.x;
    const double max0x = ( po_s0.x > po_e0.x ) ? po_s0.x : po_e0.x;
    const double max1x = ( po_s1.x > po_e1.x ) ? po_s1.x : po_e1.x;
    const double min1x = ( po_s1.x < po_e1.x ) ? po_s1.x : po_e1.x;
    const double min0y = ( po_s0.y < po_e0.y ) ? po_s0.y : po_e0.y;
    const double max0y = ( po_s0.y > po_e0.y ) ? po_s0.y : po_e0.y;
    const double max1y = ( po_s1.y > po_e1.y ) ? po_s1.y : po_e1.y;
    const double min1y = ( po_s1.y < po_e1.y ) ? po_s1.y : po_e1.y;
    const double len = ((max0x-min0x)+(max0y-min0y)+(max1x-min1x)+(max1y-min1y))*0.0001;
    //		std::cout << len << std::endl;
    if( max1x+len < min0x ) return false;
    if( max0x+len < min1x ) return false;
    if( max1y+len < min0y ) return false;
    if( max0y+len < min1y ) return false;
  }
  const double area1 = TriArea(po_s0,po_e0,po_s1);
  const double area2 = TriArea(po_s0,po_e0,po_e1);
  const double area3 = TriArea(po_s1,po_e1,po_s0);
  const double area4 = TriArea(po_s1,po_e1,po_e0);  
  //	std::cout << area1 << " " << area2 << " " << area3 << " " << area4 << std::endl;
  const double a12 = area1*area2; if( a12 > 0 ) return false;
  const double a34 = area3*area4; if( a34 > 0 ) return false;
  return true;
}

double GetDist_LineSeg_LineSeg
(const CVector2& po_s0, const CVector2& po_e0,
 const CVector2& po_s1, const CVector2& po_e1)
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
double SquareCircumradius
(const CVector2& p0,
 const CVector2& p1,
 const CVector2& p2 )
{
	const double area = TriArea(p0,p1,p2);
  
	const double dtmp0 = SquareLength(p1,p2);
	const double dtmp1 = SquareLength(p0,p2);
	const double dtmp2 = SquareLength(p0,p1);
  
	return dtmp0*dtmp1*dtmp2/(16.0*area*area);
}

//! center of the circumcircle
bool CenterCircumcircle
(const CVector2& p0,
 const CVector2& p1,
 const CVector2& p2,
 CVector2& center)
{
  const double area = TriArea(p0,p1,p2);
  if( fabs(area) < 1.0e-10 ){ return false; }
  const double tmp_val = 1.0/(area*area*16.0);
  
  const double dtmp0 = SquareLength(p1,p2);
  const double dtmp1 = SquareLength(p0,p2);
  const double dtmp2 = SquareLength(p0,p1);
  
  const double etmp0 = tmp_val*dtmp0*(dtmp1+dtmp2-dtmp0);
  const double etmp1 = tmp_val*dtmp1*(dtmp0+dtmp2-dtmp1);
  const double etmp2 = tmp_val*dtmp2*(dtmp0+dtmp1-dtmp2);
  
  center.x = etmp0*p0.x + etmp1*p1.x + etmp2*p2.x;
  center.y = etmp0*p0.y + etmp1*p1.y + etmp2*p2.y;
  return true;
}


// check if Delaunay condition satisfied
// 0 : p3 is inside circum circle on the p0,p1,p2
// 1 :       on
// 2 :       outsdie
int DetDelaunay
(const CVector2& p0,
 const CVector2& p1,
 const CVector2& p2,
 const CVector2& p3)
{
	const double area = TriArea(p0,p1,p2);
	if( fabs(area) < 1.0e-10 ){
		return 3;
	}
	const double tmp_val = 1.0/(area*area*16.0);
  
	const double dtmp0 = SquareLength(p1,p2);
	const double dtmp1 = SquareLength(p0,p2);
	const double dtmp2 = SquareLength(p0,p1);
  
	const double etmp0 = tmp_val*dtmp0*(dtmp1+dtmp2-dtmp0);
	const double etmp1 = tmp_val*dtmp1*(dtmp0+dtmp2-dtmp1);
	const double etmp2 = tmp_val*dtmp2*(dtmp0+dtmp1-dtmp2);
  
	const CVector2 out_center(etmp0*p0.x + etmp1*p1.x + etmp2*p2.x,
                            etmp0*p0.y + etmp1*p1.y + etmp2*p2.y );
  
	const double qradius = SquareLength(out_center,p0);
	const double qdistance = SquareLength(out_center,p3);
  
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

CVector2 pointCurve_BezierCubic
(double t,
 const CVector2& p1, const CVector2& p2, const CVector2& p3, const CVector2& p4)
{
  double tp = 1.0-t;
  return t*t*t*p4 + 3*t*t*tp*p3 + 3*t*tp*tp*p2 + tp*tp*tp*p1;
}


// ----------------------------------------------------------------------------------
// std::vector starts from here

//! Area of the Triangle (3 indexes and vertex array)
double TriArea
(const int iv1, const int iv2, const int iv3,
 const std::vector<CVector2>& point )
{
  return TriArea(point[iv1],point[iv2],point[iv3]);
}

void Polyline_CubicBezierCurve
(std::vector<CVector2>& aP,
 const int n, 
 const std::vector<CVector2>& aCP)
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

void Polyline_BezierCubic
(std::vector<CVector2>& aP,
 const unsigned int n,
 const CVector2& p1, const CVector2& p2, const CVector2& p3, const CVector2& p4)
{
  aP.resize(n);
  for(unsigned int i=0;i<n;++i){
    const double t = (double)i/(n-1.0);
    aP[i] = pointCurve_BezierCubic(t,
                                   p1, p2, p3, p4);
  }
}

double Length_Polygon
(const std::vector<CVector2>& aP)
{
  double len = 0;
  for(unsigned int ip0=0;ip0<aP.size()-1;ip0++){
    int ip1 = ip0+1;
    len += (aP[ip0]-aP[ip1]).Length();
  }
  return len;
}


double Area_Polygon(const std::vector<CVector2>& aP)
{
  CVector2 vtmp(0,0);
  const int ne = aP.size();
  double area_loop = 0.0;
  for(int ie=0;ie<ne;ie++){
    area_loop += TriArea(vtmp, aP[(ie+0)%ne], aP[(ie+1)%ne]);
  }
  return area_loop;
}

void MeanValueCoordinate
(std::vector<double>& aW,
 CVector2& p,
 std::vector<CVector2>& aVtx)
{
  const int nv = (int)aVtx.size();
  aW.assign(nv,0.0);
  double sum = 0;
  for(int ie=0;ie<nv;++ie){
    int iv0 = (ie+0)%nv;
    int iv1 = (ie+1)%nv;
    int iv2 = (ie+2)%nv;
    CVector2 v0 = aVtx[iv0]-p;
    CVector2 v1 = aVtx[iv1]-p;
    CVector2 v2 = aVtx[iv2]-p;
    double c01 = (v0*v1)/(v0.Length()*v1.Length());
    double c12 = (v1*v2)/(v1.Length()*v2.Length());
    double t01 = sqrt((1-c01)/(1+c01));
    double t12 = sqrt((1-c12)/(1+c12));
    double w1 =  (t01+t12)/v1.Length();
    aW[iv1] = w1;
    sum += w1;
  }
  for(int iv=0;iv<nv;++iv){
    aW[iv] /= sum;
  }
}



void makeRandomLoop
(unsigned int nCV,
 std::vector<double>& aCV)
{
  aCV.clear();
  for(unsigned int icv=0;icv<nCV;icv++){
    /*
     {
     aCV.push_back((double)rand()/(RAND_MAX+1.0));
     aCV.push_back((double)rand()/(RAND_MAX+1.0));
     }
     */
    { // polar coordinate random position
      double tht = icv*3.1415*2.0/nCV;
      double r = (double)rand()/(RAND_MAX+1.0);
      double px = r*sin(tht);
      double py = r*cos(tht);
      aCV.push_back(px);
      aCV.push_back(py);
    }
  }
}


void Translate
(std::vector<CVector2>& aP,
 double dx, double dy)
{
  for(unsigned int ip=0;ip<aP.size();++ip){
    aP[ip].x += dx;
    aP[ip].y += dy;
  }
}

void Rotate
(std::vector<CVector2>& aP,
 double dt)
{
  for(unsigned int ip=0;ip<aP.size();++ip){
    double x0 = aP[ip].x;
    double y0 = aP[ip].y;
    aP[ip].x = cos(dt)*x0 - sin(dt)*y0;
    aP[ip].y = sin(dt)*x0 + cos(dt)*y0;
  }
}

std::vector<CVector2> Polyline_Resample_Polyline
(const std::vector<CVector2>& stroke0,
 double l)
{
  if( stroke0.empty() ){
    std::vector<CVector2> a;
    return a;
  }
  std::vector<CVector2> stroke;
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

std::vector<CVector2> Polygon_Resample_Polygon
(const std::vector<CVector2>& stroke0,
 double l)
{
  std::vector<CVector2> stroke;
  if( stroke0.size() == 0 ) return stroke;
  stroke.push_back( stroke0[0] );
  int jcur = 0;
  double rcur = 0;
  double lcur = l;
  for(;;){
    if( jcur > (int)stroke0.size()-1 ) break;
    int ip0 = jcur;
    int ip1 = (jcur+1)%stroke0.size();
    double lenj  = (stroke0[ip1]-stroke0[ip0]).Length();
    double lenjr = lenj*(1.0-rcur);
    if( lenjr > lcur ){ // put point in this segment
      rcur += lcur/lenj;
      stroke.push_back( (1-rcur)*stroke0[ip0] + rcur*stroke0[ip1] );
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


void SecondMomentOfArea_Polygon
(CVector2& cg,  double& area,
 CVector2& pa1, double& I1,
 CVector2& pa2, double& I2,
 const std::vector<CVector2>& aVec2D)
{
  area = 0;
  const unsigned int nseg = aVec2D.size();
  cg = CVector2(0.0, 0.0);
  for(unsigned int iseg=0;iseg<nseg;iseg++){
    int ip0 = (iseg+0)%nseg;
    int ip1 = (iseg+1)%nseg;
    double x0 = aVec2D[ip0].x;
    double y0 = aVec2D[ip0].y;
    double x1 = aVec2D[ip1].x;
    double y1 = aVec2D[ip1].y;
    double ai = x0*y1 - x1*y0;
    area += ai;
    cg.x += ai*(x0+x1)/3.0;
    cg.y += ai*(y0+y1)/3.0;
  }
  cg.x /= area;
  cg.y /= area;
  area *= 0.5;
  ////////
  double Ix=0, Iy=0, Ixy=0;
  for(unsigned int iseg=0;iseg<nseg;iseg++){
    int ip0 = (iseg+0)%nseg;
    int ip1 = (iseg+1)%nseg;
    double x0 = aVec2D[ip0].x-cg.x;
    double y0 = aVec2D[ip0].y-cg.y;
    double x1 = aVec2D[ip1].x-cg.x;
    double y1 = aVec2D[ip1].y-cg.y;
    double ai = x0*y1 - x1*y0;
    Ix  += ai*(y0*y0 + y0*y1 + y1*y1)/12.0;
    Iy  += ai*(x0*x0 + x0*x1 + x1*x1)/12.0;
    Ixy += ai*(x0*y0 + 0.5*x0*y1 + 0.5*x1*y0 + x1*y1)/12.0;
  }
  if( fabs(Ix-Iy)+fabs(Ixy) < 1.0e-20 ){
    pa1 = CVector2(1,0);
    pa2 = CVector2(0,1);
    return;
  }
  double phi = 0.5*atan2(-2*Ixy,Ix-Iy);
  pa1 = CVector2(+cos(phi), +sin(phi));
  pa2 = CVector2(-sin(phi), +cos(phi));
  I1 = 0.5*(Ix+Iy)+0.5*sqrt( (Ix-Iy)*(Ix-Iy) + 4*Ixy*Ixy );
  I2 = 0.5*(Ix+Iy)-0.5*sqrt( (Ix-Iy)*(Ix-Iy) + 4*Ixy*Ixy );
}




void JArray_FromVecVec_XY
(std::vector<int>& aIndXYs,
 std::vector<int>& loopIP0,
 std::vector<CVector2>& aXY,
 const std::vector< std::vector<double> >& aVecAry0)
{
  const int nloop = (int)aVecAry0.size();
  aIndXYs.resize(nloop+1);
  aIndXYs[0] = 0;
  int npo_sum = 0;
  for(int iloop=0;iloop<(int)nloop;iloop++){
    const int npo = (int)aVecAry0[iloop].size()/2;
    aIndXYs[iloop+1] = aIndXYs[iloop]+npo;
    npo_sum += npo;
  }
  aXY.resize(npo_sum);
  npo_sum = 0;
  for(int iloop=0;iloop<(int)nloop;iloop++){
    const int nxys = (int)aVecAry0[iloop].size()/2;
    for(int ixys=0;ixys<nxys;ixys++){
      aXY[npo_sum].x = aVecAry0[iloop][ixys*2+0];
      aXY[npo_sum].y = aVecAry0[iloop][ixys*2+1];
      npo_sum++;
    }
  }
  loopIP0.resize(aXY.size());
  for(unsigned int ip=0;ip<aXY.size();++ip){ loopIP0[ip] = ip; }
}


void ResamplingLoop
(std::vector<int>& loopIP1_ind,
 std::vector<int>& loopIP1,
 std::vector<CVector2>& aVec2,
 double max_edge_length)
{
  assert( aVec2.size() == loopIP1.size() );
  const std::vector<int> loopIP0_ind = loopIP1_ind;
  const std::vector<int> loopIP0 = loopIP1;
  const int nloop = loopIP0_ind.size()-1;
  std::vector< std::vector<int> > aPoInEd(loopIP0.size());
  {
    for(int iloop=0;iloop<nloop;++iloop){
      const int np = loopIP0_ind[iloop+1]-loopIP0_ind[iloop];
      for(int ip=0;ip<np;ip++){
        const int iipo0 = loopIP0_ind[iloop]+(ip+0)%np; assert( iipo0>=0 && iipo0<(int)loopIP0.size() );
        const int iipo1 = loopIP0_ind[iloop]+(ip+1)%np; assert( iipo1>=0 && iipo1<(int)loopIP0.size() );
        const int ipo0 = loopIP0[iipo0]; assert(ipo0>=0&&ipo0<(int)aVec2.size());
        const int ipo1 = loopIP0[iipo1]; assert(ipo1>=0&&ipo1<(int)aVec2.size());
        const CVector2 po0 = aVec2[ipo0]; // never use reference here because aVec2 will resize afterward
        const CVector2 po1 = aVec2[ipo1]; // never use reference here because aVec2 will resize afterward
        const int nadd = (int)( Distance( po0, po1 ) / max_edge_length);
        if( nadd == 0 ) continue;
        for(int iadd=0;iadd<nadd;++iadd){
          double r2 = (double)(iadd+1)/(nadd+1);
          CVector2 v2 = (1-r2)*po0 + r2*po1;
          const int ipo2 = aVec2.size();
          aVec2.push_back(v2);
          assert( iipo0>=0 && iipo0<(int)aPoInEd.size() );
          aPoInEd[ iipo0 ].push_back(ipo2);
        }
      }
    }
  }
  ////
  loopIP1_ind.resize(nloop+1);
  loopIP1_ind[0] = 0;
  for(int iloop=0;iloop<nloop;++iloop){
    const int nbar0 = loopIP0_ind[iloop+1]-loopIP0_ind[iloop];
    int nbar1 = nbar0;
    for(int ibar=0;ibar<nbar0;ibar++){
      const int iip_loop = loopIP0_ind[iloop]+ibar;
      nbar1 += aPoInEd[iip_loop].size();
    }
    loopIP1_ind[iloop+1] = loopIP1_ind[iloop] + nbar1;
  }
  // adding new vertices on the outline
  loopIP1.resize(loopIP1_ind[nloop]);
  unsigned int ivtx0 = 0;
  for(int iloop=0;iloop<nloop;iloop++){
    for(int iip_loop=loopIP0_ind[iloop];iip_loop<loopIP0_ind[iloop+1];iip_loop++){
      const int ip_loop = loopIP0[iip_loop];
      loopIP1[ivtx0] = ip_loop;
      ivtx0++;
      for(unsigned int iadd=0;iadd<aPoInEd[ip_loop].size();iadd++){
        loopIP1[ivtx0] = aPoInEd[iip_loop][iadd];
        ivtx0++;
      }
    }
  }
  assert( loopIP1.size() == aVec2.size() );
  assert( loopIP1.size() == ivtx0 );
}

void MakeMassMatrixTri
(double M[9],
 double rho,
 const unsigned int aIP[3],
 const std::vector<CVector2>& aVec2)
{
  assert( aIP[0]<aVec2.size() );
  assert( aIP[1]<aVec2.size() );
  assert( aIP[2]<aVec2.size() );
  const double P[3][2] = {
    {aVec2[aIP[0]].x,aVec2[aIP[0]].y},
    {aVec2[aIP[1]].x,aVec2[aIP[1]].y},
    {aVec2[aIP[2]].x,aVec2[aIP[2]].y} };
  const double Area = TriArea2D(P[0], P[1], P[2]);
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

bool IsInclude_Loop
(const double co[],
 const int ixy_stt, const int ixy_end,
 const std::vector<CVector2>& aXY)
{
  int inum_cross = 0;
  for(int itr=0;itr<10;itr++){
    const double dir[2] = { cos((itr+1)*23.0), sin((itr+1)*23.0) }; // random direction
    const double codir[2] = { co[0]+dir[0], co[1]+dir[1] };
    bool is_fail = false;
    inum_cross = 0;
    for(int ixys=ixy_stt;ixys<ixy_end;ixys++){
      const int ipo0 = ixys;
      int ipo1 = ixys+1;
      if( ipo1 == ixy_end ){ ipo1 = ixy_stt; }
      const double p0[2] = {aXY[ipo0].x,aXY[ipo0].y};
      const double p1[2] = {aXY[ipo1].x,aXY[ipo1].y};
      const double area0 = TriArea2D(co,codir,p0);
      const double area1 = TriArea2D(co,p1,codir);
      double r1 =  area0 / (area0+area1);
      double r0 =  area1 / (area0+area1);
      if( fabs(area0+area1) < 1.0e-20 ){
        is_fail = true;
        break;
      }
      if( fabs(r0) < 1.0e-3 || fabs(r1) < 1.0e-3 ){
        is_fail = true;
        break;
      }
      if( r0*r1 < 0 ){ continue; }
      double po2[2] = {
        r0*aXY[ipo0].x+r1*aXY[ipo1].x,
        r0*aXY[ipo0].y+r1*aXY[ipo1].y };
      //      const double area2 = TriArea2D(co,codir,po2);
      double d = (po2[0]-co[0])*dir[0] + (po2[1]-co[1])*dir[1];
      if( d > 0 ) inum_cross++;
    }
    if( is_fail ){ continue; }
    if( inum_cross%2 == 0 ){ return false; }
    else if( inum_cross%2 == 1 ){ return true; }
  }
  return false;
}




bool CheckInputBoundaryForTriangulation
(const std::vector<int>& loopIP_ind,
 const std::vector<CVector2>& aXY)
{
  ////////////////////////////////
  // enter Input check section
  
  const int nloop = loopIP_ind.size()-1;
  
  { // make sure every loop has at least 3 points
    for(int iloop=0;iloop<nloop;iloop++){
      if( loopIP_ind[iloop+1]-loopIP_ind[iloop] < 3 ) return false;
    }
  }
  {
    ////////////////////////////////
    // check inclusion of loops
    for(int iloop=1;iloop<nloop;iloop++){
      for(int ipo=loopIP_ind[iloop];ipo<loopIP_ind[iloop+1];ipo++){
        const double pi[2] = {aXY[ipo].x,aXY[ipo].y};
        if( !IsInclude_Loop(pi,
                            loopIP_ind[0],loopIP_ind[1],
                            aXY) ) return false;
      }
    }
    // check inclusion
    for(int iloop=1;iloop<nloop;iloop++){
      for(int jloop=0;jloop<nloop;jloop++){
        if( iloop == jloop ) continue;
        for(int jpo=loopIP_ind[jloop];jpo<loopIP_ind[jloop+1];jpo++){
          const double pj[2] = {aXY[jpo].x,aXY[jpo].y};
          if( IsInclude_Loop(pj,
                             loopIP_ind[iloop],loopIP_ind[iloop+1],
                             aXY) ) return false;
        }
      }
    }
  }
  { // check intersection
    bool is_intersect = false;
    for(int iloop=0;iloop<nloop;iloop++){
      const int nei = loopIP_ind[iloop+1]-loopIP_ind[iloop];
      for(int ie=0;ie<nei;ie++){
        const int i0 = loopIP_ind[iloop] + (ie+0)%nei;
        const int i1 = loopIP_ind[iloop] + (ie+1)%nei;
        const double pi0[2] = {aXY[i0].x,aXY[i0].y};
        const double pi1[2] = {aXY[i1].x,aXY[i1].y};
        const double xmax_i = ( pi1[0] > pi0[0] ) ? pi1[0] : pi0[0];
        const double xmin_i = ( pi1[0] < pi0[0] ) ? pi1[0] : pi0[0];
        const double ymax_i = ( pi1[1] > pi0[1] ) ? pi1[1] : pi0[1];
        const double ymin_i = ( pi1[1] < pi0[1] ) ? pi1[1] : pi0[1];
        for(int je=ie+1;je<nei;je++){
          const int j0 = loopIP_ind[iloop] + (je+0)%nei;
          const int j1 = loopIP_ind[iloop] + (je+1)%nei;
          if( i0 == j0 || i0 == j1 || i1 == j0 || i1 == j1 ){ continue; }
          const double pj0[2] = {aXY[j0].x,aXY[j0].y};
          const double pj1[2] = {aXY[j1].x,aXY[j1].y};
          const double xmax_j = ( pj1[0] > pj0[0] ) ? pj1[0] : pj0[0];
          const double xmin_j = ( pj1[0] < pj0[0] ) ? pj1[0] : pj0[0];
          const double ymax_j = ( pj1[1] > pj0[1] ) ? pj1[1] : pj0[1];
          const double ymin_j = ( pj1[1] < pj0[1] ) ? pj1[1] : pj0[1];
          if( xmin_j > xmax_i || xmax_j < xmin_i ){ continue; }
          if( ymin_j > ymax_i || ymax_j < ymin_i ){ continue; }
          if( IsCrossLines(pi0,pi1,  pj0,pj1) ){
            is_intersect = true;
            break;
          }
        }
        if( is_intersect ) break;
        for(int jloop=iloop+1;jloop<nloop;jloop++){
          const int nbar_j = loopIP_ind[jloop+1]-loopIP_ind[jloop];
          for(int jbar=0;jbar<nbar_j;jbar++){
            const int jpo0 = loopIP_ind[jloop] + jbar;
            int jpo1 = loopIP_ind[jloop] + jbar+1;
            if( jbar == nbar_j-1 ){ jpo1 = loopIP_ind[jloop]; }
            const double pj0[2] = {aXY[jpo0].x,aXY[jpo0].y};
            const double pj1[2] = {aXY[jpo1].y,aXY[jpo1].y};
            const double xmax_j = ( pj1[0] > pj0[0] ) ? pj1[0] : pj0[0];
            const double xmin_j = ( pj1[0] < pj0[0] ) ? pj1[0] : pj0[0];
            const double ymax_j = ( pj1[1] > pj0[1] ) ? pj1[1] : pj0[1];
            const double ymin_j = ( pj1[1] < pj0[1] ) ? pj1[1] : pj0[1];
            if( xmin_j > xmax_i || xmax_j < xmin_i ) continue;  // åçˆÇ™Ç†ÇËÇ¶Ç»Ç¢ÉpÉ^Å[ÉìÇèúäO
            if( ymin_j > ymax_i || ymax_j < ymin_i ) continue;  // è„Ç…ìØÇ∂
            if( IsCrossLines(pi0,pi1,  pj0,pj1) ){
              is_intersect = true;
              break;
            }
          }
          if( is_intersect ) break;
        }
        if( is_intersect ) break;
      }
      if( is_intersect ) break;
    }
    if( is_intersect ) return false;
  }
  // end of input check section
  ////////////////////////////////////////////////
  return true;
}

std::vector<CVector2>
Polygon_Invert(const std::vector<CVector2>& aP)
{
  std::vector<CVector2> res;
  for(int ip=aP.size()-1;ip>=0;--ip){
    res.push_back(aP[ip]);
  }
  return res;
}

std::vector<double>
XY_Polygon(const std::vector<CVector2>& aP)
{
  std::vector<double> res;
  res.reserve(aP.size()*2);
  for(unsigned int ip=0;ip<aP.size();++ip){
    res.push_back(aP[ip].x);
    res.push_back(aP[ip].y);
  }
  return res;
}



void FixLoopOrientation
(std::vector<int>& loopIP,
 const std::vector<int>& loopIP_ind,
 const std::vector<CVector2>& aXY)
{
  const std::vector<int> loop_old = loopIP;
  const int nloop = loopIP_ind.size()-1;
  int ivtx0 = 0;
  for(int iloop=0;iloop<nloop;iloop++){
    double area_loop = 0;
    { // area of this loop
      CVector2 vtmp(0,0);
      const int nbar = loopIP_ind[iloop+1]-loopIP_ind[iloop];
      for(int ibar=0;ibar<nbar;ibar++){
        const int iipo0 = loopIP_ind[iloop]+(ibar+0)%nbar;
        const int iipo1 = loopIP_ind[iloop]+(ibar+1)%nbar;
        const int ipo0 = loop_old[iipo0];
        const int ipo1 = loop_old[iipo1];
        area_loop += TriArea(vtmp, aXY[ipo0], aXY[ipo1]);
      }
    }
    const int nbar0 = loopIP_ind[iloop+1]-loopIP_ind[iloop];
    if( (area_loop > 0) == (iloop == 0) ){ // outer loop
      for(int ibar=0;ibar<nbar0;ibar++){
        const int iipo = loopIP_ind[iloop] + ibar;
        const int ipo = loop_old[iipo];
        loopIP[ivtx0] = ipo;
        ivtx0++;
      }
    }
    else{
      for(int ibar=0;ibar<nbar0;ibar++){ // inner loop
        const int iipo = loopIP_ind[iloop+1] - 1 - ibar;
        const int ipo = loop_old[iipo];
        loopIP[ivtx0] = ipo;
        ivtx0++;
      }
    }
  }
}




