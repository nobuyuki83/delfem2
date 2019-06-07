/*
 * Copyright (c) 2019 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */


#include <cstdlib>
#include <cmath>
#include <map>
#include <stack>

#include "delfem2/vec3.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

static bool isnan_vector3(double x) { return x!=x; }

double ScalarTripleProduct3D(const double a[], const double b[], const double c[]){
  return a[0]*(b[1]*c[2] - b[2]*c[1])
  +a[1]*(b[2]*c[0] - b[0]*c[2])
  +a[2]*(b[0]*c[1] - b[1]*c[0]);
}

double Dot3D(const double a[], const double b[]){
  return a[0]*b[0]+a[1]*b[1]+a[2]*b[2];
}

void Cross3D(double r[3], const double v1[3], const double v2[3]){
  r[0] = v1[1]*v2[2] - v2[1]*v1[2];
  r[1] = v1[2]*v2[0] - v2[2]*v1[0];
  r[2] = v1[0]*v2[1] - v2[0]*v1[1];
}

double Length3D(const double v[3]){
  return sqrt( v[0]*v[0] + v[1]*v[1] + v[2]*v[2] );
}

void Normalize3D(double v[3]){
  double len = sqrt( v[0]*v[0] + v[1]*v[1] + v[2]*v[2] );
  v[0] /= len;
  v[1] /= len;
  v[2] /= len;
}

double SquareLength3D(const double v[3]){
  return v[0]*v[0] + v[1]*v[1] + v[2]*v[2];
}

double SquareDistance3D(const double p0[3], const double p1[3]){
  return (p1[0]-p0[0])*(p1[0]-p0[0]) + (p1[1]-p0[1])*(p1[1]-p0[1]) + (p1[2]-p0[2])*(p1[2]-p0[2]);
}

double Distance3D(const double p0[3], const double p1[3]){
  return sqrt( (p1[0]-p0[0])*(p1[0]-p0[0]) + (p1[1]-p0[1])*(p1[1]-p0[1]) + (p1[2]-p0[2])*(p1[2]-p0[2]) );
}

double TriArea3D(const double v1[3], const double v2[3], const double v3[3]){
  double x, y, z;
  x = ( v2[1] - v1[1] )*( v3[2] - v1[2] ) - ( v3[1] - v1[1] )*( v2[2] - v1[2] );
  y = ( v2[2] - v1[2] )*( v3[0] - v1[0] ) - ( v3[2] - v1[2] )*( v2[0] - v1[0] );
  z = ( v2[0] - v1[0] )*( v3[1] - v1[1] ) - ( v3[0] - v1[0] )*( v2[1] - v1[1] );
  return 0.5*sqrt( x*x + y*y + z*z );
}
double TetVolume3D
(const double v1[3],
 const double v2[3],
 const double v3[3],
 const double v4[3])
{
  return
  ((v2[0]-v1[0])*((v3[1]-v1[1])*(v4[2]-v1[2])-(v4[1]-v1[1])*(v3[2]-v1[2]))
   -(v2[1]-v1[1])*((v3[0]-v1[0])*(v4[2]-v1[2])-(v4[0]-v1[0])*(v3[2]-v1[2]))
   +(v2[2]-v1[2])*((v3[0]-v1[0])*(v4[1]-v1[1])-(v4[0]-v1[0])*(v3[1]-v1[1]))
   ) * 0.16666666666666666666666666666667;
}

void  UnitNormalAreaTri3D(double n[3], double& a, const double v1[3], const double v2[3], const double v3[3]){
  n[0] = ( v2[1] - v1[1] )*( v3[2] - v1[2] ) - ( v3[1] - v1[1] )*( v2[2] - v1[2] );
  n[1] = ( v2[2] - v1[2] )*( v3[0] - v1[0] ) - ( v3[2] - v1[2] )*( v2[0] - v1[0] );
  n[2] = ( v2[0] - v1[0] )*( v3[1] - v1[1] ) - ( v3[0] - v1[0] )*( v2[1] - v1[1] );
  a = sqrt(n[0]*n[0]+n[1]*n[1]+n[2]*n[2])*0.5;
  const double invlen = 0.5/a;
  n[0]*=invlen;	n[1]*=invlen;	n[2]*=invlen;
}

void NormalTri3D(double n[3], const double v1[3], const double v2[3], const double v3[3]){
  n[0] = ( v2[1] - v1[1] )*( v3[2] - v1[2] ) - ( v3[1] - v1[1] )*( v2[2] - v1[2] );
  n[1] = ( v2[2] - v1[2] )*( v3[0] - v1[0] ) - ( v3[2] - v1[2] )*( v2[0] - v1[0] );
  n[2] = ( v2[0] - v1[0] )*( v3[1] - v1[1] ) - ( v3[0] - v1[0] )*( v2[1] - v1[1] );
}

double VolumeTet3D(const double v1[3],
                    const double v2[3],
                    const double v3[3],
                    const double v4[3] )
{
  return
  (   ( v2[0] - v1[0] )*( ( v3[1] - v1[1] )*( v4[2] - v1[2] ) - ( v4[1] - v1[1] )*( v3[2] - v1[2] ) )
   -  ( v2[1] - v1[1] )*( ( v3[0] - v1[0] )*( v4[2] - v1[2] ) - ( v4[0] - v1[0] )*( v3[2] - v1[2] ) )
   +  ( v2[2] - v1[2] )*( ( v3[0] - v1[0] )*( v4[1] - v1[1] ) - ( v4[0] - v1[0] )*( v3[1] - v1[1] ) )
   ) * 0.16666666666666666666666666666667;
}

// t is a tmporary buffer size of 9
void InverseMat3(double Ainv[], const double A[])
{
  const double det =
  + A[0]*A[4]*A[8] + A[3]*A[7]*A[2] + A[6]*A[1]*A[5]
  - A[0]*A[7]*A[5] - A[6]*A[4]*A[2] - A[3]*A[1]*A[8];
  const double inv_det = 1.0/det;
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

// t is a tmporary buffer size of 9
void transposeMat3(double t[], const double a[])
{
  t[0] = a[0];
  t[1] = a[3];
  t[2] = a[6];
  t[3] = a[1];
  t[4] = a[4];
  t[5] = a[7];
  t[6] = a[2];
  t[7] = a[5];
  t[8] = a[8];
}

/////////////////////////////////////////

void GetNearest_LineSegPoint3D
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
void GetNearest_TrianglePoint3D
(double pn[3],
 double& r0, double& r1,
 const double ps[3], // origin point
 const double q0[3],
 const double q1[3],
 const double q2[3])
{
  double area, n012[3]; UnitNormalAreaTri3D(n012, area, q0, q1, q2);
  const double pe[3] = { ps[0]+n012[0], ps[1]+n012[1], ps[2]+n012[2] };
  const double v012 = VolumeTet3D(ps, q0, q1, q2);
  if (fabs(v012) > 1.0e-10){
    const double sign = (v012 > 0) ? +1 : -1;
    const double v0 = VolumeTet3D(ps, q1, q2, pe)*sign;
    const double v1 = VolumeTet3D(ps, q2, q0, pe)*sign;
    const double v2 = VolumeTet3D(ps, q0, q1, pe)*sign;
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



void GetVertical2Vector3D
(const double vec_n[3],
 double vec_x[3], double vec_y[3])
{
  const double vec_s[3] = {0,1,0};
  Cross3D(vec_x,vec_s,vec_n);
  const double len = Length3D(vec_x);
  if( len < 1.0e-10 ){
    const double vec_t[3] = {1,0,0};
    Cross3D(vec_x,vec_t,vec_n);  // z????
    Cross3D(vec_y,vec_n,vec_x);  // x????
  }
  else{
    const double invlen = 1.0/len;
    vec_x[0] *= invlen;
    vec_x[1] *= invlen;
    vec_x[2] *= invlen;
    Cross3D(vec_y,vec_n,vec_x);
  }
}

void GetRotMatrix_Rodrigues3D
(double rot[9],
 const double n[3], double theta)
{
  const double ct = cos(theta);
  const double st = sin(theta);
  rot[0] = ct+(1-ct)*n[0]*n[0];
  rot[1] =    (1-ct)*n[1]*n[0]-st*n[2];
  rot[2] =    (1-ct)*n[2]*n[0]+st*n[1];
  rot[3] =    (1-ct)*n[0]*n[1]+st*n[2];
  rot[4] = ct+(1-ct)*n[1]*n[1];
  rot[5] =    (1-ct)*n[2]*n[1]-st*n[0];
  rot[6] =    (1-ct)*n[0]*n[2]-st*n[1];
  rot[7] =    (1-ct)*n[1]*n[2]+st*n[0];
  rot[8] = ct+(1-ct)*n[2]*n[2];
}

void VecMat3D(const double x[3], const double m[9],  double y[3]){
  y[0] = m[0]*x[0] + m[3]*x[1] + m[6]*x[2];
  y[1] = m[1]*x[0] + m[4]*x[1] + m[7]*x[2];
  y[2] = m[2]*x[0] + m[5]*x[1] + m[8]*x[2];
}

void MatVec3(double y[3], const double m[9], const double x[3]){
  y[0] = m[0]*x[0] + m[1]*x[1] + m[2]*x[2];
  y[1] = m[3]*x[0] + m[4]*x[1] + m[5]*x[2];
  y[2] = m[6]*x[0] + m[7]*x[1] + m[8]*x[2];
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////

double Dot(const CVector3 &arg1, const CVector3 &arg2)
{
  return arg1.x*arg2.x + arg1.y*arg2.y + arg1.z*arg2.z;
}

// cross product
CVector3 Cross(const CVector3& arg1, const CVector3& arg2)
{
  CVector3 temp;
  temp.x = arg1.y*arg2.z - arg1.z*arg2.y;
  temp.y = arg1.z*arg2.x - arg1.x*arg2.z;
  temp.z = arg1.x*arg2.y - arg1.y*arg2.x;
  return temp;
}

//! add
CVector3 operator+ (const CVector3& lhs, const CVector3& rhs){
  CVector3 temp = lhs;
  temp += rhs;
  return temp;
}

//! subtract
CVector3 operator - (const CVector3& lhs, const CVector3& rhs){
  CVector3 temp = lhs;
  temp -= rhs;
  return temp;
}

//! scale
CVector3 operator* (double d, const CVector3& rhs){
  CVector3 temp = rhs;
  temp *= d;
  return temp;
}

//! scale
CVector3 operator* (const CVector3& vec, double d){
  CVector3 temp = vec;
  temp *= d;
  return temp;
}

//! mult
double operator* (const CVector3& lhs, const CVector3& rhs){
  return Dot(lhs,rhs);
}


//! divide by real number
CVector3 operator/ (const CVector3& vec, double d){
  CVector3 temp = vec;
  temp /= d;
  return temp;
}

//! mult
CVector3 operator^ (const CVector3& lhs, const CVector3& rhs){
  return Cross(lhs,rhs);
}

std::ostream &operator<<(std::ostream &output, const CVector3& v)
{
  output.setf(std::ios::scientific);
  output << v.x << " " << v.y << " " << v.z;
  return output;
}

std::istream &operator>>(std::istream &input, CVector3& v)
{
  input >> v.x >> v.y >> v.z;
  return input;
}


std::ostream &operator<<(std::ostream &output, const std::vector<CVector3>& aV){
  output<<aV.size()<<std::endl;
  for (int iv = 0; iv<(int)aV.size(); ++iv){
    output<<"  "<<aV[iv]<<std::endl;
  }
  return output;
}

std::istream &operator>>(std::istream &input, std::vector<CVector3>& aV){
  int nV;  input>>nV; aV.resize(nV);
  for (int iv = 0; iv<nV; iv++){
    input>>aV[iv];
  }
  return input;
}

CVector3 MatVec(double mat[9], const CVector3& v){
  CVector3 u;
  u.x = mat[0]*v.x + mat[1]*v.y + mat[2]*v.z;
  u.y = mat[3]*v.x + mat[4]*v.y + mat[5]*v.z;
  u.z = mat[6]*v.x + mat[7]*v.y + mat[8]*v.z;
  return u;
}

////////////////////////////////////////////

double ScalarTripleProduct(const CVector3& a, const CVector3& b, const CVector3& c){
  return a.x*(b.y*c.z - b.z*c.y) + a.y*(b.z*c.x - b.x*c.z) + a.z*(b.x*c.y - b.y*c.x);
}


bool operator== (const CVector3& lhs, const CVector3& rhs){
  if( fabs(lhs.x - rhs.x) < NEARLY_ZERO
     && fabs(lhs.y - rhs.y) < NEARLY_ZERO
     && fabs(lhs.z - rhs.z) < NEARLY_ZERO ){ return true; }
  else{ return false; }
}

bool operator!= (const CVector3& lhs, const CVector3& rhs){
  if( lhs == rhs )	return false;
  else return true;
}

void CVector3::SetNormalizedVector()
{
  double invmag = 1.0/Length();
  x *= invmag;
  y *= invmag;
  z *= invmag;
}

void CVector3::SetZero()
{
  x = 0.0;
  y = 0.0;
  z = 0.0;
}

//! Hight of a tetrahedra
double Height(const CVector3& v1, const CVector3& v2, const CVector3& v3, const CVector3& v4){
  // get normal vector
  double dtmp_x = (v2.y-v1.y)*(v3.z-v1.z)-(v2.z-v1.z)*(v3.y-v1.y);
  double dtmp_y = (v2.z-v1.z)*(v3.x-v1.x)-(v2.x-v1.x)*(v3.z-v1.z);
  double dtmp_z = (v2.x-v1.x)*(v3.y-v1.y)-(v2.y-v1.y)*(v3.x-v1.x);
  
  // normalize normal vector
  const double dtmp1 = 1.0 / sqrt( dtmp_x*dtmp_x + dtmp_y*dtmp_y + dtmp_z*dtmp_z );
  dtmp_x *= dtmp1;
  dtmp_y *= dtmp1;
  dtmp_z *= dtmp1;
  return (v4.x-v1.x)*dtmp_x+(v4.y-v1.y)*dtmp_y+(v4.z-v1.z)*dtmp_z;
}

////////////////////////////////////////////////////////////////

void GetVertical2Vector
(const CVector3& vec_n,
 CVector3& vec_x, CVector3& vec_y)
{
  vec_x = ::Cross(CVector3(0,1,0),vec_n);
  const double len = vec_x.Length();
  if( len < 1.0e-10 ){
    vec_x = ::Cross(CVector3(1,0,0),vec_n);  // z????
    vec_x.SetNormalizedVector();
    vec_y = ::Cross(vec_n,vec_x);  // x????
  }
  else{
    const double invlen = 1.0/len;
    vec_x *= invlen;
    vec_y = ::Cross(vec_n,vec_x);
  }
}

////////////////////////////////////////////////////////////////

CVector3 nearest_Line_Point
(const CVector3& p, // point
 const CVector3& s, // source
 const CVector3& d) // direction
{
  assert( Dot(d,d) > 1.0e-20 );
  const CVector3 ps = s-p;
  double a = Dot(d,d);
  double b = Dot(d,s-p);
  double t = -b/a;
  return s+t*d;
}

CVector3 nearest_Line_Point
(double& t,
 const CVector3& p, // point
 const CVector3& s, // source
 const CVector3& d) // direction
{
  if( Dot(d,d) < 1.0e-20 ){
    t = 0;
    return s;
  }
  const CVector3 ps = s-p;
  double a = Dot(d,d);
  double b = Dot(d,s-p);
  t = -b/a;
  return s+t*d;
}

CVector3 nearest_Origin_LineSeg
(const CVector3& s, // start
 const CVector3& e) // end
{
  CVector3 d = e-s;
  double a = Dot(d,d);
  if( a < 1.0e-20 ){ return (s+e)*0.5; }
  double b = Dot(d,s);
  double t = -b/a;
  if( t < 0 ) t = 0;
  if( t > 1 ) t = 1;
  return s+t*d;
}

// r0==0 -> p0==org
// r0==1 -> p1==org
CVector3 nearest_Origin_LineSeg
(double& r0,
 const CVector3& p0, // start
 const CVector3& p1) // end
{
  CVector3 d = p1-p0;
  double a = Dot(d,d);
  if( a < 1.0e-20 ){
    r0=0.5;
    return (p0+p1)*0.5;
  }
  double b = Dot(d,p0);
  r0 = -b/a;
  if( r0 < 0 ) r0 = 0;
  if( r0 > 1 ) r0 = 1;
  return (1.0-r0)*p0+r0*p1;
}

CVector3 nearest_LineSeg_Point
(const CVector3& p, // point
 const CVector3& s, // start
 const CVector3& e) // end
{
  CVector3 d = e-s;
  if( Dot(d,d) < 1.0e-20 ){
    return (s+e)*0.5;
  }
  const CVector3 ps = s-p;
  double a = Dot(d,d);
  double b = Dot(d,s-p);
  double t = -b/a;
  if( t < 0 ) t = 0;
  if( t > 1 ) t = 1;
  return s+t*d;
}

CVector3 nearest_LineSeg_Point
(double& t,
 const CVector3& p, // point
 const CVector3& s, // source
 const CVector3& e) // end
{
  CVector3 d = e-s;
  if( Dot(d,d) < 1.0e-20 ){
    t = 0.5;
    return (1-t)*s + t*e;
  }
  const CVector3 ps = s-p;
  double a = Dot(d,d);
  double b = Dot(d,s-p);
  t = -b/a;
  if( t < 0 ) t = 0;
  if( t > 1 ) t = 1;
  return s+t*d;
}

// Da,Db is nearest point scaled by D
void nearest_Line_Line
(double& D, CVector3& Da, CVector3& Db,
 const CVector3& pa_, const CVector3& va,
 const CVector3& pb_, const CVector3& vb)
{
  double xaa = va*va;
  double xab = vb*va;
  double xbb = vb*vb;
  D = (xaa*xbb-xab*xab);
  double xac = va*(pb_-pa_);
  double xbc = vb*(pb_-pa_);
  double da = xbb*xac-xab*xbc;
  double db = xab*xac-xaa*xbc;
  Da = D*pa_+da*va;
  Db = D*pb_+db*vb;
}

void nearest_Line_Line
(double& D, CVector3& Da, CVector3& Db, double& ta, double& tb,
 const CVector3& pa_, const CVector3& va,
 const CVector3& pb_, const CVector3& vb)
{
  double xaa = va*va;
  double xab = vb*va;
  double xbb = vb*vb;
  D = (xaa*xbb-xab*xab);
  double xac = va*(pb_-pa_);
  double xbc = vb*(pb_-pa_);
  ta = xbb*xac-xab*xbc;
  tb = xab*xac-xaa*xbc;
  Da = D*pa_+ta*va;
  Db = D*pb_+tb*vb;
}

CVector3 nearest_Plane_Point
(const CVector3& p, // point
 const CVector3& o, // origin
 const CVector3& n) // normal
{
  const CVector3 n0  = n.Normalize();
  return p + ((o-p)*n0)*n0;
}

CVector3 nearest_Orgin_PlaneTri
(double& r0,
 double& r1,
 const CVector3& q0,
 const CVector3& q1,
 const CVector3& q2)
{
  const CVector3 n1 = ((q1-q0)^(q2-q0)).Normalize();
  const double v0 = volume_OrgTet(q1, q2, n1);
  const double v1 = volume_OrgTet(q2, q0, n1);
  const double v2 = volume_OrgTet(q0, q1, n1);
  assert( fabs(v0+v1+v2) > 1.0e-10 );
  double vt_inv = 1.0/(v0+v1+v2);
  double r2;
  r0 = v0*vt_inv;
  r1 = v1*vt_inv;
  r2 = v2*vt_inv;
  return q0*r0 + q1*r1 + q2*r2;
}

CVector3 nearest_Origin_Tri
(double& r0,
 double& r1,
 const CVector3& q0,
 const CVector3& q1,
 const CVector3& q2)
{
  { // check on triangle
    CVector3 p012 = nearest_Orgin_PlaneTri(r0,r1, q0,q1,q2);
    if( r0>0 && r1>0 && (1-r0-r1)>0 ){ return p012; }
  }
  //////
  CVector3 p_min = q0;
  double d_min = q0.Length();
  double r2=0;
  r0=1; r1=0;
  {
    CVector3 p23 = nearest_Origin_LineSeg(r2, q1, q2);
    const double d23 = p23.Length();
    if(d23<d_min){ d_min=d23; p_min=p23; r1=1-r2; r0=0; }
  }
  {
    CVector3 p31 = nearest_Origin_LineSeg(r0, q2, q0);
    const double d31 = p31.Length();
    if(d31<d_min){ d_min=d31; p_min=p31; r2=1-r0; r1=0; }
  }
  {
    CVector3 p12 = nearest_Origin_LineSeg(r1, q0, q1);
    const double d12 = p12.Length();
    if(d12<d_min){ d_min=d12; p_min=p12; r0=1-r1; r2=0; }
  }
  /*
  {
    CVector3 p23 = nearest_Origin_LineSeg(r1, q1, q2);
    const double d23 = p23.Length();
    if(d23<d_min){ d_min=d23; p_min=p23; r2=1-r1; r0=0; }
  }
  {
    CVector3 p31 = nearest_Origin_LineSeg(r2, q2, q0);
    const double d31 = p31.Length();
    if(d31<d_min){ d_min=d31; p_min=p31; r0=1-r2; r1=0; }
  }
  {
    CVector3 p12 = nearest_Origin_LineSeg(r0, q0, q1);
    const double d12 = p12.Length();
    if(d12<d_min){ d_min=d12; p_min=p12; r1=1-r0; r2=0; }
  }
   */
  return p_min;
}

CVector3 nearst_Origin_Quad
(double& s0, double& s1,
 const CVector3& q0, const CVector3& q1, const CVector3& q2, const CVector3& q3)
{
  double dist_min = -1;
  CVector3 q_min;
  for(int ip=0;ip<5;++ip){
    double t0=0, t1=0;
    if(      ip == 0 ){ t0=0.0; t1=0.0; }
    else if( ip == 1 ){ t0=1.0; t1=0.0; }
    else if( ip == 2 ){ t0=1.0; t1=1.0; }
    else if( ip == 3 ){ t0=0.0; t1=1.0; }
    else if( ip == 4 ){ t0=0.5; t1=0.5; }
    CVector3 q;
    for(int itr=0;itr<4;++itr){
      CVector3 pq = (1-t0)*(1-t1)*q0 + t0*(1-t1)*q1 + t0*t1*q2 + (1-t0)*t1*q3;
      CVector3 dqt0 = -(1-t1)*q0 + (1-t1)*q1 + t1*q2 - t1*q3;
      CVector3 dqt1 = -(1-t0)*q0 - t0*q1 + t0*q2 + (1-t0)*q3;
      CVector3 ddqt0t1 = q0 - q1 + q2 - q3;
      double f0 = -dqt0*pq;
      double f1 = -dqt1*pq;
      double A00 = dqt0*dqt0;
      double A11 = dqt1*dqt1;
      double A01 = dqt1*dqt0 + ddqt0t1*pq;
      double det = A00*A11-A01*A01;
      double detinv = 1.0/det;
      double B00 = +A11*detinv;
      double B11 = +A00*detinv;
      double B01 = -A01*detinv;
      double d0 = B00*f0 + B01*f1;
      double d1 = B01*f0 + B11*f1;
      t0 += d0;
      t1 += d1;
    }
    double tol = 1.0e-4;
    if( t0 > -tol && t0 < 1.0+tol && t1 > -tol && t1 < 1.0+tol ){
      double d0 = q.Length();
      if( dist_min < 0 || d0 < dist_min ){
        dist_min = d0;
        s0 = t0;
        s1 = t1;
        q_min = q;
      }
    }
  }
  if( dist_min > 0 ){ return q_min; }
  ////
  const CVector3 q01 = nearest_Origin_LineSeg(q0,q1);
  const double d01  = q01.Length();
  if( dist_min < 0 || d01 < dist_min ){
    dist_min = d01;
    s0 = Distance(q01,q0)/Distance(q0,q1);
    s1 = 0.0;
    q_min = q01;
  }
  ////
  CVector3 q12 = nearest_Origin_LineSeg(q1,q2);
  const double d12  = q12.Length();
  if( dist_min < 0 || d12 < dist_min ){
    dist_min = d12;
    s0 = 1.0;
    s1 = Distance(q12,q1)/Distance(q1,q2);
    q_min = q12;
  }
  ////
  CVector3 q23 = nearest_Origin_LineSeg(q2,q3);
  const double d23  = q23.Length();
  if( dist_min < 0 || d23 < dist_min ){
    dist_min = d23;
    s0 = Distance(q23,q3)/Distance(q2,q3);
    s1 = 1.0;
    q_min = q23;
  }
  ////
  CVector3 q30 = nearest_Origin_LineSeg(q3,q0);
  const double d30  = q30.Length();
  if( dist_min < 0 || d30 < dist_min ){
    dist_min = d30;
    s0 = 0.0;
    s1 = Distance(q30,q0)/Distance(q3,q0);
    q_min = q30;
  }
  return q_min;
}

/*
CVector3 nearest_Origin_Tet
(double& r0, double& r1, double& r2,
 const CVector3& q0,
 const CVector3& q1,
 const CVector3& q2,
 const CVector3& q3)
{
  CVector3 p_min = q0;
  {
    bool res = barycentricCoord_Origin_Tet(r0, r1, r2, q0, q1, q2, q3);
    p_min = r0*q0 + r1*q1 + r2*q2 + (1-r0-r1-r2)*q3;
    if( r0>0 && r1>0 && r2>0 && (1-r0-r1-r2)>0 ){ return p_min; }
  }
  ////////////////////////
  double r3;
  { // face123
    r0 = 0;
    p_min = nearest_Orgin_PlaneTri(r1,r2,r3, q1,q2,q3);
    if( r1>0 && r2>0 && r3>0 ){ return p_min; }
  }
  { // face230
    r1 = 0;
    p_min = nearest_Orgin_PlaneTri(r2,r3,r0, q2,q3,q0);
    if( r2>0 && r3>0 && r0>0 ){ return p_min; }
  }
  { // face301
    r2 = 0;
    p_min = nearest_Orgin_PlaneTri(r3,r0,r1, q3,q0,q1);
    if( r3>0 && r0>0 && r1>0 ){ return p_min; }
  }
  { // face012
    r3 = 0;
    p_min = nearest_Orgin_PlaneTri(r0,r1,r2, q0,q1,q2);
    if( r0>0 && r1>0 && r2>0 ){ return p_min; }
  }
  ////////////////////////
  double d_min = q0.Length();
  double s0,s1,s2,s3;
  { // edge01
    CVector3 p01 = nearest_Origin_LineSeg(s0,s1,q0,q1);
    double d01 = p01.Length();
    if( d01<d_min ){ d_min=d01; p_min=p01; r0=s0; r1=s1; r2=0; r3=0; }
  }
  { // edge02
    CVector3 p02 = nearest_Origin_LineSeg(s0,s2,q0,q2);
    double d02 = p02.Length();
    if( d02<d_min ){ d_min=d02; p_min=p02; r0=s0; r1=0; r2=s2; r3=0; }
  }
  { // edge03
    CVector3 p03 = nearest_Origin_LineSeg(s0,s3,q0,q3);
    double d03 = p03.Length();
    if( d03<d_min ){ d_min=d03; p_min=p03; r0=s0; r1=0; r2=0; r3=s3; }
  }
  { // edge12
    CVector3 p12 = nearest_Origin_LineSeg(s1,s2,q1,q2);
    double d12 = p12.Length();
    if( d12<d_min ){ d_min=d12; p_min=p12; r0=0; r1=s1; r2=s2; r3=0; }
  }
  { // edge13
    CVector3 p13 = nearest_Origin_LineSeg(s1,s3,q1,q3);
    double d13 = p13.Length();
    if( d13<d_min ){ d_min=d13; p_min=p13; r0=0; r1=s1; r2=0; r3=s3; }
  }
  { // edge23
    CVector3 p23 = nearest_Origin_LineSeg(s2,s3,q2,q3);
    double d23 = p23.Length();
    if( d23<d_min ){ d_min=d23; p_min=p23; r0=0; r1=0; r2=s2; r3=s3; }
  }
  return p_min;
}
 */

void Nearest_Line_Circle
(CVector3& p0,
 CVector3& q0,
 const CVector3& src,
 const CVector3& dir,
 const CVector3& org, // center of the circle
 const CVector3& normal, // normal of the circle
 double rad)
{
  const int nitr = 4;
  ////
  CVector3 ex,ey; GetVertical2Vector(normal, ex, ey);
  double u0;
  {
    if( fabs(dir*normal)>fabs((org-src)*normal)*1.0e-4 ){
      u0 = ((org-src)*normal)/(dir*normal);
    }
    else{
      u0 = (org-src)*dir/(dir*dir);
    }
  }
  for(int itr=0;itr<nitr;++itr){
    p0 = src+u0*dir;
    double t0 = atan2(ey*(p0-org),ex*(p0-org));
    q0 = (rad*cos(t0))*ex + (rad*sin(t0))*ey + org;
    u0 = (q0-src)*dir/(dir*dir);
  }
}

////////////////////////////////////////////////////////////////

bool intersection_Plane_Line
(CVector3& p0, double& r0, double& r1, double& r2,
 double eps,
 const CVector3& src, const CVector3& dir,
 const CVector3& q0, const CVector3& q1, const CVector3& q2)
{
  r0 = volume_Tet(src, src+dir, q1, q2);
  r1 = volume_Tet(src, src+dir, q2, q0);
  r2 = volume_Tet(src, src+dir, q0, q1);
  double v012 = (r0+r1+r2);
  double v012_inv = 1.0/v012;
  r0 *= v012_inv;
  r1 *= v012_inv;
  r2 *= v012_inv;
  p0 = r0*q0 + r1*q1 + r2*q2;
  if( r0 > eps && r1 > eps && r2 > eps ){
    return true;
  }
  return false;
}

CVector3 intersection_Plane_Line
(const CVector3& o, // one point on plane
 const CVector3& n, // plane normal
 const CVector3& s, // one point on line
 const CVector3& d) // direction of line
{
  double t = ((o-s)*n)/(d*n);
  return s + t*d;
}

void iteration_intersection_Line_Quad
(double& t0, double& t1,
 const CVector3& src, const CVector3& u, const CVector3& v,
 const CVector3& q0, const CVector3& q1, const CVector3& q2, const CVector3& q3)
{
  CVector3 q = (1-t0)*(1-t1)*q0 + t0*(1-t1)*q1 + t0*t1*q2 + (1-t0)*t1*q3;
  CVector3 pq = q-src;
  CVector3 dqt0 = -(1-t1)*q0 + (1-t1)*q1 + t1*q2 - t1*q3;
  CVector3 dqt1 = -(1-t0)*q0 - t0*q1 + t0*q2 + (1-t0)*q3;
  CVector3 ddqt0t1 = q0 - q1 + q2 - q3;
  double f0 = -u*pq;
  double f1 = -v*pq;
  double A00 = u*dqt0;
  double A01 = u*dqt1;
  double A10 = v*dqt0;
  double A11 = v*dqt1;
  double det = A00*A11-A01*A10;
  double detinv = 1.0/det;
  double B00 = +A11*detinv;
  double B01 = -A01*detinv;
  double B10 = -A10*detinv;
  double B11 = +A00*detinv;
  double d0 = B00*f0 + B01*f1;
  double d1 = B10*f0 + B11*f1;
  t0 += d0;
  t1 += d1;
}

bool intersection_Point_Quad
(CVector3& psec, double& s0, double& s1,
 const CVector3& src, const CVector3& dir,
 const CVector3& q0, const CVector3& q1, const CVector3& q2, const CVector3& q3)
{
  CVector3 u,v; GetVertical2Vector(dir, u, v);
  ////
  double dist_min = -1;
  CVector3 q_min;
  for(int ip=0;ip<5;++ip){
    double t0=0, t1=0;
    if(      ip == 0 ){ t0=0.0; t1=0.0; }
    else if( ip == 1 ){ t0=1.0; t1=0.0; }
    else if( ip == 2 ){ t0=1.0; t1=1.0; }
    else if( ip == 3 ){ t0=0.0; t1=1.0; }
    else if( ip == 4 ){ t0=0.5; t1=0.5; }
    for(int itr=0;itr<4;++itr){
      iteration_intersection_Line_Quad(t0,t1,
                                       src, u, v,
                                       q0,q1,q2,q3);
    }
    CVector3 q = (1-t0)*(1-t1)*q0 + t0*(1-t1)*q1 + t0*t1*q2 + (1-t0)*t1*q3;
    double tol = 1.0e-4;
//    std::cout << t0 << " " << t1 << std::endl;
    if( t0 > -tol && t0 < 1.0+tol && t1 > -tol && t1 < 1.0+tol ){
      double d0 = (q-src).Length();
      if( dist_min < 0 || d0 < dist_min ){
        dist_min = d0;
        s0 = t0;
        s1 = t1;
        q_min = q;
      }
    }
  }
//  std::cout << dist_min << std::endl;
  if( dist_min > 0 ){
    psec = q_min;
    return true;
  }
  return false;
}

////////////////////////////////////////////////////////////////


CVector3 positionBarycentricCoord_Pyramid
(double r0,
 double r1,
 double r2,
 const CVector3& p0,
 const CVector3& p1,
 const CVector3& p2,
 const CVector3& p3,
 const CVector3& p4)
{
  return (1.0-r2)*(1.0-r0)*(1.0-r1)*p0
  + (1.0-r2)*r0*(1.0-r1)*p1
  + (1.0-r2)*r0*r1*p2
  + (1.0-r2)*(1.0-r0)*r1*p3
  + r2*p4;
}

CVector3 positionBarycentricCoord_Wedge
(double r0,
 double r1,
 double r2,
 const CVector3& p0,
 const CVector3& p1,
 const CVector3& p2,
 const CVector3& p3,
 const CVector3& p4,
 const CVector3& p5)
{
  return (1.0-r2)*r0*p0
  + (1.0-r2)*r1*p1
  + (1.0-r2)*(1.0-r0-r1)*p2
  + r2*r0*p3
  + r2*r1*p4
  + r2*(1.0-r0-r1)*p5;
}

void iteration_barycentricCoord_Origin_Solid
(double& r0,
 double& r1,
 double& r2,
 const CVector3& q, // q=positionBarycentricCoord_Wedge
 const CVector3& dpdr0,
 const CVector3& dpdr1,
 const CVector3& dpdr2,
 double damp=1.0)
{
  const double cxx = dpdr0.x*dpdr0.x + dpdr1.x*dpdr1.x + dpdr2.x*dpdr2.x;
  const double cxy = dpdr0.x*dpdr0.y + dpdr1.x*dpdr1.y + dpdr2.x*dpdr2.y;
  const double cxz = dpdr0.x*dpdr0.z + dpdr1.x*dpdr1.z + dpdr2.x*dpdr2.z;
  const double cyy = dpdr0.y*dpdr0.y + dpdr1.y*dpdr1.y + dpdr2.y*dpdr2.y;
  const double cyz = dpdr0.y*dpdr0.z + dpdr1.y*dpdr1.z + dpdr2.y*dpdr2.z;
  const double czz = dpdr0.z*dpdr0.z + dpdr1.z*dpdr1.z + dpdr2.z*dpdr2.z;
  double C[9] = {cxx,cxy,cxz, cxy,cyy,cyz, cxz,cyz,czz};
  double Cinv[9]; InverseMat3(Cinv, C);
  const CVector3 d = damp*MatVec(Cinv,q);
  r0 -= dpdr0*d;
  r1 -= dpdr1*d;
  r2 -= dpdr2*d;
}

bool barycentricCoord_Origin_Tet
(double& r0,
 double& r1,
 double& r2,
 const CVector3& p0,
 const CVector3& p1,
 const CVector3& p2,
 const CVector3& p3)
{
  double v0 = volume_OrgTet(p1, p2, p3);
  double v1 = volume_OrgTet(p2, p0, p3);
  double v2 = volume_OrgTet(p1, p3, p0);
  double v3 = volume_OrgTet(p1, p0, p2);
  double vt_inv = 1.0/(v0+v1+v2+v3);
  r0 = v0*vt_inv;
  r1 = v1*vt_inv;
  r2 = v2*vt_inv;
  return true;
}

bool barycentricCoord_Origin_Pyramid
(double& r0,
 double& r1,
 double& r2,
 const CVector3& p0,
 const CVector3& p1,
 const CVector3& p2,
 const CVector3& p3,
 const CVector3& p4)
{
  CVector3 q = positionBarycentricCoord_Pyramid(r0,r1,r2, p0,p1,p2,p3,p4);
  for(int itr=0;itr<5;++itr){
    const CVector3 dpdr0 = -(1.0-r2)*(1.0-r1)*p0 + (1.0-r2)*(1.0-r1)*p1 + (1.0-r2)*r1*p2 - (1.0-r2)*r1*p3;
    const CVector3 dpdr1 = -(1.0-r2)*(1.0-r0)*p0 - (1.0-r2)*r0*p1 + (1.0-r2)*r0*p2 + (1.0-r2)*(1.0-r0)*p3;
    const CVector3 dpdr2 = -(1.0-r0)*(1.0-r1)*p0 - r0*(1.0-r1)*p1 - r0*r1*p2 - (1.0-r0)*r1*p3 + p4;
    iteration_barycentricCoord_Origin_Solid(r0,r1,r2, q,
                                            dpdr0,dpdr1,dpdr2,
                                            1.0);
    q = positionBarycentricCoord_Pyramid(r0,r1,r2, p0,p1,p2,p3,p4);
  }
  return true;
}

bool barycentricCoord_Origin_Wedge
(double& r0,
 double& r1,
 double& r2,
 const CVector3& p0,
 const CVector3& p1,
 const CVector3& p2,
 const CVector3& p3,
 const CVector3& p4,
 const CVector3& p5)
{
  CVector3 q = positionBarycentricCoord_Wedge(r0,r1,r2, p0,p1,p2,p3,p4,p5);
  for(int itr=0;itr<5;++itr){
    const CVector3 dpdr0 = (1.0-r2)*(p0-p2)+r2*(p3-p5);
    const CVector3 dpdr1 = (1.0-r2)*(p1-p2)+r2*(p4-p5);
    const CVector3 dpdr2 = r0*(p3-p0)+r1*(p4-p1)+(1.0-r0-r1)*(p5-p2);
    iteration_barycentricCoord_Origin_Solid(r0,r1,r2, q,
                                            dpdr0,dpdr1,dpdr2,
                                            1.0);
    q = positionBarycentricCoord_Wedge(r0,r1,r2, p0,p1,p2,p3,p4,p5);
  }
  return true;
}


bool IsInside_Orgin_BoundingBoxPoint6
(const CVector3& p0,
 const CVector3& p1,
 const CVector3& p2,
 const CVector3& p3,
 const CVector3& p4,
 const CVector3& p5)
{
  if( p0.x>0 && p1.x>0 && p2.x>0 && p3.x>0 && p4.x>0 && p5.x>0 ){ return false; }
  if( p0.x<0 && p1.x<0 && p2.x<0 && p3.x<0 && p4.x<0 && p5.x<0 ){ return false; }
  if( p0.y>0 && p1.y>0 && p2.y>0 && p3.y>0 && p4.y>0 && p5.y>0 ){ return false; }
  if( p0.y<0 && p1.y<0 && p2.y<0 && p3.y<0 && p4.y<0 && p5.y<0 ){ return false; }
  if( p0.z>0 && p1.z>0 && p2.z>0 && p3.z>0 && p4.z>0 && p5.z>0 ){ return false; }
  if( p0.z<0 && p1.z<0 && p2.z<0 && p3.z<0 && p4.z<0 && p5.z<0 ){ return false; }
  return true;
}

bool IsInside_Orgin_BoundingBoxPoint5
(const CVector3& p0,
 const CVector3& p1,
 const CVector3& p2,
 const CVector3& p3,
 const CVector3& p4)
{
  if( p0.x>0 && p1.x>0 && p2.x>0 && p3.x>0 && p4.x>0 ){ return false; }
  if( p0.x<0 && p1.x<0 && p2.x<0 && p3.x<0 && p4.x<0 ){ return false; }
  if( p0.y>0 && p1.y>0 && p2.y>0 && p3.y>0 && p4.y>0 ){ return false; }
  if( p0.y<0 && p1.y<0 && p2.y<0 && p3.y<0 && p4.y<0 ){ return false; }
  if( p0.z>0 && p1.z>0 && p2.z>0 && p3.z>0 && p4.z>0 ){ return false; }
  if( p0.z<0 && p1.z<0 && p2.z<0 && p3.z<0 && p4.z<0 ){ return false; }
  return true;
}

bool IsInside_Orgin_BoundingBoxPoint4
(const CVector3& p0,
 const CVector3& p1,
 const CVector3& p2,
 const CVector3& p3)
{
  if( p0.x>0 && p1.x>0 && p2.x>0 && p3.x>0 ){ return false; }
  if( p0.x<0 && p1.x<0 && p2.x<0 && p3.x<0 ){ return false; }
  if( p0.y>0 && p1.y>0 && p2.y>0 && p3.y>0 ){ return false; }
  if( p0.y<0 && p1.y<0 && p2.y<0 && p3.y<0 ){ return false; }
  if( p0.z>0 && p1.z>0 && p2.z>0 && p3.z>0 ){ return false; }
  if( p0.z<0 && p1.z<0 && p2.z<0 && p3.z<0 ){ return false; }
  return true;
}




////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


// matrix are column major
static CVector3 mult_GlAffineMatrix
(const float* m,
 const CVector3& p)
{
  CVector3 v;
  v.x = m[0*4+0]*p.x + m[1*4+0]*p.y + m[2*4+0]*p.z + m[3*4+0];
  v.y = m[0*4+1]*p.x + m[1*4+1]*p.y + m[2*4+1]*p.z + m[3*4+1];
  v.z = m[0*4+2]*p.x + m[1*4+2]*p.y + m[2*4+2]*p.z + m[3*4+2];
  return v;
}

CVector3 solve_GlAffineMatrix(const float* m,
                              const CVector3& p)
{
  CVector3 v = p - CVector3(m[3*4+0],m[3*4+1],m[3*4+2]);
  double M[9] = {
    m[0*4+0],m[1*4+0],m[2*4+0],
    m[0*4+1],m[1*4+1],m[2*4+1],
    m[0*4+2],m[1*4+2],m[2*4+2] };
  double Minv[9];  InverseMat3(Minv, M);
  return MatVec(Minv,v);
//  CMatrix3 Minv = M.Inverse();  
//  return Minv*v;
}

CVector3 solve_GlAffineMatrixDirection(const float* m,
                                       const CVector3& v)
{
  double M[9] = {
    m[0*4+0],m[1*4+0],m[2*4+0],
    m[0*4+1],m[1*4+1],m[2*4+1],
    m[0*4+2],m[1*4+2],m[2*4+2] };
  double Minv[9];  InverseMat3(Minv, M);
  return MatVec(Minv,v);
  /*
  CMatrix3 M(m[0*4+0],m[1*4+0],m[2*4+0],
             m[0*4+1],m[1*4+1],m[2*4+1],
             m[0*4+2],m[1*4+2],m[2*4+2]);
   */
  /*
   CMatrix3 M(m[0*4+0], m[0*4+1], m[0*4+2],
   m[1*4+0], m[1*4+1], m[1*4+2],
   m[2*4+0], m[2*4+1], m[2*4+2]);
   */
//  CMatrix3 Minv = M.Inverse();
//  return Minv*v;
}

//////////////////////

CVector3 screenProjection(const CVector3& v,
                          const float* mMV, const float* mPj)
{
  CVector3 v0 = mult_GlAffineMatrix(mMV, v );
  CVector3 v1 = mult_GlAffineMatrix(mPj, v0);
  float w1 = mPj[11]*v0.z + mPj[15];
  return CVector3(v1.x/w1, v1.y/w1, 0.0);
}

CVector3 screenUnProjection(const CVector3& v,
                            const float* mMV, const float* mPj)
{
  float D = mPj[11] + mPj[15]; // z is 1 after model view
  CVector3 v0( D*v.x, D*v.y, 0.0 );
  CVector3 v1 = solve_GlAffineMatrix(mPj, v0);
  v1.z = 1;
  CVector3 v2 = solve_GlAffineMatrix(mMV, v1);
  return v2;
}

CVector3 screenUnProjectionDirection(const CVector3& v,
                                     const float* mMV, const float* mPj)
{
  CVector3 v0 = solve_GlAffineMatrixDirection(mPj, v);
  CVector3 v1 = solve_GlAffineMatrixDirection(mMV, v0);
  v1.SetNormalizedVector();
  return v1;
}

CVector3 screenDepthDirection
(const CVector3& v,
 const float* mMV,
 const float* mPj)
{
  float Dv = mPj[11] + mPj[15]; // z is 1 after model view
  CVector3 v0( Dv*v.x, Dv*v.y, 0.0 );
  CVector3 v1 = solve_GlAffineMatrix(mPj, v0);
  v1.z = 1;
  ////
  float Du = mPj[11]*2.0 + mPj[15]; // z is 1 after model view
  CVector3 u0( Du*v.x, Du*v.y, 0.0 );
  CVector3 u1 = solve_GlAffineMatrix(mPj, u0);
  u1.z = 2;
  ////
  CVector3 v2 = solve_GlAffineMatrixDirection(mMV, (v1-u1) );
  v2.SetNormalizedVector();
  return v2;
}


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


//! Volume of a tetrahedra
double volume_Tet
(const CVector3& v0,
 const CVector3& v1,
 const CVector3& v2,
 const CVector3& v3 )
{
  double v = (v1.x-v0.x)*( (v2.y-v0.y)*(v3.z-v0.z) - (v3.y-v0.y)*(v2.z-v0.z) )
           + (v1.y-v0.y)*( (v2.z-v0.z)*(v3.x-v0.x) - (v3.z-v0.z)*(v2.x-v0.x) )
           + (v1.z-v0.z)*( (v2.x-v0.x)*(v3.y-v0.y) - (v3.x-v0.x)*(v2.y-v0.y) );
  return v*0.16666666666666666666666666666667;
}

//! Volume of a tetrahedra v0 is orgin
double volume_OrgTet
(const CVector3& v1,
 const CVector3& v2,
 const CVector3& v3 )
{
  double v = v1.x*(v2.y*v3.z-v3.y*v2.z)
           + v1.y*(v2.z*v3.x-v3.z*v2.x)
           + v1.z*(v2.x*v3.y-v3.x*v2.y);
  return v*0.16666666666666666666666666666667;
}

double volume_Pyramid
(const CVector3& p0,
 const CVector3& p1,
 const CVector3& p2,
 const CVector3& p3,
 const CVector3& p4)
{
  double v0124 = volume_Tet(p0, p1, p2, p4);
  double v0234 = volume_Tet(p0, p2, p3, p4);
  double v0134 = volume_Tet(p0, p1, p3, p4);
  double v2314 = volume_Tet(p2, p3, p1, p4);
  double v0 = v0124+v0234;
  double v1 = v0134+v2314;
  return (v0+v1)*0.5;
}

double volume_Wedge
(const CVector3& p0,
 const CVector3& p1,
 const CVector3& p2,
 const CVector3& p3,
 const CVector3& p4,
 const CVector3& p5)
{
  CVector3 pm = (p0+p1+p2+p3+p4+p5)/6.0;
  double vm012 = volume_Tet(pm,p0,p1,p2);
  double vm435 = volume_Tet(pm,p4,p3,p5);
  double vp0143 = volume_Pyramid(p0,p1,p4,p3,pm);
  double vp1254 = volume_Pyramid(p1,p2,p5,p4,pm);
  double vp2035 = volume_Pyramid(p2,p2,p3,p5,pm);
  return vm012+vm435+vp0143+vp1254+vp2035;
}

////////////////////////////////////////////////////////////////////////////////

double SolidAngleTri
(const CVector3& v1,
const CVector3& v2,
const CVector3& v3)
{
  double l1 = v1.Length();
  double l2 = v2.Length();
  double l3 = v3.Length();
  double den = (v1^v2)*v3;
  double num = l1*l2*l3+(v1*v2)*l3+(v2*v3)*l1+(v3*v1)*l2;
  double tho = den/num;
  double v = atan(tho);
  if (v<0){ v += 2*M_PI; }
  v *= 2;
  return v;
}

////////////////////////////////////////////////

void Cross( CVector3& lhs, const CVector3& v1, const CVector3& v2 ){
  lhs.x = v1.y*v2.z - v2.y*v1.z;
  lhs.y = v1.z*v2.x - v2.z*v1.x;
  lhs.z = v1.x*v2.y - v2.x*v1.y;
}

double TriArea(const CVector3& v1, const CVector3& v2, const CVector3& v3)
{
  double x, y, z;
  x = ( v2.y - v1.y )*( v3.z - v1.z ) - ( v3.y - v1.y )*( v2.z - v1.z );
  y = ( v2.z - v1.z )*( v3.x - v1.x ) - ( v3.z - v1.z )*( v2.x - v1.x );
  z = ( v2.x - v1.x )*( v3.y - v1.y ) - ( v3.x - v1.x )*( v2.y - v1.y );
  return 0.5*sqrt( x*x + y*y + z*z );
}


double SquareTriArea(const CVector3& v1, const CVector3& v2, const CVector3& v3)
{
  double dtmp_x = (v2.y-v1.y)*(v3.z-v1.z)-(v2.z-v1.z)*(v3.y-v1.y);
  double dtmp_y = (v2.z-v1.z)*(v3.x-v1.x)-(v2.x-v1.x)*(v3.z-v1.z);
  double dtmp_z = (v2.x-v1.x)*(v3.y-v1.y)-(v2.y-v1.y)*(v3.x-v1.x);
  return (dtmp_x*dtmp_x + dtmp_y*dtmp_y + dtmp_z*dtmp_z)*0.25;
}

////////////////////////////////////////////////

double SquareDistance(const CVector3& ipo0, const CVector3& ipo1)
{
  return	( ipo1.x - ipo0.x )*( ipo1.x - ipo0.x ) + ( ipo1.y - ipo0.y )*( ipo1.y - ipo0.y ) + ( ipo1.z - ipo0.z )*( ipo1.z - ipo0.z );
}

double SquareLength(const CVector3& point)
{
  return	point.x*point.x + point.y*point.y + point.z*point.z;
}

////////////////////////////////////////////////

//! length of vector
double Length(const CVector3& point)
{
  return	sqrt( point.x*point.x + point.y*point.y + point.z*point.z );
}

//! distance between two points
double Distance(const CVector3& p0, const CVector3& p1)
{
  return	sqrt( SquareDistance(p0,p1) );
}

////////////////////////////////////////////////

double SqareLongestEdgeLength
(const CVector3& ipo0,
const CVector3& ipo1,
const CVector3& ipo2,
const CVector3& ipo3 )
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

////////////////////////////////////////////////

double LongestEdgeLength
(const CVector3& ipo0,
 const CVector3& ipo1,
 const CVector3& ipo2,
 const CVector3& ipo3 )
{
  return sqrt( SqareLongestEdgeLength(ipo0,ipo1,ipo2,ipo3) );
}

////////////////////////////////////////////////

double SqareShortestEdgeLength
(const CVector3& ipo0,
 const CVector3& ipo1,
 const CVector3& ipo2,
 const CVector3& ipo3 )
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

////////////////////////////////////////////////


double ShortestEdgeLength
(const CVector3& ipo0,
 const CVector3& ipo1,
 const CVector3& ipo2,
 const CVector3& ipo3 )
{
  return sqrt( SqareShortestEdgeLength(ipo0,ipo1,ipo2,ipo3) );
}

////////////////////////////////////////////////

void Normal
(CVector3& vnorm,
 const CVector3& v1,
 const CVector3& v2,
 const CVector3& v3)
{
  vnorm.x = (v2.y-v1.y)*(v3.z-v1.z)-(v2.z-v1.z)*(v3.y-v1.y);
  vnorm.y = (v2.z-v1.z)*(v3.x-v1.x)-(v2.x-v1.x)*(v3.z-v1.z);
  vnorm.z = (v2.x-v1.x)*(v3.y-v1.y)-(v2.y-v1.y)*(v3.x-v1.x);
}

CVector3 Normal
(const CVector3& v1,
const CVector3& v2,
const CVector3& v3)
{
  CVector3 vnorm;
  vnorm.x = (v2.y-v1.y)*(v3.z-v1.z)-(v2.z-v1.z)*(v3.y-v1.y);
  vnorm.y = (v2.z-v1.z)*(v3.x-v1.x)-(v2.x-v1.x)*(v3.z-v1.z);
  vnorm.z = (v2.x-v1.x)*(v3.y-v1.y)-(v2.y-v1.y)*(v3.x-v1.x);
  return vnorm;
}


////////////////////////////////////////////////

void UnitNormal
(CVector3& vnorm,
 const CVector3& v1,
 const CVector3& v2,
 const CVector3& v3)
{
  vnorm.x = (v2.y-v1.y)*(v3.z-v1.z)-(v2.z-v1.z)*(v3.y-v1.y);
  vnorm.y = (v2.z-v1.z)*(v3.x-v1.x)-(v2.x-v1.x)*(v3.z-v1.z);
  vnorm.z = (v2.x-v1.x)*(v3.y-v1.y)-(v2.y-v1.y)*(v3.x-v1.x);
  const double dtmp1 = 1.0 / Length(vnorm);
  vnorm.x *= dtmp1;
  vnorm.y *= dtmp1;
  vnorm.z *= dtmp1;
}

CVector3 UnitNormal
(const CVector3& v1,
const CVector3& v2,
const CVector3& v3)
{
  CVector3 vnorm;
  vnorm.x = (v2.y-v1.y)*(v3.z-v1.z)-(v2.z-v1.z)*(v3.y-v1.y);
  vnorm.y = (v2.z-v1.z)*(v3.x-v1.x)-(v2.x-v1.x)*(v3.z-v1.z);
  vnorm.z = (v2.x-v1.x)*(v3.y-v1.y)-(v2.y-v1.y)*(v3.x-v1.x);
  const double dtmp1 = 1.0/Length(vnorm);
  vnorm.x *= dtmp1;
  vnorm.y *= dtmp1;
  vnorm.z *= dtmp1;
  return vnorm;
}

////////////////////////////////////////////////

double SquareCircumradius
(const CVector3& ipo0,
 const CVector3& ipo1,
 const CVector3& ipo2,
 const CVector3& ipo3)
{
  double base[3][3] = {
    { ipo1.x-ipo0.x, ipo1.y-ipo0.y, ipo1.z-ipo0.z },
    { ipo2.x-ipo0.x, ipo2.y-ipo0.y, ipo2.z-ipo0.z },
    { ipo3.x-ipo0.x, ipo3.y-ipo0.y, ipo3.z-ipo0.z }
  };
  double s[6] = {
    base[0][0]*base[0][0]+base[0][1]*base[0][1]+base[0][2]*base[0][2],
    base[1][0]*base[1][0]+base[1][1]*base[1][1]+base[1][2]*base[1][2],
    base[2][0]*base[2][0]+base[2][1]*base[2][1]+base[2][2]*base[2][2],
    base[1][0]*base[2][0]+base[1][1]*base[2][1]+base[1][2]*base[2][2],
    base[2][0]*base[0][0]+base[2][1]*base[0][1]+base[2][2]*base[0][2],
    base[0][0]*base[1][0]+base[0][1]*base[1][1]+base[0][2]*base[1][2],
  };
  const double vol = volume_Tet(ipo0,ipo1,ipo2,ipo3)*6.0;
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
   vec1.x = base[0][0]*u[0]+base[1][0]*u[1]+base[2][0]*u[2] + ipo0.x;
   vec1.y = base[0][1]*u[0]+base[1][1]*u[1]+base[2][1]*u[2] + ipo0.y;
   vec1.z = base[0][2]*u[0]+base[1][2]*u[1]+base[2][2]*u[2] + ipo0.z;
   std::cout << square_radius << " ";
   std::cout << SquareLength(vec1,ipo0) << " ";
   std::cout << SquareLength(vec1,ipo1) << " ";
   std::cout << SquareLength(vec1,ipo2) << " ";
   std::cout << SquareLength(vec1,ipo3) << std::endl;;
   return square_radius;
   */
}

CVector3 CircumCenter
(const CVector3& ipo0,
 const CVector3& ipo1,
 const CVector3& ipo2,
 const CVector3& ipo3)
{
  
  double base[3][3] = {
    { ipo1.x-ipo0.x, ipo1.y-ipo0.y, ipo1.z-ipo0.z },
    { ipo2.x-ipo0.x, ipo2.y-ipo0.y, ipo2.z-ipo0.z },
    { ipo3.x-ipo0.x, ipo3.y-ipo0.y, ipo3.z-ipo0.z }
  };
  double s[6] = {
    base[0][0]*base[0][0]+base[0][1]*base[0][1]+base[0][2]*base[0][2],
    base[1][0]*base[1][0]+base[1][1]*base[1][1]+base[1][2]*base[1][2],
    base[2][0]*base[2][0]+base[2][1]*base[2][1]+base[2][2]*base[2][2],
    base[1][0]*base[2][0]+base[1][1]*base[2][1]+base[1][2]*base[2][2],
    base[2][0]*base[0][0]+base[2][1]*base[0][1]+base[2][2]*base[0][2],
    base[0][0]*base[1][0]+base[0][1]*base[1][1]+base[0][2]*base[1][2],
  };
  const double vol = volume_Tet(ipo0,ipo1,ipo2,ipo3)*6.0;
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
  CVector3 vec1;
  vec1.x = base[0][0]*u[0]+base[1][0]*u[1]+base[2][0]*u[2] + ipo0.x;
  vec1.y = base[0][1]*u[0]+base[1][1]*u[1]+base[2][1]*u[2] + ipo0.y;
  vec1.z = base[0][2]*u[0]+base[1][2]*u[1]+base[2][2]*u[2] + ipo0.z;
  return vec1;
}

void MeanValueCoordinate
(double w[3],
 const CVector3& v0,
 const CVector3& v1,
 const CVector3& v2)
{
  double eps  = 1.0e-5;
  double d0 = v0.Length();
  double d1 = v1.Length();
  double d2 = v2.Length();
  const CVector3 u0 = v0/d0;
  const CVector3 u1 = v1/d1;
  const CVector3 u2 = v2/d2;
  double l0 = (u1-u2).Length();
  double l1 = (u2-u0).Length();
  double l2 = (u0-u1).Length();
  if( l0<eps || l1<eps || l2<eps ){
    w[0] = 0;
    w[1] = 0;
    w[2] = 0;
    return;
  }
  double t0 = 2*asin(l0*0.5);
  double t1 = 2*asin(l1*0.5);
  double t2 = 2*asin(l2*0.5);
  double h = (t0+t1+t2)*0.5;
  double c0 = 2*sin(h)*sin(h-t0)/(sin(t1)*sin(t2))-1;
  double c1 = 2*sin(h)*sin(h-t1)/(sin(t2)*sin(t0))-1;
  double c2 = 2*sin(h)*sin(h-t2)/(sin(t0)*sin(t1))-1;
  double vol012 = ScalarTripleProduct(u0,u1,u2);
  double sign = (vol012 > 0) ? 1 : -1;
  double s0 = sign*sqrt(1.0-c0*c0);
  double s1 = sign*sqrt(1.0-c1*c1);
  double s2 = sign*sqrt(1.0-c2*c2);
  if( isnan_vector3(s0) || isnan_vector3(s1) || isnan_vector3(s2) ){
    w[0] = 0;
    w[1] = 0;
    w[2] = 0;
    return;
  }
  if( fabs(d0*sin(t1)*s2)<eps || fabs(d1*sin(t2)*s0)<eps || fabs(d2*sin(t0)*s1)<eps ){
    w[0] = 0;
    w[1] = 0;
    w[2] = 0;
    return;
  }
  w[0] = (t0-c2*t1-c1*t2)/(d0*sin(t1)*s2);
  w[1] = (t1-c0*t2-c2*t0)/(d1*sin(t2)*s0);
  w[2] = (t2-c1*t0-c0*t1)/(d2*sin(t0)*s1);
}



////////////////////////////////////////////////



CVector3 ProjectPointOnTriangle
(const CVector3 &p0,
 const CVector3 &tri_p1, const CVector3 &tri_p2, const CVector3 &tri_p3)
{
  CVector3 normal = Cross(tri_p2 - tri_p1, tri_p3 - tri_p1);
  double cosAlpha = Dot(p0 - tri_p1, normal) / (Length(p0 - tri_p1) * Length(normal));
  double lenP0ProjectedP0 = Length(tri_p1 - p0) * cosAlpha;
  CVector3 p0ProjectedP0 = -1 * lenP0ProjectedP0 * normal / Length(normal);
  
  return p0 + p0ProjectedP0;
}

bool isRayIntersectingTriangle
(const CVector3 &line0, const CVector3 &line1,
 const CVector3 &tri0, const CVector3 &tri1, const CVector3 &tri2,
 CVector3 &intersectionPoint)
{
  CVector3 normal = Cross(tri1 - tri0, tri2 - tri0);
  
  // The ray is parallel to the triangle plane
  if (Dot(normal, line1 - line0) == 0)
  {
    return false;
  }
  
  double r = Dot(normal, tri0 - line0) / Dot(normal, line1 - line0);
  
  // The ray does not intersect the triangle plane
  if (r < 0)
  {
    return false;
  }
  
  // Find the intersection point
  intersectionPoint = line0 + r * (line1 - line0);
  
  if (!isPointInsideTriangle(intersectionPoint,
                             tri0, tri1, tri2))
  {
    return false;
  }
  
  return true;
}

bool isPointInsideTriangle
(const CVector3 &p0,
 const CVector3 &tri_p1, const CVector3 &tri_p2, const CVector3 &tri_p3)
{
  if (isPointSameSide(p0, tri_p1, tri_p2, tri_p3)
      && isPointSameSide(p0, tri_p2, tri_p1, tri_p3)
      && isPointSameSide(p0, tri_p3, tri_p1, tri_p2))
  {
    return true;
  } else {
    return false;
  }
}

bool isPointSameSide
(const CVector3 &p0, const CVector3 &p1,
 const CVector3 &line_p0, const CVector3 &line_p1)
{
  CVector3 crossProd1 = Cross(line_p1 - line_p0, p0 - line_p0);
  CVector3 crossProd2 = Cross(line_p1 - line_p0, p1 - line_p0);
  
  if (Dot(crossProd1, crossProd2) >= 0)
  {
    return true;
  } else {
    return false;
  }
}

////////////////////////////////

//! check if Delaunay condition satisfied
// 0 : p3 is inside circum circle on the p0,p1,p2
// 1 :       on         
// 2 :       outsdie 
int DetDelaunay
(const CVector3& p0,
const CVector3& p1,
const CVector3& p2,
const CVector3& p3)
{
  const double area = TriArea(p0, p1, p2);
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

  const CVector3 out_center = etmp0*p0+etmp1*p1+etmp2*p2;

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

/*!
 curcumradius of a tetrahedra
 */
double Circumradius
(const CVector3& ipo0,
 const CVector3& ipo1,
 const CVector3& ipo2,
 const CVector3& ipo3)
{
  return sqrt( SquareCircumradius(ipo0,ipo1,ipo2,ipo3) );
}


CVector3 RotateVector(const CVector3& vec0, const CVector3& rot )
{
  const double theta = rot.Length();
  if( theta < 1.0e-30 ){
    return vec0;
  }
  CVector3 e0 = rot;
  e0.SetNormalizedVector();
  CVector3 e2 = ::Cross(e0,vec0);
  if( e2.Length() < 1.0e-30 ){
    return vec0;
  }
  e2.SetNormalizedVector();
  CVector3 e1 = ::Cross(e2,e0);
  assert( fabs( e1.Length() - 1 ) < 1.0e-10 );
  //	assert( e2.x*vec_0.x + e2.y*vec_0.y + e2.z*vec_0.z < 1.0e-10 );
  const double dot00 = Dot(vec0,e0);
  const double dot01 = Dot(vec0,e1);
  const double cost = cos(theta);
  const double sint = sin(theta);
  CVector3 vec1;
  vec1.x = dot00*e0.x + dot01*cost*e1.x + dot01*sint*e2.x;
  vec1.y = dot00*e0.y + dot01*cost*e1.y + dot01*sint*e2.y;
  vec1.z = dot00*e0.z + dot01*cost*e1.z + dot01*sint*e2.z;
  return vec1;
}

CVector3 RandVector(){
  CVector3 r;
  r.x = (2*(double)rand()/(RAND_MAX+1.0)-1);
  r.y = (2*(double)rand()/(RAND_MAX+1.0)-1);
  r.z = (2*(double)rand()/(RAND_MAX+1.0)-1);
  return r;
}

CVector3 RandUnitVector(){
  for(int itr=0;itr<100;itr++){
    CVector3 r = RandVector();
    double l = r.Length();
    if( (l <= 1 || itr==9) && l > 1.0e-5 ){
      r.SetNormalizedVector();
      return r;
    }
  }
  return CVector3(1,0,0);
}

CVector3 RandGaussVector()
{
  double a0 = rand()/(RAND_MAX+1.0);
  double a1 = rand()/(RAND_MAX+1.0);
  double a2 = rand()/(RAND_MAX+1.0);
  double a3 = rand()/(RAND_MAX+1.0);
  
  double x = sqrt(-2.0*log(a0))*cos(3.1415*2*a1);
  double y = sqrt(-2.0*log(a0))*sin(3.1415*2*a1);
  double z = sqrt(-2.0*log(a2))*cos(3.1415*2*a3);
  return CVector3(x,y,z);
}



////////////////////////////////////////////////
// using <vector> from here





////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


void GetConstConstDiff_Bend
(double& C, CVector3 dC[4],
 const CVector3& p0,
 const CVector3& p1,
 const CVector3& p2,
 const CVector3& p3)
{
  const CVector3 v02 = p2-p0;
  const CVector3 v03 = p3-p0;
  const CVector3 v12 = p2-p1;
  const CVector3 v13 = p3-p1;
  const CVector3 v23 = p3-p2;
  ////
  const CVector3 A = v02^v03;
  const CVector3 B = v13^v12;
  const double lA = A.Length();
  const double lB = B.Length();
  const CVector3 a = A/lA;
  const CVector3 b = B/lB;
  const double ab = a*b;
  //  C = acos(ab);
  C = ab-1;
  const double sab = 1.0;//-1.0/sin(C);
  const CVector3 tmpBA = (b-a*(a*b))*(sab/lA);
  const CVector3 tmpAB = (a-b*(b*a))*(sab/lB);
  dC[0] = (tmpBA^v23);
  dC[1] = (v23^tmpAB);
  dC[2] = (v03^tmpBA) + (tmpAB^v13);
  dC[3] = (tmpBA^v02) + (v12^tmpAB);
}
/*
 {
 CVector3& p0 = aXYZ1[ip0];
 CVector3& p1 = aXYZ1[ip1];
 const double len01_ini = (aXYZ_ini[ip1]-aXYZ_ini[ip0]).length();
 const double len01 = (p1-p0).length();
 double w0 = aInvMass[ip0];
 double w1 = aInvMass[ip1];
 const float4& n01 = (p1-p0)/len01;
 p0 += +w0/(w0+w1)*(len01-len01_ini)*n01;
 p1 += -w1/(w0+w1)*(len01-len01_ini)*n01;
 }
 */

void CheckConstDiff_Bend()
{
  CVector3 p[4];
  for(int ino=0;ino<4;++ino){
    p[ino].x = (double)rand()/(RAND_MAX+1.0);
    p[ino].y = (double)rand()/(RAND_MAX+1.0);
    p[ino].z = (double)rand()/(RAND_MAX+1.0);
  }
  double C; CVector3 dC[4];
  GetConstConstDiff_Bend(C, dC, p[0],p[1],p[2],p[3]);
  for(int ino=0;ino<4;++ino){
    for(int idim=0;idim<3;++idim){
      CVector3 p1[4] = {p[0],p[1],p[2],p[3]};
      double eps = 1.0e-4;
      p1[ino][idim] += eps;
      double C1; CVector3 dC1[4];
      GetConstConstDiff_Bend(C1, dC1, p1[0],p1[1],p1[2],p1[3]);
      std::cout << "   " << ino << " " << idim << "   -->  " << (C1-C)/eps << " " << dC[ino][idim] << std::endl;
    }
  }
}


////////////////////////////////////////////////
// TODO: following should be move to mesh class?

double TriArea
(const int iv1, const int iv2, const int iv3,
 const std::vector<CVector3>& aPoint )
{
  return TriArea(aPoint[iv1],aPoint[iv2],aPoint[iv3]);
}

double volume_Tet
( int iv1, int iv2, int iv3, int iv4,
 const std::vector<CVector3>& aPoint)
{
  return volume_Tet(aPoint[iv1],aPoint[iv2],aPoint[iv3],aPoint[iv4]);
}

///////

bool IsOut
(int itri, const CVector3& v,
 const std::vector<CVector3>& aXYZ,
 const std::vector<int>& aTri)
{
  int i0 = aTri[itri*3+0];
  int i1 = aTri[itri*3+1];
  int i2 = aTri[itri*3+2];
  const CVector3& v0 = aXYZ[i0];
  const CVector3& v1 = aXYZ[i1];
  const CVector3& v2 = aXYZ[i2];
  CVector3 n; Normal(n,v0,v1,v2);
  double dot = Dot(v-v0,n);
  return dot > 0;
}

void ConvexHull
(std::vector<int>& aTri,
 const std::vector<CVector3>& aXYZ)
{
  std::vector<int> aBflg( aXYZ.size(), -1 );
  aTri.reserve(aXYZ.size()*6);
  aTri.resize(4*3);
  aTri[ 0] = 1;  aTri[ 1] = 2;  aTri[ 2] = 3; // 0
  aTri[ 3] = 0;  aTri[ 4] = 3;  aTri[ 5] = 2; // 1
  aTri[ 6] = 0;  aTri[ 7] = 1;  aTri[ 8] = 3; // 2
  aTri[ 9] = 0;  aTri[10] = 2;  aTri[11] = 1; // 3
  std::vector< std::pair<int,int> > aTriSur;
  aTriSur.resize(4*3);
  aTriSur[ 0] = std::make_pair(1,0);
  aTriSur[ 1] = std::make_pair(2,0);
  aTriSur[ 2] = std::make_pair(3,0);
  aTriSur[ 3] = std::make_pair(0,0);
  aTriSur[ 4] = std::make_pair(3,2);
  aTriSur[ 5] = std::make_pair(2,1);
  aTriSur[ 6] = std::make_pair(0,1);
  aTriSur[ 7] = std::make_pair(1,2);
  aTriSur[ 8] = std::make_pair(3,1);
  aTriSur[ 9] = std::make_pair(0,2);
  aTriSur[10] = std::make_pair(2,2);
  aTriSur[11] = std::make_pair(1,1);
  {
    double vol = ::volume_Tet(0, 1, 2, 3, aXYZ);
    if( vol < 0 ){
      aTri[ 0] = 3;  aTri[ 1] = 2;  aTri[ 2] = 1; // 0
      aTri[ 3] = 2;  aTri[ 4] = 3;  aTri[ 5] = 0; // 1
      aTri[ 6] = 3;  aTri[ 7] = 1;  aTri[ 8] = 0; // 2
      aTri[ 9] = 1;  aTri[10] = 2;  aTri[11] = 0; // 3
      aTriSur[ 0] = std::make_pair(3,2);
      aTriSur[ 1] = std::make_pair(2,2);
      aTriSur[ 2] = std::make_pair(1,2);
      aTriSur[ 3] = std::make_pair(2,1);
      aTriSur[ 4] = std::make_pair(3,0);
      aTriSur[ 5] = std::make_pair(0,2);
      aTriSur[ 6] = std::make_pair(3,1);
      aTriSur[ 7] = std::make_pair(1,0);
      aTriSur[ 8] = std::make_pair(0,1);
      aTriSur[ 9] = std::make_pair(1,1);
      aTriSur[10] = std::make_pair(2,0);
      aTriSur[11] = std::make_pair(0,0);
    }
  }
  const int triEd[3][2] = {
    { 1, 2 },
    { 2, 0 },
    { 0, 1 } };
  for(int iv=4;iv<(int)aXYZ.size();iv++){
    CVector3 v = aXYZ[iv];
    int itri_ker = -1;
    for(int itri=0;itri<(int)aTri.size()/3;itri++){
      if( IsOut(itri,v,aXYZ,aTri) ){ itri_ker = itri; break; }
    }
#ifdef DEBUG
    {
      for(int itri=0;itri<aTri.size()/3;itri++){
        for(int ied=0;ied<3;ied++){
          int ied1 = triEd[ied][0];
          int ied2 = triEd[ied][1];
          int itri_s = aTriSur[itri*3+ied].first;
          int ied_s0 = aTriSur[itri*3+ied].second;
          assert( aTriSur[itri_s*3+ied_s0].first  == itri );
          assert( aTriSur[itri_s*3+ied_s0].second == ied );
          int ied_s1 = triEd[ied_s0][0];
          int ied_s2 = triEd[ied_s0][1];
          assert( aTri[itri*3+ied1] == aTri[itri_s*3+ied_s2] );
          assert( aTri[itri*3+ied2] == aTri[itri_s*3+ied_s1] );
        }
      }
    }
#endif
    if( itri_ker == -1 ) continue; // inside
    std::vector< std::pair<int,int> > aB;
    std::vector<int> isDelTri( aTri.size()/3, -1 );
    {
      std::vector<int> isLookedEdge( aTri.size(), -1 );
      std::stack< std::pair<int,int> > sBound;
      { // initialize
        sBound.push( aTriSur[itri_ker*3+0] );
        sBound.push( aTriSur[itri_ker*3+1] );
        sBound.push( aTriSur[itri_ker*3+2] );
        isDelTri[itri_ker] = 1;
      }
      for(;;){
        if( sBound.empty() ) break;
        int itri0 = sBound.top().first;
        int ied0  = sBound.top().second;
        sBound.pop();
        if( isLookedEdge[itri0*3+ied0] == 1 ) continue;
        isLookedEdge[itri0*3+ied0] = 1;
        {
          const std::pair<int,int>& s0 = aTriSur[itri0*3+ied0];
          isLookedEdge[s0.first*3+s0.second] = 1;
        }
        isDelTri[itri0] = ( IsOut(itri0,v,aXYZ,aTri) ) ? 1 : 0;
        if( isDelTri[itri0] == 1 ){ // expand this boundary
          int ied1 = triEd[ied0][0];
          int ied2 = triEd[ied0][1];
          const std::pair<int,int>& s1 = aTriSur[itri0*3+ied1];
          const std::pair<int,int>& s2 = aTriSur[itri0*3+ied2];
          assert( aTriSur[s1.first*3+s1.second].first == itri0 );
          assert( aTriSur[s2.first*3+s2.second].first == itri0 );
          sBound.push( s1 );
          sBound.push( s2 );
        }
        else{ // this is a actuall boundary
          aB.push_back( std::make_pair(itri0,ied0) );
        }
      }
    }
    std::vector<int> aBSur( aB.size()*2, -1);
    {
      for(int ib=0;ib<(int)aB.size();ib++){
        int itri0 = aB[ib].first;
        int itn0  = aB[ib].second;
        int itn1 = triEd[itn0][0];
        int itn2 = triEd[itn0][1];
        int iv1 = aTri[itri0*3+itn1];
        int iv2 = aTri[itri0*3+itn2];
        aBflg[iv1] = -1;
        aBflg[iv2] = -1;
      }
      for(int ib=0;ib<(int)aB.size();ib++){
        int itri0 = aB[ib].first;
        int itn0  = aB[ib].second;
        int itn1 = triEd[itn0][0];
        int itn2 = triEd[itn0][1];
        int iv1 = aTri[itri0*3+itn1];
        int iv2 = aTri[itri0*3+itn2];
        if(      aBflg[iv1] == -2 ){}
        else if( aBflg[iv1] == -1 ){ aBflg[iv1] = ib; }
        else{
          assert( aBflg[iv1] >= 0 );
          int ib0 = aBflg[iv1];
          aBSur[ib *2+1] = ib0;
          aBSur[ib0*2+0] = ib;
          aBflg[iv1] = -2;
        }
        if(      aBflg[iv2] == -2 ){}
        else if( aBflg[iv2] == -1 ){ aBflg[iv2] = ib; }
        else{
          assert( aBflg[iv2] >= 0 );
          int ib0 = aBflg[iv2];
          aBSur[ib *2+0] = ib0;
          aBSur[ib0*2+1] = ib;
          aBflg[iv2] = -2;
        }
      }
    }
#ifdef DEBUG
    for(int ib=0;ib<aB.size();ib++){
      for(int inb=0;inb<2;inb++){
        int itri0 = aB[ib].first;
        int itn0  = aB[ib].second;
        int iv1 = aTri[itri0*3+triEd[itn0][0]];
        int iv2 = aTri[itri0*3+triEd[itn0][1]];
        int ib_s0 = aBSur[ib*2+inb];
        int itri_s0 = aB[ib_s0].first;
        int itn_s0  = aB[ib_s0].second;
        int iv_s1 = aTri[itri_s0*3+triEd[itn_s0][0]];
        int iv_s2 = aTri[itri_s0*3+triEd[itn_s0][1]];
        if( inb == 0 ){ assert( iv2 == iv_s1 ); }
        else{           assert( iv1 == iv_s2 ); }
      }
    }
#endif
    std::vector<int> mapOld2New( aTri.size()/3, -1 );
    std::vector<int> aTri1; aTri1.reserve(aTri.size()+60);
    std::vector< std::pair<int,int> > aTriSur1; aTriSur1.reserve(aTriSur.size()+60);
    for(int itri=0;itri<(int)aTri.size()/3;itri++){
      if( isDelTri[itri] ==  1) continue;
      assert( !IsOut(itri,v,aXYZ,aTri) );
      // itri is inside
      mapOld2New[itri] = (int)aTri1.size()/3;
      aTri1.push_back( aTri[itri*3+0] );
      aTri1.push_back( aTri[itri*3+1] );
      aTri1.push_back( aTri[itri*3+2] );
      aTriSur1.push_back( std::make_pair(-1,0) );
      aTriSur1.push_back( std::make_pair(-1,0) );
      aTriSur1.push_back( std::make_pair(-1,0) );
    }
    for(int itri=0;itri<(int)aTri.size()/3;itri++){ // set old relation
      if( isDelTri[itri] ==  1) continue;
      int jtri0 = mapOld2New[itri];
      assert( jtri0 >= 0 && (int)aTri1.size()/3 );
      for(int iet=0;iet<3;iet++){
        int itri_s = aTriSur[itri*3+iet].first;
        if( mapOld2New[itri_s] == -1 ) continue;
        aTriSur1[jtri0*3+iet].first = mapOld2New[itri_s];
        aTriSur1[jtri0*3+iet].second = aTriSur[itri*3+iet].second;
      }
    }
    const int ntri_old = (int)aTri1.size()/3;
    for(int ib=0;ib<(int)aB.size();ib++){
      int itri0 = aB[ib].first;
      int itn0  = aB[ib].second;
      int itn1 = triEd[itn0][0];
      int itn2 = triEd[itn0][1];
      assert( !IsOut(itri0,v,aXYZ,aTri) );
#ifdef DEBUG
      {
        int itri_s = aTriSur[itri0*3+itn0].first;
        assert( IsOut(itri_s,v,aXYZ,aTri) );
        int ied_s0 = aTriSur[itri0*3+itn0].second;
        assert( aTriSur[itri_s*3+ied_s0].first == itri0 );
        assert( aTriSur[itri_s*3+ied_s0].second == itn0 );
        int ied_s1 = triEd[ied_s0][0];
        int ied_s2 = triEd[ied_s0][1];
        assert( aTri[itri0*3+itn1] == aTri[itri_s*3+ied_s2] );
        assert( aTri[itri0*3+itn2] == aTri[itri_s*3+ied_s1] );
      }
#endif
      assert( isDelTri[itri0] == 0 );
      int jtri0 = mapOld2New[itri0]; assert( jtri0 != -1 );
      int jtri1 = (int)aTri1.size()/3;
      assert( jtri1 == ntri_old + ib );
      aTri1.push_back( iv );
      aTri1.push_back( aTri[itri0*3+itn2] );
      aTri1.push_back( aTri[itri0*3+itn1] );
      aTriSur1[jtri0*3+itn0] = std::make_pair(jtri1,0);
      ////
      int jtri2 = aBSur[ib*2+0] + ntri_old;
      int jtri3 = aBSur[ib*2+1] + ntri_old;
      aTriSur1.push_back( std::make_pair(jtri0,itn0) );
      aTriSur1.push_back( std::make_pair(jtri3,2) );
      aTriSur1.push_back( std::make_pair(jtri2,1) );
    }
    aTri    = aTri1;
    aTriSur = aTriSur1;
  }
}


