/*
 * Copyright (c) 2019 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */


#include <cstdlib>
#include <cstring>
#include <fstream>
#include <sstream>
#include <iostream>
#include <cassert>
#include <cmath>
#include <vector>
#include <stack>
#include <set>

#include "delfem2/msh.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

static double Length3D(const double p[3]){
  return sqrt(p[0]*p[0] + p[1]*p[1] + p[2]*p[2]);
}

static void Cross3D(double r[3], const double v1[3], const double v2[3]){
  r[0] = v1[1]*v2[2] - v2[1]*v1[2];
  r[1] = v1[2]*v2[0] - v2[2]*v1[0];
  r[2] = v1[0]*v2[1] - v2[0]*v1[1];
}

static double TriArea2D(const double p0[], const double p1[], const double p2[]){
  return 0.5*((p1[0]-p0[0])*(p2[1]-p0[1])-(p2[0]-p0[0])*(p1[1]-p0[1]));
}

static double TriArea3D(const double v1[3], const double v2[3], const double v3[3]){
  double n[3];
  n[0] = ( v2[1] - v1[1] )*( v3[2] - v1[2] ) - ( v3[1] - v1[1] )*( v2[2] - v1[2] );
  n[1] = ( v2[2] - v1[2] )*( v3[0] - v1[0] ) - ( v3[2] - v1[2] )*( v2[0] - v1[0] );
  n[2] = ( v2[0] - v1[0] )*( v3[1] - v1[1] ) - ( v3[0] - v1[0] )*( v2[1] - v1[1] );
  return sqrt(n[0]*n[0]+n[1]*n[1]+n[2]*n[2])*0.5;
}

static void UnitNormalAreaTri3D(double n[3], double& a, const double v1[3], const double v2[3], const double v3[3])
{
  n[0] = ( v2[1] - v1[1] )*( v3[2] - v1[2] ) - ( v3[1] - v1[1] )*( v2[2] - v1[2] );
  n[1] = ( v2[2] - v1[2] )*( v3[0] - v1[0] ) - ( v3[2] - v1[2] )*( v2[0] - v1[0] );
  n[2] = ( v2[0] - v1[0] )*( v3[1] - v1[1] ) - ( v3[0] - v1[0] )*( v2[1] - v1[1] );
  a = sqrt(n[0]*n[0]+n[1]*n[1]+n[2]*n[2])*0.5;
  const double invlen = 0.5/a;
  n[0]*=invlen;	n[1]*=invlen;	n[2]*=invlen;
}


static void MatVec3(const double m[9], const double x[3], double y[3]){
  y[0] = m[0]*x[0] + m[1]*x[1] + m[2]*x[2];
  y[1] = m[3]*x[0] + m[4]*x[1] + m[5]*x[2];
  y[2] = m[6]*x[0] + m[7]*x[1] + m[8]*x[2];
}

/*
static void VecMat3(const double x[3], const double m[9],  double y[3]){
  y[0] = m[0]*x[0] + m[3]*x[1] + m[6]*x[2];
  y[1] = m[1]*x[0] + m[4]*x[1] + m[7]*x[2];
  y[2] = m[2]*x[0] + m[5]*x[1] + m[8]*x[2];
}
 */
/*
static inline double Distance3D(const double p0[3], const double p1[3]){
  return sqrt( (p1[0]-p0[0])*(p1[0]-p0[0]) + (p1[1]-p0[1])*(p1[1]-p0[1]) + (p1[2]-p0[2])*(p1[2]-p0[2]) );
}
 */

static inline double Distance2D(const double p0[3], const double p1[3]){
  return sqrt( (p1[0]-p0[0])*(p1[0]-p0[0]) + (p1[1]-p0[1])*(p1[1]-p0[1]) );
}

static double Dot(const double p0[3], const double p1[3]){
  return p0[0]*p1[0] + p0[1]*p1[1] + p0[2]*p1[2];
}

static inline double largest(double x0, double x1, double x2) {
  double wmax = x0;
  wmax = (x1 > wmax) ? x1 : wmax;
  wmax = (x2 > wmax) ? x2 : wmax;
  return wmax;
}

void updateMinMaxXYZ
(double& x_min, double& x_max,
 double& y_min, double& y_max,
 double& z_min, double& z_max,
 /////
 double x, double y, double z)
{
  x_min = (x_min < x) ? x_min : x;
  x_max = (x_max > x) ? x_max : x;
  y_min = (y_min < y) ? y_min : y;
  y_max = (y_max > y) ? y_max : y;
  z_min = (z_min < z) ? z_min : z;
  z_max = (z_max > z) ? z_max : z;
}

void GetCenterWidth_MinMaxXYZ
(double& cx, double& cy, double& cz,
 double& wx, double& wy, double& wz,
 ////
 double x_min, double x_max,
 double y_min, double y_max,
 double z_min, double z_max)
{
  cx = (x_min+x_max)*0.5;
  cy = (y_min+y_max)*0.5;
  cz = (z_min+z_max)*0.5;
  wx = x_max-x_min;
  wy = y_max-y_min;
  wz = z_max-z_min;
}

////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

void RemoveUnreferencedPoints_MeshElem
(std::vector<double>& aXYZ1,
 std::vector<unsigned int>& aElem1,
 std::vector<int>& aMap01,
 int ndim,
 const std::vector<double>& aXYZ0,
 const std::vector<unsigned int>& aElem0)
{
  int np0 = aXYZ0.size()/ndim;
  aMap01.assign(np0,-2);
  for(unsigned int it=0;it<aElem0.size();++it){
    int ip = aElem0[it];
    aMap01[ip] = -1;
  }
  int npj = 0;
  for(int ip=0;ip<np0;++ip){
    if( aMap01[ip] == -2 ) continue;
    aMap01[ip] = npj;
    npj++;
  }
  aXYZ1.resize(npj*ndim);
  for(int ip=0;ip<np0;++ip){
    if( aMap01[ip] == -2 ) continue;
    int jp = aMap01[ip];
    for(int idim=0;idim<ndim;++idim){
      aXYZ1[jp*ndim+idim] = aXYZ0[ip*ndim+idim];
    }
  }
  aElem1.resize(aElem0.size());
  for(unsigned int it=0;it<aElem0.size();++it){
    int ip = aElem0[it];
    int jp = aMap01[ip];
    aElem1[it] = jp;
  }
}

void Normal_MeshTri3D
(double* aNorm,
 const double* aXYZ, int nXYZ,
 const unsigned int* aTri, int nTri)
{
  for(int i=0;i<nXYZ*3;i++){ aNorm[i] = 0; }
  for(int itri=0;itri<nTri;itri++){
    const int i0 = aTri[itri*3+0];
    const int i1 = aTri[itri*3+1];
    const int i2 = aTri[itri*3+2];
    if( i0 == -1 ){
      assert( i1 == -1 );
      assert( i2 == -1 );
      continue;
    }
    double p0[3] = {aXYZ[i0*3+0], aXYZ[i0*3+1], aXYZ[i0*3+2]};
    double p1[3] = {aXYZ[i1*3+0], aXYZ[i1*3+1], aXYZ[i1*3+2]};
    double p2[3] = {aXYZ[i2*3+0], aXYZ[i2*3+1], aXYZ[i2*3+2]};
    double un[3], area;    
    UnitNormalAreaTri3D(un,area, p0,p1,p2); 
    aNorm[i0*3+0] += un[0];  aNorm[i0*3+1] += un[1];  aNorm[i0*3+2] += un[2];
    aNorm[i1*3+0] += un[0];  aNorm[i1*3+1] += un[1];  aNorm[i1*3+2] += un[2];    
    aNorm[i2*3+0] += un[0];  aNorm[i2*3+1] += un[1];  aNorm[i2*3+2] += un[2];    
  }
  for(int ino=0;ino<nXYZ;ino++){
    const double n[3] = {aNorm[ino*3+0],aNorm[ino*3+1],aNorm[ino*3+2]};
    const double invlen = 1.0/Length3D(n);
    aNorm[ino*3+0] *= invlen;
    aNorm[ino*3+1] *= invlen;
    aNorm[ino*3+2] *= invlen;    
  }  
}



void Quality_MeshTri2D
(double& max_aspect, double& min_area,
 const double* aXY,
 const unsigned int* aTri, int nTri)
{
  max_aspect = 0;
  min_area = 0;
  for(int itri=0;itri<nTri;itri++){
    const int i0 = aTri[itri*3+0];
    const int i1 = aTri[itri*3+1];
    const int i2 = aTri[itri*3+2];
    const double* p0 = aXY+i0*2;
    const double* p1 = aXY+i1*2;
    const double* p2 = aXY+i2*2;
    const double area = TriArea2D(p0,p1,p2);
    const double len01 = Distance2D(p0,p1);
    const double len12 = Distance2D(p1,p2);
    const double len20 = Distance2D(p2,p0);
    const double len_ave = (len01+len12+len20)/3.0;
    const double aspect = len_ave * len_ave / area;
    if( itri == 0 ){
      max_aspect = aspect;
      min_area = area;
    }
    else{
      if( aspect > max_aspect ){ max_aspect = aspect; }
      if( area < min_area ){ min_area = area; }
    }
  }
}

////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////

void GetCenterWidth
(double& cx, double& cy, double& cz,
 double& wx, double& wy, double& wz,
 const int nXYZ, const double* paXYZ)
{
  if( paXYZ == 0 ){ cx=cy=cz=0; wx=wy=wz=1; return; }
  double x_min=paXYZ[0], x_max=paXYZ[0];
  double y_min=paXYZ[1], y_max=paXYZ[1];
  double z_min=paXYZ[2], z_max=paXYZ[2];
  for(int ino=0;ino<nXYZ;ino++){
    updateMinMaxXYZ(x_min,x_max, y_min,y_max, z_min,z_max,
                    paXYZ[ino*3+0], paXYZ[ino*3+1], paXYZ[ino*3+2]);
  }
  GetCenterWidth_MinMaxXYZ(cx,cy,cz, wx,wy,wz,
                           x_min,x_max, y_min,y_max, z_min,z_max);
}

void GetCenterWidth
(double& cx, double& cy, double& cz,
double& wx, double& wy, double& wz,
const std::vector<double>& aXYZ)
{
  const int np = (int)aXYZ.size()/3;
  if(np==0){ cx=cy=cz=0; wx=wy=wz=1; return; }
  double x_min=aXYZ[0], x_max=aXYZ[0];
  double y_min=aXYZ[1], y_max=aXYZ[1];
  double z_min=aXYZ[2], z_max=aXYZ[2];
  for (int ip=0; ip<np; ++ip){
    updateMinMaxXYZ(x_min,x_max, y_min,y_max, z_min,z_max,
                    aXYZ[ip*3+0], aXYZ[ip*3+1], aXYZ[ip*3+2]);
  }
  GetCenterWidth_MinMaxXYZ(cx,cy,cz, wx,wy,wz,
                           x_min,x_max, y_min,y_max, z_min,z_max);
}

void GetCenterWidth
(double cw[6],
 const std::vector<double>& aXYZ)
{
  GetCenterWidth(cw[0],cw[1],cw[2],cw[3],cw[4],cw[5],
                 aXYZ);
}

void MinMaxXYZ(double mm[6],
               const std::vector<double>& aXYZ)
{
  mm[0] = +1;
  mm[1] = -1;
  for(unsigned int ixyz=0;ixyz<aXYZ.size()/3;++ixyz){
    updateMinMaxXYZ(mm[0], mm[1], mm[2], mm[3], mm[4], mm[5],
                    aXYZ[ixyz*3+0], aXYZ[ixyz*3+1], aXYZ[ixyz*3+2]);
  }
}

void GetCenterWidthGroup
(double& cx, double& cy, double& cz,
 double& wx, double& wy, double& wz,
 ////
 const std::vector<double>& aXYZ,
 const std::vector<int>& aElem,
 const int nnoel,
 int igroup,
 const std::vector<int>& aIndGroup)
{
  const unsigned int nelem = aElem.size()/nnoel;
  assert( aElem.size() == nelem*nnoel );
  assert( aIndGroup.size() == nelem );
  bool is_ini = true;
  double x_min=0, x_max=0, y_min=0, y_max=0, z_min=0, z_max=0;
  for(unsigned int ielem=0;ielem<nelem;ielem++){
    if( aIndGroup[ielem] != igroup ){ continue; }
    for(int inotri=0;inotri<nnoel;inotri++){
      const int ip = aElem[ielem*3+inotri];
      if( is_ini ){
        x_min = x_max = aXYZ[ip*3+0];
        y_min = y_max = aXYZ[ip*3+1];
        z_min = z_max = aXYZ[ip*3+2];
        is_ini = false;
        continue;
      }
      updateMinMaxXYZ(x_min,x_max, y_min,y_max, z_min,z_max,
                      aXYZ[ip*3+0], aXYZ[ip*3+1], aXYZ[ip*3+2]);
    }
  }
  if (is_ini){ cx=cy=cz=0; wx=wy=wz=1; return; }
  GetCenterWidth_MinMaxXYZ(cx,cy,cz, wx,wy,wz,
                           x_min,x_max, y_min,y_max, z_min,z_max);
}

void GetCenterWidthGroup
(double& cx, double& cy, double& cz,
 double& wx, double& wy, double& wz,
 ////
 const std::vector<double>& aXYZ,
 const std::vector<int>& aElemInd,
 const std::vector<int>& aElem,
 int igroup,
 const std::vector<int>& aIndGroup)
{
  assert(!aElemInd.empty());
  const unsigned int nelem = aElemInd.size()-1;
  assert( aIndGroup.size() == nelem );
  bool is_ini = true;
  double x_min=0, x_max=0, y_min=0, y_max=0, z_min=0, z_max=0;
  for(unsigned int ielem=0;ielem<nelem;ielem++){
    if( aIndGroup[ielem] != igroup ){ continue; }
    for(int iip=aElemInd[ielem];iip<aElemInd[ielem+1];iip++){
      const int ip = aElem[iip];
      if( is_ini ){
        x_min = x_max = aXYZ[ip*3+0];
        y_min = y_max = aXYZ[ip*3+1];
        z_min = z_max = aXYZ[ip*3+2];
        is_ini = false;
        continue;
      }
      updateMinMaxXYZ(x_min,x_max, y_min,y_max, z_min,z_max,
                      aXYZ[ip*3+0], aXYZ[ip*3+1], aXYZ[ip*3+2]);
    }
  }
  if (is_ini){ cx=cy=cz=0; wx=wy=wz=1; return; }
  GetCenterWidth_MinMaxXYZ(cx,cy,cz, wx,wy,wz,
                           x_min,x_max, y_min,y_max, z_min,z_max);
}

void GetCenterWidth3DGroup
(double cw[6],
 ////
 const std::vector<double>& aXYZ,
 const std::vector<int>& aElemInd,
 const std::vector<int>& aElem,
 int igroup,
 const std::vector<int>& aIndGroup)
{
  GetCenterWidthGroup(cw[0],cw[1],cw[2], cw[3],cw[4],cw[5],
                      aXYZ,aElemInd,aElem, igroup, aIndGroup);
}


void GetCenterWidthLocal
(double& lcx, double& lcy, double& lcz,
 double& lwx, double& lwy, double& lwz,
 const std::vector<double>& aXYZ,
 const double lex[3],
 const double ley[3],
 const double lez[3])
{
  const int nno = (int)aXYZ.size()/3;
  if (nno==0){ lcx=lcy=lcz=0; lwx=lwy=lwz=1; return; }
  const double p0[3] = {aXYZ[0],aXYZ[1],aXYZ[2]};
  double x_min = Dot(p0,lex); double x_max = x_min;
  double y_min = Dot(p0,ley); double y_max = y_min;
  double z_min = Dot(p0,lez); double z_max = z_min;
  for (int ino = 0; ino<nno; ++ino){
    const double pi[3] = {aXYZ[ino*3+0],aXYZ[ino*3+1],aXYZ[ino*3+2]};
    updateMinMaxXYZ(x_min,x_max, y_min,y_max, z_min,z_max,
                    Dot(pi,lex), Dot(pi,ley), Dot(pi,lez));
  }
  GetCenterWidth_MinMaxXYZ(lcx,lcy,lcz, lwx,lwy,lwz,
                           x_min,x_max, y_min,y_max, z_min,z_max);
}

//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////

void Translate
(double tx, double ty, double tz,
std::vector<double>& aXYZ)
{
  const int nno = (int)aXYZ.size()/3;
  for (int ino = 0; ino<nno; ++ino){
    aXYZ[ino*3+0] += tx;
    aXYZ[ino*3+1] += ty;
    aXYZ[ino*3+2] += tz;
  }
}

void Scale
(double s,
std::vector<double>& aXYZ)
{
  const int nno = (int)aXYZ.size()/3;
  for (int ino = 0; ino<nno; ++ino){
    aXYZ[ino*3+0] *= s;
    aXYZ[ino*3+1] *= s;
    aXYZ[ino*3+2] *= s;
  }
}


void Translate
(std::vector<double>& aXYZ,
double tx, double ty, double tz)
{
  const int nXYZ = (int)aXYZ.size()/3;
  for (int ixyz = 0; ixyz<nXYZ; ixyz++){
    aXYZ[ixyz*3+0] += tx;
    aXYZ[ixyz*3+1] += ty;
    aXYZ[ixyz*3+2] += tz;
  }
}

void Rotate
(std::vector<double>& aXYZ,
double radx, double rady, double radz)
{
  const double phi = radx;
  const double theta = rady;
  const double psi = radz;
  ////
  const double mat[9] = {
    cos(psi)*cos(theta), cos(psi)*sin(theta)*sin(phi)-sin(psi)*cos(phi), cos(psi)*sin(theta)*cos(phi)+sin(psi)*sin(phi),
    sin(psi)*cos(theta), sin(psi)*sin(theta)*sin(phi)+cos(psi)*cos(phi), sin(psi)*sin(theta)*cos(phi)-cos(psi)*sin(phi),
    -sin(theta),         cos(theta)*sin(phi),                            cos(theta)*cos(phi) };
  ////  
  const int nXYZ = (int)aXYZ.size()/3;
  for (int ixyz = 0; ixyz<nXYZ; ++ixyz){
    double p[3] = { aXYZ[ixyz*3+0], aXYZ[ixyz*3+1], aXYZ[ixyz*3+2] };
    double res[3];  MatVec3(mat, p, res);
    aXYZ[ixyz*3+0] = res[0];
    aXYZ[ixyz*3+1] = res[1];
    aXYZ[ixyz*3+2] = res[2];
  }
}



double Size(const std::vector<double>& aXYZ){
  double cx, cy, cz, wx, wy, wz;
  GetCenterWidth(cx, cy, cz, wx, wy, wz, aXYZ);
  double wmax = largest(wx, wy, wz);
  return wmax;
}

void Normalize
(std::vector<double>& aXYZ,
 double s)
{
  double cx, cy, cz, wx, wy, wz;
  GetCenterWidth(cx, cy, cz, wx, wy, wz, aXYZ);
  Translate(-cx, -cy, -cz, aXYZ);
  double wmax = largest(wx, wy, wz);
  Scale(s/wmax, aXYZ);
}

void CenterOfGravity
(double& cgx, double& cgy, double& cgz,
const std::vector<double>& aXYZ)
{
  cgx = 0;
  cgy = 0;
  cgz = 0;
  int nXYZ = (int)aXYZ.size()/3;
  for (int ixyz = 0; ixyz<nXYZ; ixyz++){
    cgx += aXYZ[ixyz*3+0];
    cgy += aXYZ[ixyz*3+1];
    cgz += aXYZ[ixyz*3+2];
  }
  cgx /= nXYZ;
  cgy /= nXYZ;
  cgz /= nXYZ;
}

static double TetVolume3D(const double v1[3],
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


void CenterOfGravity_Solid
(double& cgx, double& cgy, double& cgz,
const std::vector<double>& aXYZ,
const std::vector<int>& aTri)
{ // center of gravity
  cgx = 0.0;
  cgy = 0.0;
  cgz = 0.0;
  double tw = 0;
  for (unsigned int itri = 0; itri<aTri.size()/3; itri++){
    unsigned int i1 = aTri[itri*3+0];
    unsigned int i2 = aTri[itri*3+1];
    unsigned int i3 = aTri[itri*3+2];
    const double q0[3] = { 0, 0, 0 };
    const double q1[3] = { aXYZ[i1*3+0], aXYZ[i1*3+1], aXYZ[i1*3+2] };
    const double q2[3] = { aXYZ[i2*3+0], aXYZ[i2*3+1], aXYZ[i2*3+2] };
    const double q3[3] = { aXYZ[i3*3+0], aXYZ[i3*3+1], aXYZ[i3*3+2] };
    double v = ::TetVolume3D(q0, q1, q2, q3);
    tw += v;
    cgx += (q0[0]+q1[0]+q2[0]+q3[0])*0.25*v;
    cgy += (q0[1]+q1[1]+q2[1]+q3[1])*0.25*v;
    cgz += (q0[2]+q1[2]+q2[2]+q3[2])*0.25*v;
  }
  cgx /= tw;
  cgy /= tw;
  cgz /= tw;
}

void CenterOfGravity_Shell
(double& cgx, double& cgy, double& cgz,
 const std::vector<double>& aXYZ,
 const std::vector<int>& aTri)
{ // center of gravity
  cgx = 0.0;
  cgy = 0.0;
  cgz = 0.0;
  double tw = 0;
  for (unsigned int itri = 0; itri<aTri.size()/3; itri++){
    unsigned int i1 = aTri[itri*3+0];
    unsigned int i2 = aTri[itri*3+1];
    unsigned int i3 = aTri[itri*3+2];
    const double q1[3] = { aXYZ[i1*3+0], aXYZ[i1*3+1], aXYZ[i1*3+2] };
    const double q2[3] = { aXYZ[i2*3+0], aXYZ[i2*3+1], aXYZ[i2*3+2] };
    const double q3[3] = { aXYZ[i3*3+0], aXYZ[i3*3+1], aXYZ[i3*3+2] };
    double a = TriArea3D(q1, q2, q3);
    tw += a;
    cgx += (q1[0]+q2[0]+q3[0])*0.333333*a;
    cgy += (q1[1]+q2[1]+q3[1])*0.333333*a;
    cgz += (q1[2]+q2[2]+q3[2])*0.333333*a;
  }
  cgx /= tw;
  cgy /= tw;
  cgz /= tw;
}


double CenterOfGravity_TriMsh3DFlg_Shell
(double& cgx, double& cgy, double& cgz,
 const std::vector<double>& aXYZ,
 const std::vector<int>& aTri,
 int iflg,
 const std::vector<int>& aFlg)
{
  cgx = 0.0;
  cgy = 0.0;
  cgz = 0.0;
  double tw = 0;
  for (unsigned int itri = 0; itri<aTri.size()/3; itri++){
    if( aFlg[itri] != iflg ) continue;
    const unsigned int i1 = aTri[itri*3+0];
    const unsigned int i2 = aTri[itri*3+1];
    const unsigned int i3 = aTri[itri*3+2];
    const double q1[3] = { aXYZ[i1*3+0], aXYZ[i1*3+1], aXYZ[i1*3+2] };
    const double q2[3] = { aXYZ[i2*3+0], aXYZ[i2*3+1], aXYZ[i2*3+2] };
    const double q3[3] = { aXYZ[i3*3+0], aXYZ[i3*3+1], aXYZ[i3*3+2] };
    double a = TriArea3D(q1, q2, q3);
    tw += a;
    cgx += (q1[0]+q2[0]+q3[0])*0.333333*a;
    cgy += (q1[1]+q2[1]+q3[1])*0.333333*a;
    cgz += (q1[2]+q2[2]+q3[2])*0.333333*a;
  }
  cgx /= tw;
  cgy /= tw;
  cgz /= tw;
  return tw;
}

void CenterOfGravity_Tri
(double& cgx, double& cgy, double& cgz,
 int itri,
 const std::vector<double>& aXYZ,
 const std::vector<int>& aTri)
{ // center of gravity
  assert( itri >= 0 && itri < (int)aTri.size()/3 );
  const int i1 = aTri[itri*3+0];
  const int i2 = aTri[itri*3+1];
  const int i3 = aTri[itri*3+2];
  const double q1[3] = { aXYZ[i1*3+0], aXYZ[i1*3+1], aXYZ[i1*3+2] };
  const double q2[3] = { aXYZ[i2*3+0], aXYZ[i2*3+1], aXYZ[i2*3+2] };
  const double q3[3] = { aXYZ[i3*3+0], aXYZ[i3*3+1], aXYZ[i3*3+2] };
//  double a = TriArea3D(q1, q2, q3);
  cgx = (q1[0]+q2[0]+q3[0])*0.333333;
  cgy = (q1[1]+q2[1]+q3[1])*0.333333;
  cgz = (q1[2]+q2[2]+q3[2])*0.333333;
}


///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////



void Translate
(double tx, double ty, double tz,
 const unsigned int nnode_, double* pXYZs_)
{
  for(unsigned int ino=0;ino<nnode_;ino++){
    pXYZs_[ino*3+0] += tx;
    pXYZs_[ino*3+1] += ty;
    pXYZs_[ino*3+2] += tz;
  }
}

void Scale
(double s,
 const unsigned int nnode_, double* pXYZs_)
{
  for(unsigned int ino=0;ino<nnode_;ino++){
    pXYZs_[ino*3+0] *= s;
    pXYZs_[ino*3+1] *= s;
    pXYZs_[ino*3+2] *= s;
  }  
}


/////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////

void MeshQuad2D_Grid
(std::vector<double>& aXYZ,
 std::vector<unsigned int>& aQuad,
 int nx, int ny)
{
  int np = (nx+1)*(ny+1);
  aXYZ.resize(np*2);
  for(int iy=0;iy<ny+1;++iy){
    for(int ix=0;ix<nx+1;++ix){
      int ip = iy*(nx+1)+ix;
      aXYZ[ip*2+0] = ix;
      aXYZ[ip*2+1] = iy;
    }
  }
  aQuad.resize(nx*ny*4);
  for(int iy=0;iy<ny;++iy){
    for(int ix=0;ix<nx;++ix){
      int iq = iy*nx+ix;
      aQuad[iq*4+0] = (iy+0)*(nx+1)+(ix+0);
      aQuad[iq*4+1] = (iy+0)*(nx+1)+(ix+1);
      aQuad[iq*4+2] = (iy+1)*(nx+1)+(ix+1);
      aQuad[iq*4+3] = (iy+1)*(nx+1)+(ix+0);
    }
  }
}

void MeshTri3D_Disk
(std::vector<double>& aXYZ, 
 std::vector<int>& aTri,
 double r, int nr, int nth)
{
  aXYZ.clear();
  aTri.clear();
  const double pi = 3.1415926535;
  { // make coordinates
    const int npo = 1+nr*nth;
    double dr = r/nr;
    double dth = 2.0*pi/nth;
    aXYZ.reserve(npo*3);
    aXYZ.push_back(0.0);
    aXYZ.push_back(0.0);
    aXYZ.push_back(0.0);
    for(int ir=1;ir<=nr;ir++){
      double ri = dr*ir;
      for(int ith=0;ith<nth;ith++){
        aXYZ.push_back(ri*cos(ith*dth));
        aXYZ.push_back(0);
        aXYZ.push_back(ri*sin(ith*dth));
      }
    }
  }
  int ntri = nth*(nr-1)*2+nth;
  aTri.reserve(ntri*3);
  for (int ith = 0; ith<nth; ith++){
    aTri.push_back(0);
    aTri.push_back((ith+1)%nth+1);
    aTri.push_back((ith+0)%nth+1);
  }
  for (int ir = 0; ir<nr-1; ir++){
    for (int ith = 0; ith<nth; ith++){
      int i1 = (ir+0)*nth+1+(ith+0)%nth;
      int i2 = (ir+0)*nth+1+(ith+1)%nth;
      int i3 = (ir+1)*nth+1+(ith+1)%nth;
      int i4 = (ir+1)*nth+1+(ith+0)%nth;
      aTri.push_back(i3);
      aTri.push_back(i1);
      aTri.push_back(i2);
      aTri.push_back(i4);
      aTri.push_back(i1);
      aTri.push_back(i3);
    }
  }
}



void MeshTri3D_OpenCylinder
(std::vector<double>& aXYZ, 
std::vector<int>& aTri,
double r, double l,
int nr, int nl)
{
  aXYZ.clear();
  aTri.clear();
  const double pi = 3.1415926535;
  double dl = l/nl;
  double dr = 2.0*pi/nr;
  const int npo = (nl+1)*nr;
  aXYZ.reserve(npo*3);
  for (int il = 0; il<nl+1; il++){
    double y0 = -0.5*l+il*dl;
    for (int ir = 0; ir<nr; ir++){
      double x0 = r*cos(dr*ir);
      double z0 = r*sin(dr*ir);
      aXYZ.push_back(x0);
      aXYZ.push_back(y0);
      aXYZ.push_back(z0);
    }
  }
  /////
  const int ntri = nl*nr*2;
  aTri.reserve(ntri*3);
  for (int il = 0; il<nl; il++){
    for (int ir = 0; ir<nr; ir++){
      const int i1 = (il+0)*nr+(ir+0)%nr;
      const int i2 = (il+0)*nr+(ir+1)%nr;
      const int i3 = (il+1)*nr+(ir+1)%nr;
      const int i4 = (il+1)*nr+(ir+0)%nr;
//      std::cout<<i1<<" "<<i2<<" "<<i3<<" "<<i4<<" "<<npo<<std::endl;
      aTri.push_back(i3);
      aTri.push_back(i2);
      aTri.push_back(i1);
      aTri.push_back(i4);
      aTri.push_back(i3);
      aTri.push_back(i1);
    }
  }
}

void MeshTri3D_ClosedCylinder
(std::vector<double>& aXYZ, 
std::vector<unsigned int>& aTri,
double r, double l,
int nlo, int nl)
{
  int nla = nl+2;
  aXYZ.clear();
  aTri.clear();
  if (nla<=1||nlo<=2){ return; }
  const double pi = 3.1415926535;
  double dl = l/nl;
  double dr = 2.0*pi/nlo;
  aXYZ.reserve((nlo*(nla-1)+2)*3);
  for (int ila = 0; ila<nla+1; ila++){
    double y0 = -0.5*l+dl*(ila-1);
    if (ila==0  ){ y0 = -0.5*l; }
    if (ila==nla){ y0 = +0.5*l; }
    for (int ilo = 0; ilo<nlo; ilo++){
      double x0 = r*cos(dr*ilo);
      double z0 = r*sin(dr*ilo);
      if (ila==0){
        aXYZ.push_back(0);
        aXYZ.push_back(y0);
        aXYZ.push_back(0);
        break;
      }
      else if(ila==nla){ 
        aXYZ.push_back(0);
        aXYZ.push_back(y0);
        aXYZ.push_back(0);
        break; 
      }
      else{
        aXYZ.push_back(x0);
        aXYZ.push_back(y0);
        aXYZ.push_back(z0);
      }
    }
  }
  /////
  int ntri = nlo*(nla-1)*2+nlo*2;
  aTri.reserve(ntri*3);
  for (int ilo = 0; ilo<nlo; ilo++){
    aTri.push_back(0);
    aTri.push_back((ilo+0)%nlo+1);
    aTri.push_back((ilo+1)%nlo+1);
  }
  for (int ila = 0; ila<nla-2; ila++){
    for (int ilo = 0; ilo<nlo; ilo++){
      int i1 = (ila+0)*nlo+1+(ilo+0)%nlo;
      int i2 = (ila+0)*nlo+1+(ilo+1)%nlo;
      int i3 = (ila+1)*nlo+1+(ilo+1)%nlo;
      int i4 = (ila+1)*nlo+1+(ilo+0)%nlo;
      aTri.push_back(i3);
      aTri.push_back(i2);
      aTri.push_back(i1);
      aTri.push_back(i4);
      aTri.push_back(i3);
      aTri.push_back(i1);
    }
  }
  for (int ilo = 0; ilo<nlo; ilo++){
    aTri.push_back(nlo*(nla-1)+1);
    aTri.push_back((nla-2)*nlo+1+(ilo+1)%nlo);
    aTri.push_back((nla-2)*nlo+1+(ilo+0)%nlo);
  }
/*
  for(int itri=0;itri<aTri.size()/3;itri++){
    for(int inotri=0;inotri<3;++inotri){
      const int i0 = aTri[itri*3+inotri];
      assert( i0 >=0 && i0 < aXYZ.size()/3 );
    }
  }
*/ 
}

void MeshTri3D_Sphere
(std::vector<double>& aXYZ,
 std::vector<unsigned int>& aTri,
 double r,
 int nla, int nlo)
{
  aXYZ.clear();
  aTri.clear();
  if( nla <= 1 || nlo <= 2 ){ return; }
  const double pi = 3.1415926535;
  double dl = pi/nla;
  double dr = 2.0*pi/nlo;
  aXYZ.reserve( (nlo*(nla-1)+2)*3 );
  for(int ila=0;ila<nla+1;ila++){
    double y0 = cos(dl*ila);
    double r0 = sin(dl*ila);
    for(int ilo=0;ilo<nlo;ilo++){
      double x0 = r0*sin(dr*ilo);
      double z0 = r0*cos(dr*ilo);
      aXYZ.push_back(r*x0);
      aXYZ.push_back(r*y0);
      aXYZ.push_back(r*z0);
      if( ila == 0 || ila == nla ){ break; }
    }
  }
  /////
  int ntri = nlo*(nla-1)*2+nlo*2;
  aTri.reserve(ntri*3);
  for(int ilo=0;ilo<nlo;ilo++){
    aTri.push_back(0);
    aTri.push_back((ilo+0)%nlo+1);
    aTri.push_back((ilo+1)%nlo+1);
  }
  for(int ila=0;ila<nla-2;ila++){
    for(int ilo=0;ilo<nlo;ilo++){
      int i1 = (ila+0)*nlo+1+(ilo+0)%nlo;
      int i2 = (ila+0)*nlo+1+(ilo+1)%nlo;
      int i3 = (ila+1)*nlo+1+(ilo+1)%nlo;
      int i4 = (ila+1)*nlo+1+(ilo+0)%nlo;
      aTri.push_back(i3);
      aTri.push_back(i2);
      aTri.push_back(i1);
      aTri.push_back(i4);
      aTri.push_back(i3);
      aTri.push_back(i1);
    }
  }
  for(int ilo=0;ilo<nlo;ilo++){
    aTri.push_back(nlo*(nla-1)+1);
    aTri.push_back((nla-2)*nlo+1+(ilo+1)%nlo);
    aTri.push_back((nla-2)*nlo+1+(ilo+0)%nlo);
  }
}

// p0: -x, -y, -z
// p1: +x, -y, -z
// p2: -x, +y, -z
// p3: +x, +y, -z
// p4: -x, -y, +z
// p5: +x, -y, +z
// p6: -x, +y, +z
// p7: +x, +y, +z
// f0: -x
// f1: +x
// f2: -y
// f3: +y
// f4: -z
// f5: +z
void SetTopoQuad_CubeVox(std::vector<unsigned int>& aQuad)
{
  aQuad.resize(6*4);
  aQuad[0*4+0] = 0;    aQuad[0*4+1] = 4;   aQuad[0*4+2] = 6;   aQuad[0*4+3] = 2;
  aQuad[1*4+0] = 1;    aQuad[1*4+1] = 3;   aQuad[1*4+2] = 7;   aQuad[1*4+3] = 5;
  aQuad[2*4+0] = 0;    aQuad[2*4+1] = 1;   aQuad[2*4+2] = 5;   aQuad[2*4+3] = 4;
  aQuad[3*4+0] = 2;    aQuad[3*4+1] = 6;   aQuad[3*4+2] = 7;   aQuad[3*4+3] = 3;
  aQuad[4*4+0] = 0;    aQuad[4*4+1] = 2;   aQuad[4*4+2] = 3;   aQuad[4*4+3] = 1;
  aQuad[5*4+0] = 4;    aQuad[5*4+1] = 5;   aQuad[5*4+2] = 7;   aQuad[5*4+3] = 6;
}

void MeshQuad3D_CubeVox
(std::vector<double>& aXYZ, std::vector<unsigned int>& aQuad,
 double x_min, double x_max,
 double y_min, double y_max,
 double z_min, double z_max)
{
  aXYZ.resize(0);
  aXYZ.reserve(8*3);
  aXYZ.push_back(x_min);    aXYZ.push_back(y_min);    aXYZ.push_back(z_min);
  aXYZ.push_back(x_max);    aXYZ.push_back(y_min);    aXYZ.push_back(z_min);
  aXYZ.push_back(x_min);    aXYZ.push_back(y_max);    aXYZ.push_back(z_min);
  aXYZ.push_back(x_max);    aXYZ.push_back(y_max);    aXYZ.push_back(z_min);
  aXYZ.push_back(x_min);    aXYZ.push_back(y_min);    aXYZ.push_back(z_max);
  aXYZ.push_back(x_max);    aXYZ.push_back(y_min);    aXYZ.push_back(z_max);
  aXYZ.push_back(x_min);    aXYZ.push_back(y_max);    aXYZ.push_back(z_max);
  aXYZ.push_back(x_max);    aXYZ.push_back(y_max);    aXYZ.push_back(z_max);
  SetTopoQuad_CubeVox(aQuad);
}

void MeshTri3D_Cube
(std::vector<double>& aXYZ,
 std::vector<unsigned int>& aTri,
 int n)
{
  aXYZ.clear();
  aTri.clear();
  if( n < 1 ){ return; }
  double r = 1.0/n;
  const int np = 4*n*(n+1)+(n-1)*(n-1)*2;
  aXYZ.reserve( np*3 );
  for(int iz=0;iz<n+1;++iz){ // height
    for(int ix=0;ix<n;++ix){
      aXYZ.push_back(-0.5+r*ix);
      aXYZ.push_back(-0.5);
      aXYZ.push_back(-0.5+r*iz);
    }
    for(int iy=0;iy<n;++iy){
      aXYZ.push_back(+0.5);
      aXYZ.push_back(-0.5+r*iy);
      aXYZ.push_back(-0.5+r*iz);
    }
    for(int ix=n;ix>0;--ix){
      aXYZ.push_back(-0.5+r*ix);
      aXYZ.push_back(+0.5);
      aXYZ.push_back(-0.5+r*iz);
    }
    for(int iy=n;iy>0;--iy){
      aXYZ.push_back(-0.5);
      aXYZ.push_back(-0.5+r*iy);
      aXYZ.push_back(-0.5+r*iz);
    }
  }
  for(int iy=1;iy<n;++iy){
    for(int ix=1;ix<n;++ix){
      aXYZ.push_back(-0.5+r*ix);
      aXYZ.push_back(-0.5+r*iy);
      aXYZ.push_back(-0.5);
    }
  }
  for(int iy=1;iy<n;++iy){
   for(int ix=1;ix<n;++ix){
      aXYZ.push_back(-0.5+r*ix);
      aXYZ.push_back(-0.5+r*iy);
      aXYZ.push_back(+0.5);
    }
  }
  /////
  int ntri = n*n*6*2;
  aTri.reserve(ntri*3);
  for(int iz=0;iz<n;++iz){
    for(int ixy=0;ixy<4*n;++ixy){
      int i0 = ixy          +4*n*iz;
      int i1 = (ixy+1)%(4*n)+4*n*iz;
      int i2 = (ixy+1)%(4*n)+4*n*(iz+1);
      int i3 = ixy          +4*n*(iz+1);
      aTri.push_back(i0);
      aTri.push_back(i1);
      aTri.push_back(i2);
      ///
      aTri.push_back(i2);
      aTri.push_back(i3);
      aTri.push_back(i0);
    }
  }
  // bottom
  for(int ix=0;ix<n;++ix){
    for(int iy=0;iy<n;++iy){
      int i0, i1, i2, i3;
      i0 = 4*n*(n+1) + (iy-1)*(n-1)+(ix-1);
      i1 = 4*n*(n+1) + (iy-1)*(n-1)+(ix+0);
      i2 = 4*n*(n+1) + (iy+0)*(n-1)+(ix+0);
      i3 = 4*n*(n+1) + (iy+0)*(n-1)+(ix-1);
      if( ix==0 ){
        i0 = (iy==0) ? 0 : 4*n-iy;
        i3 = 4*n-iy-1;
      }
      if( ix==n-1 ){
        i1 = n+iy;
        i2 = n+iy+1;
      }
      if( iy==0 ){
        i0 = ix;
        i1 = ix+1;
      }
      if( iy==n-1 ){
        i2 = 3*n-ix-1;
        i3 = 3*n-ix+0;
      }
      aTri.push_back(i1);
      aTri.push_back(i0);
      aTri.push_back(i2);
      ///
      aTri.push_back(i3);
      aTri.push_back(i2);
      aTri.push_back(i0);
    }
  }
  // top
  int nps  = 4*n*(n+1); // side vertex
  int nps0 = 4*n*n; // side vertex
  for(int ix=0;ix<n;++ix){
    for(int iy=0;iy<n;++iy){
      int i0, i1, i2, i3;
      i0 = nps + (n-1)*(n-1) + (iy-1)*(n-1)+(ix-1);
      i1 = nps + (n-1)*(n-1) + (iy-1)*(n-1)+(ix+0);
      i2 = nps + (n-1)*(n-1) + (iy+0)*(n-1)+(ix+0);
      i3 = nps + (n-1)*(n-1) + (iy+0)*(n-1)+(ix-1);
      if( ix==0 ){
        i0 = (iy==0) ? nps0 : nps0+4*n-iy;
        i3 = nps0+4*n-iy-1;
      }
      if( ix==n-1 ){
        i1 = nps0+n+iy;
        i2 = nps0+n+iy+1;
      }
      if( iy==0 ){
        i0 = nps0+ix;
        i1 = nps0+ix+1;
      }
      if( iy==n-1 ){
        i2 = nps0+3*n-ix-1;
        i3 = nps0+3*n-ix+0;
      }
      aTri.push_back(i0);
      aTri.push_back(i1);
      aTri.push_back(i2);
      ///
      aTri.push_back(i2);
      aTri.push_back(i3);
      aTri.push_back(i0);
    }
  }
}


void MeshTri3D_Icosahedron
(std::vector<double>& aXYZ,
 std::vector<unsigned int>& aTri)
{
  double p = (1+sqrt(5))*0.5;
  aXYZ.resize(12*3);
  aXYZ[ 0*3+0]= 0;    aXYZ[ 0*3+1]=-1;   aXYZ[ 0*3+2]=-p;
  aXYZ[ 1*3+0]= 0;    aXYZ[ 1*3+1]=-1;   aXYZ[ 1*3+2]=+p;
  aXYZ[ 2*3+0]= 0;    aXYZ[ 2*3+1]=+1;   aXYZ[ 2*3+2]=-p;
  aXYZ[ 3*3+0]= 0;    aXYZ[ 3*3+1]=+1;   aXYZ[ 3*3+2]=+p;
  aXYZ[ 4*3+0]=-p;    aXYZ[ 4*3+1]= 0;   aXYZ[ 4*3+2]=-1;
  aXYZ[ 5*3+0]=+p;    aXYZ[ 5*3+1]= 0;   aXYZ[ 5*3+2]=-1;
  aXYZ[ 6*3+0]=-p;    aXYZ[ 6*3+1]= 0;   aXYZ[ 6*3+2]=+1;
  aXYZ[ 7*3+0]=+p;    aXYZ[ 7*3+1]= 0;   aXYZ[ 7*3+2]=+1;
  aXYZ[ 8*3+0]=-1;    aXYZ[ 8*3+1]=-p;   aXYZ[ 8*3+2]= 0;
  aXYZ[ 9*3+0]=-1;    aXYZ[ 9*3+1]=+p;   aXYZ[ 9*3+2]= 0;
  aXYZ[10*3+0]=+1;    aXYZ[10*3+1]=-p;   aXYZ[10*3+2]= 0;
  aXYZ[11*3+0]=+1;    aXYZ[11*3+1]=+p;   aXYZ[11*3+2]= 0;
  /////
  aTri.resize(20*3);
  aTri[ 0*3+0]= 7; aTri[ 0*3+1]=11; aTri[ 0*3+2]= 3;
  aTri[ 1*3+0]=11; aTri[ 1*3+1]= 9; aTri[ 1*3+2]= 3;
  aTri[ 2*3+0]= 9; aTri[ 2*3+1]= 6; aTri[ 2*3+2]= 3;
  aTri[ 3*3+0]= 6; aTri[ 3*3+1]= 1; aTri[ 3*3+2]= 3;
  aTri[ 4*3+0]= 1; aTri[ 4*3+1]= 7; aTri[ 4*3+2]= 3;
  /////
  aTri[ 5*3+0]= 2; aTri[ 5*3+1]= 5; aTri[ 5*3+2]= 0;
  aTri[ 6*3+0]= 4; aTri[ 6*3+1]= 2; aTri[ 6*3+2]= 0;
  aTri[ 7*3+0]= 8; aTri[ 7*3+1]= 4; aTri[ 7*3+2]= 0;
  aTri[ 8*3+0]=10; aTri[ 8*3+1]= 8; aTri[ 8*3+2]= 0;
  aTri[ 9*3+0]= 5; aTri[ 9*3+1]=10; aTri[ 9*3+2]= 0;
  /////
  aTri[10*3+0]=11; aTri[10*3+1]= 7; aTri[10*3+2]= 5;
  aTri[11*3+0]= 9; aTri[11*3+1]=11; aTri[11*3+2]= 2;
  aTri[12*3+0]= 6; aTri[12*3+1]= 9; aTri[12*3+2]= 4;
  aTri[13*3+0]= 1; aTri[13*3+1]= 6; aTri[13*3+2]= 8;
  aTri[14*3+0]= 7; aTri[14*3+1]= 1; aTri[14*3+2]=10;
  /////
  aTri[15*3+0]= 5; aTri[15*3+1]= 2; aTri[15*3+2]=11;
  aTri[16*3+0]= 2; aTri[16*3+1]= 4; aTri[16*3+2]= 9;
  aTri[17*3+0]= 4; aTri[17*3+1]= 8; aTri[17*3+2]= 6;
  aTri[18*3+0]= 8; aTri[18*3+1]=10; aTri[18*3+2]= 1;
  aTri[19*3+0]=10; aTri[19*3+1]= 5; aTri[19*3+2]= 7;
}

void SetTopology_ExtrudeTri2Tet
(unsigned int* aTet,
 int nXY,
 const unsigned int* aTri, int nTri,
 int nlayer)
{
  for(int il=0;il<nlayer;++il){
    for(int itri=0;itri<nTri;++itri){
      int ip0=-1, ip1=-1, ip2=-1;
      {
        const int i0 = aTri[itri*3+0];
        const int i1 = aTri[itri*3+1];
        const int i2 = aTri[itri*3+2];
        if( i0 > i1 && i0 > i2 ){ ip0=i0; ip1=i1; ip2=i2; }
        if( i1 > i0 && i1 > i2 ){ ip0=i1; ip1=i2; ip2=i0; }
        if( i2 > i0 && i2 > i1 ){ ip0=i2; ip1=i0; ip2=i1; }
        assert(ip0!=-1);
      }
      const int aIQ[6] = {
        (il+0)*nXY+ip0, (il+1)*nXY+ip0,
        (il+0)*nXY+ip1, (il+1)*nXY+ip1,
        (il+0)*nXY+ip2, (il+1)*nXY+ip2 };
      aTet[il*nTri*12+itri*12+4*0+0] = aIQ[0];
      aTet[il*nTri*12+itri*12+4*0+1] = aIQ[2];
      aTet[il*nTri*12+itri*12+4*0+2] = aIQ[4];
      aTet[il*nTri*12+itri*12+4*0+3] = aIQ[1];
      if( ip1 > ip2 ){
        aTet[il*nTri*12+itri*12+4*1+0] = aIQ[1];
        aTet[il*nTri*12+itri*12+4*1+1] = aIQ[3];
        aTet[il*nTri*12+itri*12+4*1+2] = aIQ[2];
        aTet[il*nTri*12+itri*12+4*1+3] = aIQ[4];
        aTet[il*nTri*12+itri*12+4*2+0] = aIQ[1];
        aTet[il*nTri*12+itri*12+4*2+1] = aIQ[3];
        aTet[il*nTri*12+itri*12+4*2+2] = aIQ[4];
        aTet[il*nTri*12+itri*12+4*2+3] = aIQ[5];
      }
      else{
        aTet[il*nTri*12+itri*12+4*1+0] = aIQ[1];
        aTet[il*nTri*12+itri*12+4*1+1] = aIQ[2];
        aTet[il*nTri*12+itri*12+4*1+2] = aIQ[5];
        aTet[il*nTri*12+itri*12+4*1+3] = aIQ[3];
        aTet[il*nTri*12+itri*12+4*2+0] = aIQ[1];
        aTet[il*nTri*12+itri*12+4*2+1] = aIQ[2];
        aTet[il*nTri*12+itri*12+4*2+2] = aIQ[4];
        aTet[il*nTri*12+itri*12+4*2+3] = aIQ[5];
      }
    }
  }
}

void ExtrudeTri2Tet
(int nlayer, double h,
 std::vector<double>& aXYZ,
 std::vector<unsigned int>& aTet,
 const std::vector<double>& aXY,
 const std::vector<unsigned int>& aTri)
{
  const int nXY = (int)aXY.size()/2;
  const int nTri = (int)aTri.size()/3;
  aXYZ.resize(nXY*(nlayer+1)*3);
  for(int il=0;il<nlayer+1;++il){
    for(int ixy=0;ixy<nXY;ixy++){
      aXYZ[il*nXY*3+ixy*3+0] = aXY[ixy*2+0];
      aXYZ[il*nXY*3+ixy*3+1] = aXY[ixy*2+1];
      aXYZ[il*nXY*3+ixy*3+2] = il*h;
    }
  }
  aTet.resize(nTri*nlayer*3*4);
  SetTopology_ExtrudeTri2Tet(aTet.data(),
                             nXY,aTri.data(),nTri,nlayer);
}

//////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////

void LaplacianSmoothing
(std::vector<double>& aXYZ,
 const std::vector<int>& aTri,
 const std::vector<int>& elsup_ind,
 const std::vector<int>& elsup)
{
  for(unsigned int ip=0;ip<aXYZ.size()/3;++ip){
    double sum_area = 0.0;
    double pcnt[3] = {0,0,0};
    for(int ielsup=elsup_ind[ip];ielsup<elsup_ind[ip+1];++ielsup){
      assert( ielsup < (int)elsup.size() );
      int iel = elsup[ielsup];
      assert( iel>=0 && iel<(int)aTri.size()/3 );
      int i0 = aTri[iel*3+0];
      int i1 = aTri[iel*3+1];
      int i2 = aTri[iel*3+2];
      double aP[3][3] = {
        { aXYZ[i0*3+0],aXYZ[i0*3+1],aXYZ[i0*3+2] },
        { aXYZ[i1*3+0],aXYZ[i1*3+1],aXYZ[i1*3+2] },
        { aXYZ[i2*3+0],aXYZ[i2*3+1],aXYZ[i2*3+2] } };
      double area = TriArea3D(aP[0],aP[1],aP[2]);
      sum_area += area;
      pcnt[0] += area*(aP[0][0]+aP[1][0]+aP[2][0])/3.0;
      pcnt[1] += area*(aP[0][1]+aP[1][1]+aP[2][1])/3.0;
      pcnt[2] += area*(aP[0][2]+aP[1][2]+aP[2][2])/3.0;
    }
    pcnt[0] /= sum_area;
    pcnt[1] /= sum_area;
    pcnt[2] /= sum_area;
    aXYZ[ip*3+0] = pcnt[0];
    aXYZ[ip*3+1] = pcnt[1];
    aXYZ[ip*3+2] = pcnt[2];
  }
}

/*
void LaplacianSmoothing_Cotan
(std::vector<double>& aXYZ,
 const std::vector<int>& aTri,
 const std::vector<int>& elsup_ind,
 const std::vector<int>& elsup)
{
  const int np = (int)aXYZ.size()/3;
  std::vector<double> aTmp(np,0.0);
  //  std::vector<double> aDist = aXYZ;
  for(int ip=0;ip<np;++ip){
    const CVector3 p0(aXYZ[ip*3+0],aXYZ[ip*3+1],aXYZ[ip*3+2]);
    double volo_area = 0;
    ////////////////////////////////////////////////////////
    for(int ielsup=elsup_ind[ip];ielsup<elsup_ind[ip+1];++ielsup){
      assert( ielsup < elsup.size() );
      int iel = elsup[ielsup];
      assert( iel>=0 && iel<aTri.size()/3 );
      const int i0 = aTri[iel*3+0];
      const int i1 = aTri[iel*3+1];
      const int i2 = aTri[iel*3+2];
      const int aIP[3] = {i0,i1,i2};
      const CVector3 aP[3] = {
        CVector3(aXYZ[i0*3+0],aXYZ[i0*3+1],aXYZ[i0*3+2]),
        CVector3(aXYZ[i1*3+0],aXYZ[i1*3+1],aXYZ[i1*3+2]),
        CVector3(aXYZ[i2*3+0],aXYZ[i2*3+1],aXYZ[i2*3+2]) };
      const double area = TriArea(aP[0],aP[1],aP[2]);
      for(int ie=0;ie<3;++ie){
        if( aIP[ie] == ip ) continue;
        const int ie1 = (ie+1)%3;
        const int ie2 = (ie+2)%3;
        assert( aIP[ie1] == ip || aIP[ie2] == ip );
        const CVector3 v01 = aP[ie1]-aP[ie];
        const CVector3 v02 = aP[ie2]-aP[ie];
        const CVector3 v12 = aP[ie2]-aP[ie1];
        const double cot102 = (v01*v02)/area*2;
        if( v01*v02>0 && v02*v12>0 && v01*v12<0 ){
          volo_area += v12.DLength()*cot102/0.125;
        }
        else if( v01*v02>0 ){
          volo_area += area*0.25;
        }
        else{
          volo_area += area*0.5;
        }
        aTmp[ aIP[ie1] ] += cot102;
        aTmp[ aIP[ie2] ] += cot102;
        //        std::cout << ip << " " << iel << " " << ie << " " << cot102 << std::endl;
      }
    }
    ////////////////////////////////////////////////////////
    CVector3 pcnt(0,0,0);
    double sum = 0.0;
    for(int ielsup=elsup_ind[ip];ielsup<elsup_ind[ip+1];++ielsup){
      assert( ielsup < elsup.size() );
      int iel = elsup[ielsup];
      assert( iel>=0 && iel<aTri.size()/3 );
      const int i0 = aTri[iel*3+0];
      const int i1 = aTri[iel*3+1];
      const int i2 = aTri[iel*3+2];
      const int aIP[3] = {i0,i1,i2};
      for(int ie=0;ie<3;++ie){
        if( aIP[ie] != ip ) continue;
        const int ie1 = (ie+1)%3;
        const int ie2 = (ie+2)%3;
        const int j1 = aIP[ie1];
        const int j2 = aIP[ie2];
        assert( j1 != ip );
        assert( j2 != ip );
        pcnt += aTmp[j1]*(CVector3(aXYZ[j1*3+0],aXYZ[j1*3+1],aXYZ[j1*3+2])-p0)*0.5;
        pcnt += aTmp[j2]*(CVector3(aXYZ[j2*3+0],aXYZ[j2*3+1],aXYZ[j2*3+2])-p0)*0.5;
      }
    }
    pcnt /= volo_area*2;
    double eps = 0.01;
    aXYZ[ip*3+0] += pcnt.x*eps;
    aXYZ[ip*3+1] += pcnt.y*eps;
    aXYZ[ip*3+2] += pcnt.z*eps;
    ////////////////////////////////////////////////////////
    for(int ielsup=elsup_ind[ip];ielsup<elsup_ind[ip+1];++ielsup){
      assert( ielsup < elsup.size() );
      const int iel = elsup[ielsup];
      assert( iel>=0 && iel<aTri.size()/3 );
      const int i0 = aTri[iel*3+0];
      const int i1 = aTri[iel*3+1];
      const int i2 = aTri[iel*3+2];
      aTmp[i0] = 0.0;
      aTmp[i1] = 0.0;
      aTmp[i2] = 0.0;
    }
  }
  //  aXYZ = aDist;
}
 */


double SolidAngleTri3D
(const double v1[3],
 const double v2[3],
 const double v3[3])
{
  double l1 = Length3D(v1);
  double l2 = Length3D(v2);
  double l3 = Length3D(v3);
  double crs_v1_v2[3]; Cross3D(crs_v1_v2,v1,v2);
  double den = Dot(crs_v1_v2,v3);
  double num = l1*l2*l3+(Dot(v1,v2))*l3+(Dot(v2,v3))*l1+(Dot(v3,v1))*l2;
  double tho = den/num;
  double v = atan(tho);
  if (v<0){ v += 2*M_PI; }
  v *= 2;
  return v;
}

void makeSolidAngle
(std::vector<double>& aSolidAngle,
 const std::vector<double>& aXYZ,
 const std::vector<int>& aTri,
 const std::vector<double>& aNorm, 
 std::vector<int>& elsup_ind,
 std::vector<int>& elsup)
{
  const int nXYZ = (int)aXYZ.size()/3;
  /*
  std::vector<double> aNorm;
  MakeNormal(aNorm, aXYZ, aTri);
  std::vector<int> elsup_ind,elsup;
  makeTriSurroundingPoint(elsup_ind,elsup,
                          aTri,(int)aXYZ.size()/3);
   */
  aSolidAngle.resize(nXYZ);
  for(int ip=0;ip<nXYZ;++ip){
    const double n0[3] = {aNorm[ip*3+0], aNorm[ip*3+1], aNorm[ip*3+2]};
    const double p0[3] = {aXYZ[ip*3+0],aXYZ[ip*3+1],aXYZ[ip*3+2]};
    double sa = 0;
    for(int ielsup=elsup_ind[ip];ielsup<elsup_ind[ip+1];++ielsup){
      int itri0 = elsup[ielsup];
      assert( itri0 >= 0 && itri0 < (int)aTri.size()/3 );
      int inotri0 = -1;
      for(int i=0;i<3;++i){
        if( aTri[itri0*3+i] != ip ){ continue; }
        inotri0 = i;
        break;
      }
      int inotri1 = (inotri0+1)%3;
      int inotri2 = (inotri0+2)%3;
      int ip1 = aTri[itri0*3+inotri1];
      int ip2 = aTri[itri0*3+inotri2];
      const double p1[3] = {aXYZ[ip1*3+0],aXYZ[ip1*3+1],aXYZ[ip1*3+2]};
      const double p2[3] = {aXYZ[ip2*3+0],aXYZ[ip2*3+1],aXYZ[ip2*3+2]};
      const double p10[3] = {p1[0]-p0[0],p1[1]-p0[1],p1[2]-p0[2]};
      const double p20[3] = {p2[0]-p0[0],p2[1]-p0[1],p2[2]-p0[2]};
      sa += SolidAngleTri3D(p10,p20,n0);
    }
    if( elsup_ind[ip+1]-elsup_ind[ip] == 0 ){
      sa = -1.0; // floting point: negative
    }
    aSolidAngle[ip] = sa;
  }
}


void MassLumped_Tet3D
(double* aMassMatrixLumped,
 double rho,
 const double* aXYZ, int nXYZ,
 const unsigned int* aTet, int nTet)
{
  for(int i=0;i<nXYZ;++i){ aMassMatrixLumped[i] = 0.0; }
  for(int it=0;it<nTet;++it){
    const int i0 = aTet[it*4+0]; assert(i0>=0&&i0<nXYZ);
    const int i1 = aTet[it*4+1]; assert(i1>=0&&i1<nXYZ);
    const int i2 = aTet[it*4+2]; assert(i2>=0&&i2<nXYZ);
    const int i3 = aTet[it*4+3]; assert(i3>=0&&i3<nXYZ);
    const double* p0 = aXYZ+i0*3;
    const double* p1 = aXYZ+i1*3;
    const double* p2 = aXYZ+i2*3;
    const double* p3 = aXYZ+i3*3;
    const double v0123 = TetVolume3D(p0, p1, p2, p3);
    aMassMatrixLumped[i0] += 0.25*rho*v0123;
    aMassMatrixLumped[i1] += 0.25*rho*v0123;
    aMassMatrixLumped[i2] += 0.25*rho*v0123;
    aMassMatrixLumped[i3] += 0.25*rho*v0123;
  }
}

// TODO: make this handle open surface (average face & edge independently)
void SubdivisionPoints_QuadCatmullClark
(std::vector<double>& aXYZ1,
 ///
 const std::vector<unsigned int>& aQuad1,
 const std::vector<int>& aEdgeFace0,
 const std::vector<int>& psupIndQuad0,
 const std::vector<int>& psupQuad0,
 const unsigned int* aQuad0, int nQuad0,
 const double* aXYZ0, int nXYZ0)
{
  /*
  std::vector<int> aEdgeFace0;
  std::vector<int> psupIndQuad0, psupQuad0;

  QuadSubdiv(aQuad1,
             psupIndQuad0,psupQuad0, aEdgeFace0,
             aQuad0, nv0);
   */
  const int nv0 = nXYZ0;
  const int ne0 = (int)psupQuad0.size();
  const int nq0 = nQuad0;
  assert( (int)aEdgeFace0.size() == ne0*4 );
  aXYZ1.resize((nv0+ne0+nq0)*3);
  std::vector<int> aW(nv0,0);
  for(int iv=0;iv<nv0;++iv){
    aXYZ1[iv*3+0] = 0;
    aXYZ1[iv*3+1] = 0;
    aXYZ1[iv*3+2] = 0;
  }
  for(int iq=0;iq<nq0;++iq){
    const int iv0 = aQuad0[iq*4+0];
    const int iv1 = aQuad0[iq*4+1];
    const int iv2 = aQuad0[iq*4+2];
    const int iv3 = aQuad0[iq*4+3];
    double p0x = (aXYZ0[iv0*3+0] + aXYZ0[iv1*3+0] + aXYZ0[iv2*3+0] + aXYZ0[iv3*3+0])*0.25;
    double p0y = (aXYZ0[iv0*3+1] + aXYZ0[iv1*3+1] + aXYZ0[iv2*3+1] + aXYZ0[iv3*3+1])*0.25;
    double p0z = (aXYZ0[iv0*3+2] + aXYZ0[iv1*3+2] + aXYZ0[iv2*3+2] + aXYZ0[iv3*3+2])*0.25;
    aXYZ1[(nv0+ne0+iq)*3+0] = p0x;
    aXYZ1[(nv0+ne0+iq)*3+1] = p0y;
    aXYZ1[(nv0+ne0+iq)*3+2] = p0z;
    const int aIV[4] = { iv0, iv1, iv2, iv3 };
    for(int iiv=0;iiv<4;++iiv){
      int jv0 = aIV[iiv];
      aXYZ1[jv0*3+0] += p0x;
      aXYZ1[jv0*3+1] += p0y;
      aXYZ1[jv0*3+2] += p0z;
      aW[jv0] += 1;
    }
  }
  for(int ie=0;ie<ne0;++ie){
    int iv0 = aEdgeFace0[ie*4+0];
    int iv1 = aEdgeFace0[ie*4+1];
    int iq0 = aEdgeFace0[ie*4+2];
    int iq1 = aEdgeFace0[ie*4+3];
    aXYZ1[(nv0+ie)*3+0] = (aXYZ0[iv0*3+0]+aXYZ0[iv1*3+0]+aXYZ1[(nv0+ne0+iq0)*3+0]+aXYZ1[(nv0+ne0+iq1)*3+0])*0.25;
    aXYZ1[(nv0+ie)*3+1] = (aXYZ0[iv0*3+1]+aXYZ0[iv1*3+1]+aXYZ1[(nv0+ne0+iq0)*3+1]+aXYZ1[(nv0+ne0+iq1)*3+1])*0.25;
    aXYZ1[(nv0+ie)*3+2] = (aXYZ0[iv0*3+2]+aXYZ0[iv1*3+2]+aXYZ1[(nv0+ne0+iq0)*3+2]+aXYZ1[(nv0+ne0+iq1)*3+2])*0.25;
    aXYZ1[iv0*3+0] += aXYZ0[iv0*3+0] + aXYZ0[iv1*3+0];
    aXYZ1[iv0*3+1] += aXYZ0[iv0*3+1] + aXYZ0[iv1*3+1];
    aXYZ1[iv0*3+2] += aXYZ0[iv0*3+2] + aXYZ0[iv1*3+2];
    aXYZ1[iv1*3+0] += aXYZ0[iv0*3+0] + aXYZ0[iv1*3+0];
    aXYZ1[iv1*3+1] += aXYZ0[iv0*3+1] + aXYZ0[iv1*3+1];
    aXYZ1[iv1*3+2] += aXYZ0[iv0*3+2] + aXYZ0[iv1*3+2];
  }
  for(int iv=0;iv<nv0;++iv){
    const int iw = aW[iv];
    if( iw == 0 ){ continue; }
    const double tmp0 = 1.0/(iw*iw);
    aXYZ1[iv*3+0] *= tmp0;
    aXYZ1[iv*3+1] *= tmp0;
    aXYZ1[iv*3+2] *= tmp0;
    const double tmp1 = (iw-3.0)/(iw);
    aXYZ1[iv*3+0] += tmp1*aXYZ0[iv*3+0];
    aXYZ1[iv*3+1] += tmp1*aXYZ0[iv*3+1];
    aXYZ1[iv*3+2] += tmp1*aXYZ0[iv*3+2];
  }
}

void SubdivisionPoints_Quad
(std::vector<double>& aXYZ1,
 ///
 const std::vector<int>& aQuad1,
 const std::vector<int>& aEdgeFace0,
 const std::vector<int>& psupIndQuad0,
 const std::vector<int>& psupQuad0,
 const std::vector<int>& aQuad0,
 const std::vector<double>& aXYZ0)
{
  /*
   std::vector<int> aEdgeFace0;
   std::vector<int> psupIndQuad0, psupQuad0;
   
   QuadSubdiv(aQuad1,
   psupIndQuad0,psupQuad0, aEdgeFace0,
   aQuad0, nv0);
   */
  const int nv0 = (int)aXYZ0.size()/3;
  const int ne0 = (int)psupQuad0.size();
  const int nq0 = (int)aQuad0.size()/4;
  assert( (int)aEdgeFace0.size() == ne0*4 );
  aXYZ1.resize((nv0+ne0+nq0)*3);
  for(int iv=0;iv<nv0;++iv){
    aXYZ1[iv*3+0] = aXYZ0[iv*3+0];
    aXYZ1[iv*3+1] = aXYZ0[iv*3+1];
    aXYZ1[iv*3+2] = aXYZ0[iv*3+2];
  }
  for(int ie=0;ie<ne0;++ie){
    const int iv0 = aEdgeFace0[ie*4+0];
    const int iv1 = aEdgeFace0[ie*4+1];
    aXYZ1[(nv0+ie)*3+0] = (aXYZ0[iv0*3+0] + aXYZ0[iv1*3+0])*0.5;
    aXYZ1[(nv0+ie)*3+1] = (aXYZ0[iv0*3+1] + aXYZ0[iv1*3+1])*0.5;
    aXYZ1[(nv0+ie)*3+2] = (aXYZ0[iv0*3+2] + aXYZ0[iv1*3+2])*0.5;
  }
  for(int iq=0;iq<nq0;++iq){
    const int iv0 = aQuad0[iq*4+0];
    const int iv1 = aQuad0[iq*4+1];
    const int iv2 = aQuad0[iq*4+2];
    const int iv3 = aQuad0[iq*4+3];
    aXYZ1[(nv0+ne0+iq)*3+0] = (aXYZ0[iv0*3+0] + aXYZ0[iv1*3+0] + aXYZ0[iv2*3+0] + aXYZ0[iv3*3+0])*0.25;
    aXYZ1[(nv0+ne0+iq)*3+1] = (aXYZ0[iv0*3+1] + aXYZ0[iv1*3+1] + aXYZ0[iv2*3+1] + aXYZ0[iv3*3+1])*0.25;
    aXYZ1[(nv0+ne0+iq)*3+2] = (aXYZ0[iv0*3+2] + aXYZ0[iv1*3+2] + aXYZ0[iv2*3+2] + aXYZ0[iv3*3+2])*0.25;
  }
}


void SubdivisionPoints_Hex
(std::vector<double>& aXYZ1,
 ///
 const std::vector<int>& psupIndHex0,
 const std::vector<int>& psupHex0,
 const std::vector<unsigned int>& aQuadHex0,
 const unsigned int* aHex0, int nHex0,
 const double* aXYZ0, int nXYZ0)
{
  const int nv0 = nXYZ0;
  const int ne0 = (int)psupHex0.size();
  const int nq0 = (int)aQuadHex0.size()/4;
  const int nh0 = nHex0;
  aXYZ1.resize((nv0+ne0+nq0+nh0)*3);
  for(int iv=0;iv<nv0;++iv){
    aXYZ1[iv*3+0] = aXYZ0[iv*3+0];
    aXYZ1[iv*3+1] = aXYZ0[iv*3+1];
    aXYZ1[iv*3+2] = aXYZ0[iv*3+2];
  }
  for(int iv=0;iv<nv0;++iv){
    for(int ipsup=psupIndHex0[iv];ipsup<psupIndHex0[iv+1];++ipsup){
      int jv = psupHex0[ipsup];
      aXYZ1[(nv0+ipsup)*3+0] = (aXYZ0[iv*3+0] + aXYZ0[jv*3+0])*0.5;
      aXYZ1[(nv0+ipsup)*3+1] = (aXYZ0[iv*3+1] + aXYZ0[jv*3+1])*0.5;
      aXYZ1[(nv0+ipsup)*3+2] = (aXYZ0[iv*3+2] + aXYZ0[jv*3+2])*0.5;
    }
  }
  for(int iq=0;iq<nq0;++iq){
    const int iv0 = aQuadHex0[iq*4+0];
    const int iv1 = aQuadHex0[iq*4+1];
    const int iv2 = aQuadHex0[iq*4+2];
    const int iv3 = aQuadHex0[iq*4+3];
    aXYZ1[(nv0+ne0+iq)*3+0] = (aXYZ0[iv0*3+0] + aXYZ0[iv1*3+0] + aXYZ0[iv2*3+0] + aXYZ0[iv3*3+0])*0.25;
    aXYZ1[(nv0+ne0+iq)*3+1] = (aXYZ0[iv0*3+1] + aXYZ0[iv1*3+1] + aXYZ0[iv2*3+1] + aXYZ0[iv3*3+1])*0.25;
    aXYZ1[(nv0+ne0+iq)*3+2] = (aXYZ0[iv0*3+2] + aXYZ0[iv1*3+2] + aXYZ0[iv2*3+2] + aXYZ0[iv3*3+2])*0.25;
  }
  for(int ih=0;ih<nh0;++ih){
    const int iv0 = aHex0[ih*8+0];
    const int iv1 = aHex0[ih*8+1];
    const int iv2 = aHex0[ih*8+2];
    const int iv3 = aHex0[ih*8+3];
    const int iv4 = aHex0[ih*8+4];
    const int iv5 = aHex0[ih*8+5];
    const int iv6 = aHex0[ih*8+6];
    const int iv7 = aHex0[ih*8+7];
    aXYZ1[(nv0+ne0+nq0+ih)*3+0] = (aXYZ0[iv0*3+0]+aXYZ0[iv1*3+0]+aXYZ0[iv2*3+0]+aXYZ0[iv3*3+0]+aXYZ0[iv4*3+0]+aXYZ0[iv5*3+0]+aXYZ0[iv6*3+0]+aXYZ0[iv7*3+0])*0.125;
    aXYZ1[(nv0+ne0+nq0+ih)*3+1] = (aXYZ0[iv0*3+1]+aXYZ0[iv1*3+1]+aXYZ0[iv2*3+1]+aXYZ0[iv3*3+1]+aXYZ0[iv4*3+1]+aXYZ0[iv5*3+1]+aXYZ0[iv6*3+1]+aXYZ0[iv7*3+1])*0.125;
    aXYZ1[(nv0+ne0+nq0+ih)*3+2] = (aXYZ0[iv0*3+2]+aXYZ0[iv1*3+2]+aXYZ0[iv2*3+2]+aXYZ0[iv3*3+2]+aXYZ0[iv4*3+2]+aXYZ0[iv5*3+2]+aXYZ0[iv6*3+2]+aXYZ0[iv7*3+2])*0.125;
  }
}






void CenterOfGravity_Tet
(double& v_tot,
 double& cgx, double& cgy, double& cgz,
 const std::vector<double>& aXYZC,
 const std::vector<int>& aTetC)
{
  cgx = 0.0;
  cgy = 0.0;
  cgz = 0.0;
  v_tot = 0;
  const double* pXYZ = aXYZC.data();
  for(unsigned int it=0;it<aTetC.size()/4;++it){
    const double* p0 = pXYZ+aTetC[it*4+0]*3;
    const double* p1 = pXYZ+aTetC[it*4+1]*3;
    const double* p2 = pXYZ+aTetC[it*4+2]*3;
    const double* p3 = pXYZ+aTetC[it*4+3]*3;
    double v = TetVolume3D(p0, p1, p2, p3);
    v_tot += v;
    cgx += v*(p0[0]+p1[0]+p2[0]+p3[0])*0.25;
    cgy += v*(p0[1]+p1[1]+p2[1]+p3[1])*0.25;
    cgz += v*(p0[2]+p1[2]+p2[2]+p3[2])*0.25;
  }
  cgx /= v_tot;
  cgy /= v_tot;
  cgz /= v_tot;
}
