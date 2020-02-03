/*
 * Copyright (c) 2019 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */

#include <cassert>
#include <cmath>
#include <vector>
#include "delfem2/mshmisc.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

namespace dfm2 = delfem2;

// ------------------------------------------------

static double Length3(const double p[3]){
  return sqrt(p[0]*p[0] + p[1]*p[1] + p[2]*p[2]);
}

static double Length3(const float p[3]){
  return sqrtf(p[0]*p[0] + p[1]*p[1] + p[2]*p[2]);
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

template <typename T>
inline static void UnitNormalAreaTri3
 (T n[3], T& a,
  const T v1[3], const T v2[3], const T v3[3])
{
  n[0] = ( v2[1] - v1[1] )*( v3[2] - v1[2] ) - ( v3[1] - v1[1] )*( v2[2] - v1[2] );
  n[1] = ( v2[2] - v1[2] )*( v3[0] - v1[0] ) - ( v3[2] - v1[2] )*( v2[0] - v1[0] );
  n[2] = ( v2[0] - v1[0] )*( v3[1] - v1[1] ) - ( v3[0] - v1[0] )*( v2[1] - v1[1] );
  a = sqrt(n[0]*n[0]+n[1]*n[1]+n[2]*n[2])*0.5;
  const T invlen = 0.5/a;
  n[0]*=invlen;	n[1]*=invlen;	n[2]*=invlen;
}


template <typename T>
static void MatVec3(T y[3],
                    const T m[9], const T x[3]){
  y[0] = m[0]*x[0] + m[1]*x[1] + m[2]*x[2];
  y[1] = m[3]*x[0] + m[4]*x[1] + m[5]*x[2];
  y[2] = m[6]*x[0] + m[7]*x[1] + m[8]*x[2];
}

static double Distance3(const double p0[3], const double p1[3]) {
  return sqrt( (p1[0]-p0[0])*(p1[0]-p0[0]) + (p1[1]-p0[1])*(p1[1]-p0[1]) + (p1[2]-p0[2])*(p1[2]-p0[2]) );
}

static float Distance3(const float p0[3], const float p1[3]) {
  return sqrtf( (p1[0]-p0[0])*(p1[0]-p0[0]) + (p1[1]-p0[1])*(p1[1]-p0[1]) + (p1[2]-p0[2])*(p1[2]-p0[2]) );
}


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

template <typename T>
static T TetVolume3D
 (const T v1[3],
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

static void Mat3_Bryant(double m[9],
                        double rx, double ry, double rz)
{
  m[0] = cos(rz)*cos(ry);
  m[1] = cos(rz)*sin(ry)*sin(rx)-sin(rz)*cos(rx);
  m[2] = cos(rz)*sin(ry)*cos(rx)+sin(rz)*sin(rx);
  m[3] = sin(rz)*cos(ry);
  m[4] = sin(rz)*sin(ry)*sin(rx)+cos(rz)*cos(rx);
  m[5] = sin(rz)*sin(ry)*cos(rx)-cos(rz)*sin(rx);
  m[6] = -sin(ry);
  m[7] = cos(ry)*sin(rx);
  m[8] = cos(ry)*cos(rx);
}

static void Mat3_Bryant(float m[9],
                        float rx, float ry, float rz)
{
  m[0] = cosf(rz)*cosf(ry);
  m[1] = cosf(rz)*sinf(ry)*sinf(rx)-sinf(rz)*cosf(rx);
  m[2] = cosf(rz)*sinf(ry)*cosf(rx)+sinf(rz)*sinf(rx);
  m[3] = sinf(rz)*cosf(ry);
  m[4] = sinf(rz)*sinf(ry)*sinf(rx)+cosf(rz)*cosf(rx);
  m[5] = sinf(rz)*sinf(ry)*cosf(rx)-cosf(rz)*sinf(rx);
  m[6] = -sinf(ry);
  m[7] = cosf(ry)*sinf(rx);
  m[8] = cosf(ry)*cosf(rx);
}

// static function above
// -----------------------------------------------------------------
// exposed function below

template <typename T>
void CenterWidth_MinMaxXYZ
(T& cx, T& cy, T& cz,
 T& wx, T& wy, T& wz,
 //
 T x_min, T x_max,
 T y_min, T y_max,
 T z_min, T z_max)
{
  cx = (x_min+x_max)*0.5;
  cy = (y_min+y_max)*0.5;
  cz = (z_min+z_max)*0.5;
  wx = x_max-x_min;
  wy = y_max-y_min;
  wz = z_max-z_min;
}

// -----------------------------------------------------------------------------

template<typename T>
void dfm2::updateMinMaxXYZ(
    T& x_min, T& x_max,
    T& y_min, T& y_max,
    T& z_min, T& z_max,
    T x, T y, T z)
{
  if( x_min > x_max ){
    x_min = x_max = x;
    y_min = y_max = y;
    z_min = z_max = z;
    return;
  }
  x_min = (x_min < x) ? x_min : x;
  x_max = (x_max > x) ? x_max : x;
  y_min = (y_min < y) ? y_min : y;
  y_max = (y_max > y) ? y_max : y;
  z_min = (z_min < z) ? z_min : z;
  z_max = (z_max > z) ? z_max : z;
}
template void dfm2::updateMinMaxXYZ(float& x_min, float& x_max,
                                    float& y_min, float& y_max,
                                    float& z_min, float& z_max,
                                    float X, float Y, float Z);
template void dfm2::updateMinMaxXYZ(double& x_min, double& x_max,
                                    double& y_min, double& y_max,
                                    double& z_min, double& z_max,
                                    double X, double Y, double Z);

// -----------------------------

template<typename T>
void dfm2::Min3Max3_Points3(
    T min3[3],
    T max3[3],
    const T* aXYZ,
    const unsigned int nXYZ)
{
  min3[0] = +1;
  max3[0] = -1;
  for(unsigned int ixyz=0;ixyz<nXYZ;++ixyz){
    updateMinMaxXYZ(min3[0], max3[0], min3[1], max3[1], min3[2], max3[2],
                    aXYZ[ixyz*3+0], aXYZ[ixyz*3+1], aXYZ[ixyz*3+2]);
  }
}
template void dfm2::Min3Max3_Points3(double min3[3], double max3[3],
                                     const double* aXYZ, const unsigned int nXYZ);
template void dfm2::Min3Max3_Points3(float min3[3], float max3[3],
                                     const float* aXYZ, const unsigned int nXYZ);

// --------------------------------------------------------------------------------

template <typename T>
void dfm2::CenterWidth_Point3
(T& cx, T& cy, T& cz,
 T& wx, T& wy, T& wz,
 const T* paXYZ, const unsigned int nXYZ)
{
  if( paXYZ == 0 ){ cx=cy=cz=0; wx=wy=wz=1; return; }
  T x_min=paXYZ[0], x_max=paXYZ[0];
  T y_min=paXYZ[1], y_max=paXYZ[1];
  T z_min=paXYZ[2], z_max=paXYZ[2];
  for(unsigned int ino=0;ino<nXYZ;ino++){
    updateMinMaxXYZ(x_min,x_max, y_min,y_max, z_min,z_max,
                    paXYZ[ino*3+0], paXYZ[ino*3+1], paXYZ[ino*3+2]);
  }
  CenterWidth_MinMaxXYZ(cx,cy,cz, wx,wy,wz,
                           x_min,x_max, y_min,y_max, z_min,z_max);
}
template void dfm2::CenterWidth_Point3 (float& cx, float& cy, float& cz,
                                        float& wx, float& wy, float& wz,
                                        const float* paXYZ, const unsigned int nXYZ);
template void dfm2::CenterWidth_Point3 (double& cx, double& cy, double& cz,
                                        double& wx, double& wy, double& wz,
                                        const double* paXYZ, const unsigned int nXYZ);


// ---------------------------------

template <typename T>
void dfm2::CenterWidth_Points3
(T& cx, T& cy, T& cz,
 T& wx, T& wy, T& wz,
 const std::vector<T>& aXYZ)
{
  const int np = (int)aXYZ.size()/3;
  if(np==0){ cx=cy=cz=0; wx=wy=wz=1; return; }
  T x_min=aXYZ[0], x_max=aXYZ[0];
  T y_min=aXYZ[1], y_max=aXYZ[1];
  T z_min=aXYZ[2], z_max=aXYZ[2];
  for (int ip=0; ip<np; ++ip){
    updateMinMaxXYZ(x_min,x_max, y_min,y_max, z_min,z_max,
                    aXYZ[ip*3+0], aXYZ[ip*3+1], aXYZ[ip*3+2]);
  }
  CenterWidth_MinMaxXYZ(cx,cy,cz, wx,wy,wz,
                        x_min,x_max, y_min,y_max, z_min,z_max);
}
template void dfm2::CenterWidth_Points3 (float& cx, float& cy, float& cz,
                                         float& wx, float& wy, float& wz,
                                         const std::vector<float>& aXYZ);
template void dfm2::CenterWidth_Points3 (double& cx, double& cy, double& cz,
                                         double& wx, double& wy, double& wz,
                                         const std::vector<double>& aXYZ);

// -------------------------------------

template <typename T>
void dfm2::CenterWidth_Points3(T c[3],
                               T w[3],
                               const std::vector<T>& aXYZ)
{
  dfm2::CenterWidth_Points3(c[0],c[1],c[2],
                            w[0],w[1],w[2],
                            aXYZ);
}
template void dfm2::CenterWidth_Points3(float c[3],
                                        float w[3],
                                        const std::vector<float>& aXYZ);
template void dfm2::CenterWidth_Points3(double c[3],
                                        double w[3],
                                        const std::vector<double>& aXYZ);


// -------------------------------------

void dfm2::GetCenterWidthGroup
(double& cx, double& cy, double& cz,
 double& wx, double& wy, double& wz,
 // ----------
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
  CenterWidth_MinMaxXYZ(cx,cy,cz, wx,wy,wz,
                           x_min,x_max, y_min,y_max, z_min,z_max);
}

void dfm2::GetCenterWidthGroup
 (double& cx, double& cy, double& cz,
  double& wx, double& wy, double& wz,
  //
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
  CenterWidth_MinMaxXYZ(cx,cy,cz, wx,wy,wz,
                           x_min,x_max, y_min,y_max, z_min,z_max);
}

void dfm2::GetCenterWidth3DGroup
 (double cw[6],
  //
  const std::vector<double>& aXYZ,
  const std::vector<int>& aElemInd,
  const std::vector<int>& aElem,
  int igroup,
  const std::vector<int>& aIndGroup)
{
  GetCenterWidthGroup(cw[0],cw[1],cw[2], cw[3],cw[4],cw[5],
                      aXYZ,aElemInd,aElem, igroup, aIndGroup);
}


void dfm2::GetCenterWidthLocal(
    double& lcx, double& lcy, double& lcz,
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
  CenterWidth_MinMaxXYZ(lcx,lcy,lcz, lwx,lwy,lwz,
                           x_min,x_max, y_min,y_max, z_min,z_max);
}

// -------------------------------------

template <typename T>
T dfm2::CentsMaxRad_MeshTri3(
    std::vector<T>& aXYZ_c0,
    const std::vector<T>& aXYZ,
    const std::vector<unsigned int>& aTri)
{
  T max_rad0 = -1;
  const unsigned int nTri = aTri.size()/3;
  aXYZ_c0.resize(nTri*3);
  for(std::size_t itri=0;itri<nTri;++itri) {
    const unsigned int i0 = aTri[itri*3+0];
    const unsigned int i1 = aTri[itri*3+1];
    const unsigned int i2 = aTri[itri*3+2];
    const T pc[3] = {
        (aXYZ[i0*3+0] + aXYZ[i1*3+0] + aXYZ[i2*3+0])/3,
        (aXYZ[i0*3+1] + aXYZ[i1*3+1] + aXYZ[i2*3+1])/3,
        (aXYZ[i0*3+2] + aXYZ[i1*3+2] + aXYZ[i2*3+2])/3 };
    aXYZ_c0[itri*3+0] = pc[0];
    aXYZ_c0[itri*3+1] = pc[1];
    aXYZ_c0[itri*3+2] = pc[2];
    const T l0 = Distance3(pc,aXYZ.data()+i0*3);
    const T l1 = Distance3(pc,aXYZ.data()+i1*3);
    const T l2 = Distance3(pc,aXYZ.data()+i2*3);
    if( max_rad0 < 0 || l0 > max_rad0 ){ max_rad0 = l0; }
    if( max_rad0 < 0 || l1 > max_rad0 ){ max_rad0 = l1; }
    if( max_rad0 < 0 || l2 > max_rad0 ){ max_rad0 = l2; }
  }
  return max_rad0;
}
template float dfm2::CentsMaxRad_MeshTri3(std::vector<float>& aXYZ_c0,
                                          const std::vector<float>& aXYZ,
                                          const std::vector<unsigned int>& aTri);
template double dfm2::CentsMaxRad_MeshTri3(std::vector<double>& aXYZ_c0,
                                           const std::vector<double>& aXYZ,
                                           const std::vector<unsigned int>& aTri);

// -------------------------------------

template <typename T>
void dfm2::Scale_PointsX
(std::vector<T>& aXYZ,
 T s)
{
  const unsigned int n = aXYZ.size();
  for (unsigned int i = 0; i<n; ++i){ aXYZ[i] *= s; }
}
template void dfm2::Scale_PointsX(std::vector<float>& aXYZ, float s);
template void dfm2::Scale_PointsX(std::vector<double>& aXYZ, double s);

// --------------

template <typename T>
void dfm2::Scale_Points3
 (T* pXYZs_,
  const unsigned int nnode_,
  T s)
{
  for(unsigned int ino=0;ino<nnode_;ino++){
    pXYZs_[ino*3+0] *= s;
    pXYZs_[ino*3+1] *= s;
    pXYZs_[ino*3+2] *= s;
  }
}
template void dfm2::Scale_Points3(float* pXYZ, const unsigned int n, float s);
template void dfm2::Scale_Points3(double* pXYZ, const unsigned int n, double s);

// --------------

template <typename T>
void dfm2::Translate_Points3
(std::vector<T>& aXYZ,
T tx, T ty, T tz)
{
  const unsigned int nXYZ = aXYZ.size()/3;
  for (unsigned int ixyz = 0; ixyz<nXYZ; ixyz++){
    aXYZ[ixyz*3+0] += tx;
    aXYZ[ixyz*3+1] += ty;
    aXYZ[ixyz*3+2] += tz;
  }
}
template void dfm2::Translate_Points3(std::vector<float>& aXYZ, float tx, float ty, float tz);
template void dfm2::Translate_Points3(std::vector<double>& aXYZ, double tx, double ty, double tz);

// --------------

template <typename T>
void dfm2::Translate_Points3
(T* pXYZs_,
 const unsigned int nnode_,
 // ----
 T tx, T ty, T tz)
{
  for(unsigned int ino=0;ino<nnode_;ino++){
    pXYZs_[ino*3+0] += tx;
    pXYZs_[ino*3+1] += ty;
    pXYZs_[ino*3+2] += tz;
  }
}
template void dfm2::Translate_Points3(float* pXYZ, unsigned int nN, float tx, float ty, float tz);
template void dfm2::Translate_Points3(double* pXYZ, unsigned int nN, double tx, double ty, double tz);

// --------------

template <typename T>
void dfm2::Translate_Points2
 (std::vector<T>& aXY,
  T tx, T ty)
{
  const unsigned int np = aXY.size()/2;
  for (unsigned int ip = 0; ip<np; ip++){
    aXY[ip*2+0] += tx;
    aXY[ip*2+1] += ty;
  }
}
template void dfm2::Translate_Points2(std::vector<float>& aXYZ, float tx, float ty);
template void dfm2::Translate_Points2(std::vector<double>& aXYZ, double tx, double ty);

// --------------

template <typename T>
void dfm2::Rotate_Points3
(std::vector<T>& aXYZ,
T radx, T rady, T radz)
{
  T mat[9]; Mat3_Bryant(mat, radx, rady, radz);
  T* pXYZ = aXYZ.data();
  const unsigned int nXYZ = aXYZ.size()/3;
  for (unsigned int ixyz = 0; ixyz<nXYZ; ++ixyz){
    const T p[3] = { aXYZ[ixyz*3+0], aXYZ[ixyz*3+1], aXYZ[ixyz*3+2] };
    MatVec3(pXYZ+ixyz*3,  mat, p);
  }
}
template void dfm2::Rotate_Points3 (std::vector<float>& aXYZ,
                                    float radx, float rady, float radz);
template void dfm2::Rotate_Points3 (std::vector<double>& aXYZ,
                                    double radx, double rady, double radz);

// -----------------------------------------

double dfm2::Size_Points3D_LongestAABBEdge
 (const std::vector<double>& aXYZ)
{
  double c[3], w[3];
  CenterWidth_Points3(c, w,
                      aXYZ);
  return largest(w[0], w[1], w[2]);
}

void dfm2::Normalize_Points3D
(std::vector<double>& aXYZ,
 double s)
{
  double c[3], w[3];
  CenterWidth_Points3(c,w,
                      aXYZ);
  Translate_Points3(aXYZ,
                    -c[0], -c[1], -c[2]);
  double wmax = largest(w[0], w[1], w[2]);
  Scale_PointsX(aXYZ,
                s/wmax);
}



// ---------------------------------------

template <typename T>
void dfm2::CG_Point3
 (T cg[3],
  const std::vector<T>& aXYZ)
{
  cg[0] = cg[1] = cg[2] = 0;
  unsigned int nXYZ = aXYZ.size()/3;
  for (unsigned int ixyz = 0; ixyz<nXYZ; ixyz++){
    cg[0] += aXYZ[ixyz*3+0];
    cg[1] += aXYZ[ixyz*3+1];
    cg[2] += aXYZ[ixyz*3+2];
  }
  cg[0] /= nXYZ;
  cg[1] /= nXYZ;
  cg[2] /= nXYZ;
}
template void dfm2::CG_Point3(float cg[3], const std::vector<float>& aXYZ);
template void dfm2::CG_Point3(double cg[3], const std::vector<double>& aXYZ);

// above: points3
// --------------------------------------------------------------------------------------------------------------------
// below: mesh


void dfm2::RemoveUnreferencedPoints_MeshElem
 (std::vector<double>& aXYZ1,
  std::vector<unsigned int>& aElem1,
  std::vector<int>& aMap01,
  unsigned int ndim,
  const std::vector<double>& aXYZ0,
  const std::vector<unsigned int>& aElem0)
{
  unsigned int np0 = aXYZ0.size()/ndim;
  aMap01.assign(np0,-2);
  for(int ip : aElem0){
    aMap01[ip] = -1;
  }
  int npj = 0;
  for(unsigned int ip=0;ip<np0;++ip){
    if( aMap01[ip] == -2 ) continue;
    aMap01[ip] = npj;
    npj++;
  }
  aXYZ1.resize(npj*ndim);
  for(unsigned int ip=0;ip<np0;++ip){
    if( aMap01[ip] == -2 ) continue;
    int jp = aMap01[ip];
    for(unsigned int idim=0;idim<ndim;++idim){
      aXYZ1[jp*ndim+idim] = aXYZ0[ip*ndim+idim];
    }
  }
  aElem1.resize(aElem0.size());
  for(std::size_t it=0;it<aElem0.size();++it){
    int ip = aElem0[it];
    int jp = aMap01[ip];
    aElem1[it] = jp;
  }
}

void dfm2::Normal_MeshTri3D
(double* aNorm,
 const double* aXYZ,
 unsigned int nXYZ,
 const unsigned int* aTri,
 unsigned int nTri)
{
  for(unsigned int i=0;i<nXYZ*3;i++){ aNorm[i] = 0; }
  for(unsigned int itri=0;itri<nTri;itri++){
    const unsigned int i0 = aTri[itri*3+0];
    const unsigned int i1 = aTri[itri*3+1];
    const unsigned int i2 = aTri[itri*3+2];
    const double* p0 = aXYZ+i0*3;
    const double* p1 = aXYZ+i1*3;
    const double* p2 = aXYZ+i2*3;
    double un[3], area;
    UnitNormalAreaTri3(un,area, p0,p1,p2);
    aNorm[i0*3+0] += un[0];  aNorm[i0*3+1] += un[1];  aNorm[i0*3+2] += un[2];
    aNorm[i1*3+0] += un[0];  aNorm[i1*3+1] += un[1];  aNorm[i1*3+2] += un[2];
    aNorm[i2*3+0] += un[0];  aNorm[i2*3+1] += un[1];  aNorm[i2*3+2] += un[2];
  }
  for(unsigned int ino=0;ino<nXYZ;ino++){
    const double n[3] = {aNorm[ino*3+0],aNorm[ino*3+1],aNorm[ino*3+2]};
    const double invlen = 1.0/Length3(n);
    aNorm[ino*3+0] *= invlen;
    aNorm[ino*3+1] *= invlen;
    aNorm[ino*3+2] *= invlen;
  }
}

template <typename REAL>
void dfm2::Normal_MeshQuad3
 (std::vector<REAL>& aNorm,
  const std::vector<REAL>& aXYZ,
  const std::vector<unsigned int>& aQuad)
{
  const unsigned int nXYZ = aXYZ.size()/3;
  const unsigned int nQuad = aQuad.size()/4;
  aNorm.resize(nXYZ*3);
  // -------
  for(unsigned int i=0;i<nXYZ*3;i++){ aNorm[i] = 0; }
  for(unsigned int iquad=0;iquad<nQuad;++iquad){
    const unsigned int i0 = aQuad[iquad*4+0];
    const unsigned int i1 = aQuad[iquad*4+1];
    const unsigned int i2 = aQuad[iquad*4+2];
    const unsigned int i3 = aQuad[iquad*4+3];
    const REAL* p0 = aXYZ.data()+i0*3;
    const REAL* p1 = aXYZ.data()+i1*3;
    const REAL* p2 = aXYZ.data()+i2*3;
    const REAL* p3 = aXYZ.data()+i3*3;
    REAL un0[3], a0; UnitNormalAreaTri3(un0,a0, p3,p0,p1);
    REAL un1[3], a1; UnitNormalAreaTri3(un1,a1, p0,p1,p2);
    REAL un2[3], a2; UnitNormalAreaTri3(un2,a2, p1,p2,p3);
    REAL un3[3], a3; UnitNormalAreaTri3(un3,a3, p2,p3,p0);
    aNorm[i0*3+0] += un0[0];  aNorm[i0*3+1] += un0[1];  aNorm[i0*3+2] += un0[2];
    aNorm[i1*3+0] += un1[0];  aNorm[i1*3+1] += un1[1];  aNorm[i1*3+2] += un1[2];
    aNorm[i2*3+0] += un2[0];  aNorm[i2*3+1] += un2[1];  aNorm[i2*3+2] += un2[2];
    aNorm[i3*3+0] += un3[0];  aNorm[i3*3+1] += un3[1];  aNorm[i3*3+2] += un3[2];
  }
  for(unsigned int ino=0;ino<nXYZ;ino++){
    const REAL n[3] = {aNorm[ino*3+0],aNorm[ino*3+1],aNorm[ino*3+2]};
    const REAL invlen = 1.0/Length3(n);
    aNorm[ino*3+0] *= invlen;
    aNorm[ino*3+1] *= invlen;
    aNorm[ino*3+2] *= invlen;
  }
}
template void dfm2::Normal_MeshQuad3(std::vector<float>& aNorm,
                                     const std::vector<float>& aXYZ,
                                     const std::vector<unsigned int>& aQuad);
template void dfm2::Normal_MeshQuad3(std::vector<double>& aNorm,
                                     const std::vector<double>& aXYZ,
                                     const std::vector<unsigned int>& aQuad);



void dfm2::Quality_MeshTri2D
 (double& max_aspect, double& min_area,
  const double* aXY,
  const unsigned int* aTri, unsigned int nTri)
{
  max_aspect = 0;
  min_area = 0;
  for(unsigned int itri=0;itri<nTri;itri++){
    const unsigned int i0 = aTri[itri*3+0];
    const unsigned int i1 = aTri[itri*3+1];
    const unsigned int i2 = aTri[itri*3+2];
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


template <typename T>
void dfm2::CG_MeshTri3_Solid(
    T cg[3],
    const std::vector<T>& aXYZ,
    const std::vector<unsigned int>& aTri)
{ // center of gravity
  cg[0] = cg[1] = cg[2] = 0.0;
  double tw = 0;
  const unsigned int nTri = aTri.size()/3;
  for (std::size_t itri = 0; itri<nTri; itri++){
    unsigned int i1 = aTri[itri*3+0];
    unsigned int i2 = aTri[itri*3+1];
    unsigned int i3 = aTri[itri*3+2];
    const T q0[3] = { 0, 0, 0 };
    const T q1[3] = { aXYZ[i1*3+0], aXYZ[i1*3+1], aXYZ[i1*3+2] };
    const T q2[3] = { aXYZ[i2*3+0], aXYZ[i2*3+1], aXYZ[i2*3+2] };
    const T q3[3] = { aXYZ[i3*3+0], aXYZ[i3*3+1], aXYZ[i3*3+2] };
    T v = ::TetVolume3D(q0, q1, q2, q3);
    tw += v;
    cg[0] += (q0[0]+q1[0]+q2[0]+q3[0])*0.25*v;
    cg[1] += (q0[1]+q1[1]+q2[1]+q3[1])*0.25*v;
    cg[2] += (q0[2]+q1[2]+q2[2]+q3[2])*0.25*v;
  }
  cg[0] /= tw;
  cg[1] /= tw;
  cg[2] /= tw;
}
template void dfm2::CG_MeshTri3_Solid(float cg[3],
                                      const std::vector<float>& aXYZ,
                                      const std::vector<unsigned int>& aTri);
template void dfm2::CG_MeshTri3_Solid(double cg[3],
                                      const std::vector<double>& aXYZ,
                                      const std::vector<unsigned int>& aTri);

// ----------------------------------------

template <typename T>
void dfm2::CG_MeshTri3_Shell
(T cg[3],
 const std::vector<T>& aXYZ,
 const std::vector<unsigned int>& aTri)
{ // center of gravity
  cg[0] = cg[1] = cg[2] = 0.0;
  double tw = 0;
  const unsigned int nTri = aTri.size()/3;
  for (std::size_t itri = 0; itri<nTri; itri++){
    unsigned int i1 = aTri[itri*3+0];
    unsigned int i2 = aTri[itri*3+1];
    unsigned int i3 = aTri[itri*3+2];
    const double q1[3] = { aXYZ[i1*3+0], aXYZ[i1*3+1], aXYZ[i1*3+2] };
    const double q2[3] = { aXYZ[i2*3+0], aXYZ[i2*3+1], aXYZ[i2*3+2] };
    const double q3[3] = { aXYZ[i3*3+0], aXYZ[i3*3+1], aXYZ[i3*3+2] };
    double a = TriArea3D(q1, q2, q3);
    tw += a;
    cg[0] += (q1[0]+q2[0]+q3[0])*0.333333*a;
    cg[1] += (q1[1]+q2[1]+q3[1])*0.333333*a;
    cg[2] += (q1[2]+q2[2]+q3[2])*0.333333*a;
  }
  cg[0] /= tw;
  cg[1] /= tw;
  cg[2] /= tw;
}
template void dfm2::CG_MeshTri3_Shell(float cg[3],
                                      const std::vector<float>& aXYZ,
                                      const std::vector<unsigned int>& aTri);
template void dfm2::CG_MeshTri3_Shell(double cg[3],
                                      const std::vector<double>& aXYZ,
                                      const std::vector<unsigned int>& aTri);

// ------------------------------------------

template <typename T>
T dfm2::CG_TriMsh3Flg_Shell
(T cg[3],
 const std::vector<T>& aXYZ,
 const std::vector<unsigned int>& aTri,
 int iflg,
 const std::vector<int>& aFlg)
{
  cg[0] = cg[1] = cg[2] = 0.0;
  double tw = 0;
  const std::size_t nTri = aTri.size()/3;
  for (std::size_t itri = 0; itri<nTri; itri++){
    if( aFlg[itri] != iflg ) continue;
    const unsigned int i1 = aTri[itri*3+0];
    const unsigned int i2 = aTri[itri*3+1];
    const unsigned int i3 = aTri[itri*3+2];
    const double q1[3] = { aXYZ[i1*3+0], aXYZ[i1*3+1], aXYZ[i1*3+2] };
    const double q2[3] = { aXYZ[i2*3+0], aXYZ[i2*3+1], aXYZ[i2*3+2] };
    const double q3[3] = { aXYZ[i3*3+0], aXYZ[i3*3+1], aXYZ[i3*3+2] };
    double a = TriArea3D(q1, q2, q3);
    tw += a;
    cg[0] += (q1[0]+q2[0]+q3[0])*0.333333*a;
    cg[1] += (q1[1]+q2[1]+q3[1])*0.333333*a;
    cg[2] += (q1[2]+q2[2]+q3[2])*0.333333*a;
  }
  cg[0] /= tw;
  cg[1] /= tw;
  cg[2] /= tw;
  return tw;
}

template float dfm2::CG_TriMsh3Flg_Shell(float cg[3],
                                         const std::vector<float>& aXYZ,
                                         const std::vector<unsigned int>& aTri,
                                         int iflg,
                                         const std::vector<int>& aFlg);

template double dfm2::CG_TriMsh3Flg_Shell(double cg[3],
                                          const std::vector<double>& aXYZ,
                                          const std::vector<unsigned int>& aTri,
                                          int iflg,
                                          const std::vector<int>& aFlg);

// ------------------------------------------

void dfm2::CG_Tri
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

template <typename T>
void dfm2::CG_MeshTet3
 (T& v_tot,
  T cg[3],
  const std::vector<T>& aXYZC,
  const std::vector<unsigned int>& aTet)
{
  v_tot = cg[0] = cg[1] = cg[2] = 0.0;
  const T* pXYZ = aXYZC.data();
  const std::size_t nTet = aTet.size()/4;
  for(std::size_t it=0;it<nTet;++it){
    const T* p0 = pXYZ+aTet[it*4+0]*3;
    const T* p1 = pXYZ+aTet[it*4+1]*3;
    const T* p2 = pXYZ+aTet[it*4+2]*3;
    const T* p3 = pXYZ+aTet[it*4+3]*3;
    const double v = TetVolume3D(p0, p1, p2, p3);
    v_tot += v;
    cg[0] += v*(p0[0]+p1[0]+p2[0]+p3[0])*0.25;
    cg[1] += v*(p0[1]+p1[1]+p2[1]+p3[1])*0.25;
    cg[2] += v*(p0[2]+p1[2]+p2[2]+p3[2])*0.25;
  }
  cg[0] /= v_tot;
  cg[1] /= v_tot;
  cg[2] /= v_tot;
}
template void dfm2::CG_MeshTet3(float& v_tot, float cg[3], const std::vector<float>& aXYZ, const std::vector<unsigned int>& aTet);
template void dfm2::CG_MeshTet3(double& v_tot, double cg[3], const std::vector<double>& aXYZ, const std::vector<unsigned int>& aTet);

// -----------------------------------------------------------------------------------------

void dfm2::SetTopology_ExtrudeTri2Tet
(unsigned int* aTet,
 int nXY,
 const unsigned int* aTri, int nTri,
 int nlayer)
{
  for(int il=0;il<nlayer;++il){
    for(int itri=0;itri<nTri;++itri){
      unsigned int ip0=0, ip1=0, ip2=0;
      {
        const unsigned int i0 = aTri[itri*3+0];
        const unsigned int i1 = aTri[itri*3+1];
        const unsigned int i2 = aTri[itri*3+2];
        assert( i0 != i1 && i1 != i2 );
        if( i0 > i1 && i0 > i2 ){ ip0=i0; ip1=i1; ip2=i2; }
        if( i1 > i0 && i1 > i2 ){ ip0=i1; ip1=i2; ip2=i0; }
        if( i2 > i0 && i2 > i1 ){ ip0=i2; ip1=i0; ip2=i1; }
      }
      const unsigned int aIQ[6] = {
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

void dfm2::ExtrudeTri2Tet
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

// -----------------------------------------------------------------------

void dfm2::LaplacianSmoothing(
    std::vector<double>& aXYZ,
    const std::vector<int>& aTri,
    const std::vector<int>& elsup_ind,
    const std::vector<int> elsup)
{
  for(std::size_t ip=0;ip<aXYZ.size()/3;++ip){
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
  double l1 = Length3(v1);
  double l2 = Length3(v2);
  double l3 = Length3(v3);
  double crs_v1_v2[3]; Cross3D(crs_v1_v2,v1,v2);
  double den = Dot(crs_v1_v2,v3);
  double num = l1*l2*l3+(Dot(v1,v2))*l3+(Dot(v2,v3))*l1+(Dot(v3,v1))*l2;
  double tho = den/num;
  double v = atan(tho);
  if (v<0){ v += 2*M_PI; }
  v *= 2;
  return v;
}

void dfm2::makeSolidAngle
(std::vector<double>& aSolidAngle,
 const std::vector<double>& aXYZ,
 const std::vector<unsigned int>& aTri,
 const std::vector<double>& aNorm, 
 std::vector<int>& elsup_ind,
 std::vector<int>& elsup)
{
  const unsigned int nXYZ = aXYZ.size()/3;
  aSolidAngle.resize(nXYZ);
  for(unsigned int ip=0;ip<nXYZ;++ip){
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


void dfm2::MassPoint_Tet3D
(double* aMassMatrixLumped,
 double rho,
 const double* aXYZ, unsigned int nXYZ,
 const unsigned int* aTet, unsigned int nTet)
{
  for(unsigned int i=0;i<nXYZ;++i){ aMassMatrixLumped[i] = 0.0; }
  for(unsigned int it=0;it<nTet;++it){
    const unsigned int i0 = aTet[it*4+0]; assert(i0<nXYZ);
    const unsigned int i1 = aTet[it*4+1]; assert(i1<nXYZ);
    const unsigned int i2 = aTet[it*4+2]; assert(i2<nXYZ);
    const unsigned int i3 = aTet[it*4+3]; assert(i3<nXYZ);
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

void dfm2::MassPoint_Tri2D
(double* aMassMatrixLumped,
 double rho,
 const double* aXY, unsigned int nXY,
 const unsigned int* aTri, unsigned int nTri)
{
  for(unsigned int i=0;i<nXY;++i){ aMassMatrixLumped[i] = 0.0; }
  for(unsigned int it=0;it<nTri;++it){
    const unsigned int i0 = aTri[it*3+0]; assert(i0<nXY);
    const unsigned int i1 = aTri[it*3+1]; assert(i1<nXY);
    const unsigned int i2 = aTri[it*3+2]; assert(i2<nXY);
    const double* p0 = aXY+i0*2;
    const double* p1 = aXY+i1*2;
    const double* p2 = aXY+i2*2;
    const double a012 = TriArea2D(p0, p1, p2);
    aMassMatrixLumped[i0] += rho*a012/3.0;
    aMassMatrixLumped[i1] += rho*a012/3.0;
    aMassMatrixLumped[i2] += rho*a012/3.0;
  }
}

// TODO: make this handle open surface (average face & edge independently)
void dfm2::SubdivisionPoints_QuadCatmullClark
(std::vector<double>& aXYZ1,
 // ------------------------
 const std::vector<unsigned int>& aQuad1,
 const std::vector<int>& aEdgeFace0,
 const std::vector<unsigned int> &psupIndQuad0,
 const std::vector<unsigned int> &psupQuad0,
 const unsigned int* aQuad0,
 unsigned int nQuad0,
 const double* aXYZ0,
 unsigned int nXYZ0)
{
  const unsigned int nv0 = nXYZ0;
  const unsigned int ne0 = psupQuad0.size();
  const unsigned int nq0 = nQuad0;
  assert( aEdgeFace0.size() == ne0*4 );
  aXYZ1.resize((nv0+ne0+nq0)*3);
  std::vector<int> aW(nv0,0);
  for(unsigned int iv=0;iv<nv0;++iv){
    aXYZ1[iv*3+0] = 0;
    aXYZ1[iv*3+1] = 0;
    aXYZ1[iv*3+2] = 0;
  }
  for(unsigned int iq=0;iq<nq0;++iq){
    const unsigned int iv0 = aQuad0[iq*4+0];
    const unsigned int iv1 = aQuad0[iq*4+1];
    const unsigned int iv2 = aQuad0[iq*4+2];
    const unsigned int iv3 = aQuad0[iq*4+3];
    double p0x = (aXYZ0[iv0*3+0] + aXYZ0[iv1*3+0] + aXYZ0[iv2*3+0] + aXYZ0[iv3*3+0])*0.25;
    double p0y = (aXYZ0[iv0*3+1] + aXYZ0[iv1*3+1] + aXYZ0[iv2*3+1] + aXYZ0[iv3*3+1])*0.25;
    double p0z = (aXYZ0[iv0*3+2] + aXYZ0[iv1*3+2] + aXYZ0[iv2*3+2] + aXYZ0[iv3*3+2])*0.25;
    aXYZ1[(nv0+ne0+iq)*3+0] = p0x;
    aXYZ1[(nv0+ne0+iq)*3+1] = p0y;
    aXYZ1[(nv0+ne0+iq)*3+2] = p0z;
    const unsigned int aIV[4] = { iv0, iv1, iv2, iv3 };
    for(unsigned int jv0 : aIV){
      aXYZ1[jv0*3+0] += p0x;
      aXYZ1[jv0*3+1] += p0y;
      aXYZ1[jv0*3+2] += p0z;
      aW[jv0] += 1;
    }
  }
  for(unsigned int ie=0;ie<ne0;++ie){
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
  for(unsigned int iv=0;iv<nv0;++iv){
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

void dfm2::SubdivPoints3_MeshQuad
(std::vector<double>& aXYZ1,
 // ------------
 const std::vector<int>& aEdgeFace0,
 const std::vector<unsigned int>& aQuad0,
 const std::vector<double>& aXYZ0)
{
  const unsigned int nv0 = aXYZ0.size()/3;
  const unsigned int ne0 = aEdgeFace0.size()/4;
  const unsigned int nq0 = aQuad0.size()/4;
  assert( aEdgeFace0.size() == ne0*4 );
  aXYZ1.resize((nv0+ne0+nq0)*3);
  for(unsigned int iv=0;iv<nv0;++iv){
    aXYZ1[iv*3+0] = aXYZ0[iv*3+0];
    aXYZ1[iv*3+1] = aXYZ0[iv*3+1];
    aXYZ1[iv*3+2] = aXYZ0[iv*3+2];
  }
  for(unsigned int ie=0;ie<ne0;++ie){
    const int iv0 = aEdgeFace0[ie*4+0];
    const int iv1 = aEdgeFace0[ie*4+1];
    aXYZ1[(nv0+ie)*3+0] = (aXYZ0[iv0*3+0] + aXYZ0[iv1*3+0])*0.5;
    aXYZ1[(nv0+ie)*3+1] = (aXYZ0[iv0*3+1] + aXYZ0[iv1*3+1])*0.5;
    aXYZ1[(nv0+ie)*3+2] = (aXYZ0[iv0*3+2] + aXYZ0[iv1*3+2])*0.5;
  }
  for(unsigned int iq=0;iq<nq0;++iq){
    const int iv0 = aQuad0[iq*4+0];
    const int iv1 = aQuad0[iq*4+1];
    const int iv2 = aQuad0[iq*4+2];
    const int iv3 = aQuad0[iq*4+3];
    aXYZ1[(nv0+ne0+iq)*3+0] = (aXYZ0[iv0*3+0] + aXYZ0[iv1*3+0] + aXYZ0[iv2*3+0] + aXYZ0[iv3*3+0])*0.25;
    aXYZ1[(nv0+ne0+iq)*3+1] = (aXYZ0[iv0*3+1] + aXYZ0[iv1*3+1] + aXYZ0[iv2*3+1] + aXYZ0[iv3*3+1])*0.25;
    aXYZ1[(nv0+ne0+iq)*3+2] = (aXYZ0[iv0*3+2] + aXYZ0[iv1*3+2] + aXYZ0[iv2*3+2] + aXYZ0[iv3*3+2])*0.25;
  }
}


void dfm2::SubdivisionPoints_Hex
(std::vector<double>& aXYZ1,
 // -------------
 const std::vector<unsigned int> &psupIndHex0,
 const std::vector<unsigned int> &psupHex0,
 const std::vector<unsigned int>& aQuadHex0,
 const unsigned int* aHex0, unsigned int nHex0,
 const double* aXYZ0, unsigned int nXYZ0)
{
  const unsigned int nv0 = nXYZ0;
  const unsigned int ne0 = psupHex0.size();
  const unsigned int nq0 = aQuadHex0.size()/4;
  const unsigned int nh0 = nHex0;
  aXYZ1.resize((nv0+ne0+nq0+nh0)*3);
  for(unsigned int iv=0;iv<nv0;++iv){
    aXYZ1[iv*3+0] = aXYZ0[iv*3+0];
    aXYZ1[iv*3+1] = aXYZ0[iv*3+1];
    aXYZ1[iv*3+2] = aXYZ0[iv*3+2];
  }
  for(unsigned int iv=0;iv<nv0;++iv){
    for(unsigned int ipsup=psupIndHex0[iv];ipsup<psupIndHex0[iv+1];++ipsup){
      int jv = psupHex0[ipsup];
      aXYZ1[(nv0+ipsup)*3+0] = (aXYZ0[iv*3+0] + aXYZ0[jv*3+0])*0.5;
      aXYZ1[(nv0+ipsup)*3+1] = (aXYZ0[iv*3+1] + aXYZ0[jv*3+1])*0.5;
      aXYZ1[(nv0+ipsup)*3+2] = (aXYZ0[iv*3+2] + aXYZ0[jv*3+2])*0.5;
    }
  }
  for(unsigned int iq=0;iq<nq0;++iq){
    const int iv0 = aQuadHex0[iq*4+0];
    const int iv1 = aQuadHex0[iq*4+1];
    const int iv2 = aQuadHex0[iq*4+2];
    const int iv3 = aQuadHex0[iq*4+3];
    aXYZ1[(nv0+ne0+iq)*3+0] = (aXYZ0[iv0*3+0] + aXYZ0[iv1*3+0] + aXYZ0[iv2*3+0] + aXYZ0[iv3*3+0])*0.25;
    aXYZ1[(nv0+ne0+iq)*3+1] = (aXYZ0[iv0*3+1] + aXYZ0[iv1*3+1] + aXYZ0[iv2*3+1] + aXYZ0[iv3*3+1])*0.25;
    aXYZ1[(nv0+ne0+iq)*3+2] = (aXYZ0[iv0*3+2] + aXYZ0[iv1*3+2] + aXYZ0[iv2*3+2] + aXYZ0[iv3*3+2])*0.25;
  }
  for(unsigned int ih=0;ih<nh0;++ih){
    const unsigned int iv0 = aHex0[ih*8+0];
    const unsigned int iv1 = aHex0[ih*8+1];
    const unsigned int iv2 = aHex0[ih*8+2];
    const unsigned int iv3 = aHex0[ih*8+3];
    const unsigned int iv4 = aHex0[ih*8+4];
    const unsigned int iv5 = aHex0[ih*8+5];
    const unsigned int iv6 = aHex0[ih*8+6];
    const unsigned int iv7 = aHex0[ih*8+7];
    aXYZ1[(nv0+ne0+nq0+ih)*3+0] = (aXYZ0[iv0*3+0]+aXYZ0[iv1*3+0]+aXYZ0[iv2*3+0]+aXYZ0[iv3*3+0]+aXYZ0[iv4*3+0]+aXYZ0[iv5*3+0]+aXYZ0[iv6*3+0]+aXYZ0[iv7*3+0])*0.125;
    aXYZ1[(nv0+ne0+nq0+ih)*3+1] = (aXYZ0[iv0*3+1]+aXYZ0[iv1*3+1]+aXYZ0[iv2*3+1]+aXYZ0[iv3*3+1]+aXYZ0[iv4*3+1]+aXYZ0[iv5*3+1]+aXYZ0[iv6*3+1]+aXYZ0[iv7*3+1])*0.125;
    aXYZ1[(nv0+ne0+nq0+ih)*3+2] = (aXYZ0[iv0*3+2]+aXYZ0[iv1*3+2]+aXYZ0[iv2*3+2]+aXYZ0[iv3*3+2]+aXYZ0[iv4*3+2]+aXYZ0[iv5*3+2]+aXYZ0[iv6*3+2]+aXYZ0[iv7*3+2])*0.125;
  }
}


