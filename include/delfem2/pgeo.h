/*
 * Copyright (c) 2019 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */

#ifndef DFM2_PGEO_H
#define DFM2_PGEO_H

#include <cassert>
#include <cstdio>
#include <vector>

#include "delfem2/dfm2_inline.h"

namespace delfem2 {

inline void FlatKnot(
    std::vector<double>& aKnotFlat,
    const std::vector<int>& aKnotMulti,
    const std::vector<double>& aKnot)
{
  assert(aKnot.size()==aKnotMulti.size());
  aKnotFlat.clear();
  for(size_t ik=0;ik<aKnot.size();++ik){
    for(int iik=0;iik<aKnotMulti[ik];++iik){
      aKnotFlat.push_back(aKnot[ik]);
    }
  }
}

 
template <typename T>
T QuadBilinear(
    int iq, double r0, double r1,
    std::vector<int>& aQuad,
    std::vector<T>& aPoint)
{
  int i0 = aQuad[iq*4+0];
  int i1 = aQuad[iq*4+1];
  int i2 = aQuad[iq*4+2];
  int i3 = aQuad[iq*4+3];
  return (1-r0)*(1-r1)*aPoint[i0] + r0*(1-r1)*aPoint[i1] + r0*r1*aPoint[i2] + (1-r0)*r1*aPoint[i3];
}

template <typename T>
T getPointCoonsQuad_CubicBezier(
    double u, double v,
    T aP[12])
{
  T p01u = getPointCubicBezierCurve(u, aP[0], aP[1], aP[2], aP[3]);
  T p32u = getPointCubicBezierCurve(u, aP[9], aP[8], aP[7], aP[6]);
  T p = (1-v)*p01u + v*p32u;
  //
  T q03v = getPointCubicBezierCurve(v, aP[0], aP[11], aP[10], aP[9]);
  T q12v = getPointCubicBezierCurve(v, aP[3], aP[ 4], aP[ 5], aP[6]);
  T q = (1-u)*q03v + u*q12v;
  //
  T r = (1-u)*(1-v)*aP[0] + u*(1-v)*aP[3] + u*v*aP[6] + (1-u)*v*aP[9];
  return p+q-r;
}
 
template <typename T>
void getCubicBezierSurface(
    const int n, // number of segment
    std::vector<T>& aP,
    const std::vector<T>& aCP)
{
  aP.resize((n+1)*(n+1));
  for (int i = 0; i<(n+1); ++i){
    for (int j = 0; j<(n+1); ++j){
      double u = (double)i/n;
      double v = (double)j/n;
      aP[i*(n+1)+j] = getPointSurfaceBezierCubic(
          u,v,
          aCP[ 0], aCP[ 1], aCP[ 2], aCP[ 3],
          aCP[ 4], aCP[ 5], aCP[ 6], aCP[ 7],
          aCP[ 8], aCP[ 9], aCP[10], aCP[11],
          aCP[12], aCP[13], aCP[14], aCP[15]);
    }
  }
}

template <typename T>
T getPointCoonsQuad_CubicBezierEdge(
    double u, double v,
    T aP[12])
{
  T p01u = getPointCubicBezierCurve(u, aP[0], aP[1], aP[2], aP[3]);
  T p32u = getPointCubicBezierCurve(u, aP[9], aP[8], aP[7], aP[6]);
  T p = (1-v)*p01u + v*p32u;
  //
  T q03v = getPointCubicBezierCurve(v, aP[0], aP[11], aP[10], aP[9]);
  T q12v = getPointCubicBezierCurve(v, aP[3], aP[ 4], aP[ 5], aP[6]);
  T q = (1-u)*q03v + u*q12v;
  T r = (1-u)*(1-v)*aP[0] + u*(1-v)*aP[3] + u*v*aP[6] + (1-u)*v*aP[9];
  return p+q-r;
}

template <typename T>
T getPointCoonsTri_CubicBezierEdge(
    double u, double v, double w,
    T aP[9])
{
  T peu = PointOnCubicBezierCurve(w/(1-u), aP[3], aP[4], aP[5], aP[6]);
  T pev = PointOnCubicBezierCurve(u/(1-v), aP[6], aP[7], aP[8], aP[0]);
  T pew = PointOnCubicBezierCurve(v/(1-w), aP[0], aP[1], aP[2], aP[3]);
  T pu = (1-u)*peu + u*aP[0];
  T pv = (1-v)*pev + v*aP[3];
  T pw = (1-w)*pew + w*aP[6];
  T pl = u*aP[0] + v*aP[3] + w*aP[6];
  return pu+pv+pw-2*pl;
}

template <typename T>
T getPointHermetianQuad(
    double u, double v,
    T aP[12])
{
  double u0 = +2*u*u*u-3*u*u+1;
  double u1 = -2*u*u*u+3*u*u;
  double du0 = +1*u*u*u-2*u*u+u;
  double du1 = +1*u*u*u-1*u*u;
  //
  double v0 = +2*v*v*v-3*v*v+1;
  double v1 = -2*v*v*v+3*v*v;
  double dv0 = +1*v*v*v-2*v*v+v;
  double dv1 = +1*v*v*v-1*v*v;
  //
  T p = aP[0]*u0*v0 + aP[3]*u1*v0 + aP[6]*u1*v1 + aP[9]*u0*v1;
  T q = 3*(aP[ 1]-aP[0])*du0*v0 + 3*(aP[3]-aP[2])*du1*v0 + 3*(aP[6]-aP[7])*du1*v1 + 3*(aP[8]-aP[ 9])*du0*v1;
  T r = 3*(aP[11]-aP[0])*u0*dv0 + 3*(aP[4]-aP[3])*u1*dv0 + 3*(aP[6]-aP[5])*u1*dv1 + 3*(aP[9]-aP[10])*u0*dv1;
  return p+q+r;
}

template <typename VEC>
VEC PointOnCubicBezierCurve(
    double t,
    const VEC& p1,
    const VEC& p2,
    const VEC& p3,
    const VEC& p4)
{
  double tp = 1.0-t;
  return t*t*t*p4
  + 3*t*t*tp*p3
  + 3*t*tp*tp*p2
  + tp*tp*tp*p1;
}

template <typename VEC>
VEC PointOnQuadraticBezierCurve(
    double t,
    const VEC &p1,
    const VEC &p2,
    const VEC &p3) {
  double tp = 1.0 - t;
  return (t * t) * p3 + (2 * t * tp) * p2 + (tp * tp) * p1;
}

template <typename T>
T getTangentCubicBezierCurve(
    double t,
    const T& p1, const T& p2, const T& p3, const T& p4)
{
  double tp = 1.0-t;
  return 3*t*t*p4
  + 3*t*(2-3*t)*p3
  + 3*tp*(1-3*t)*p2
  - 3*tp*tp*p1;
}

template <typename T>
bool getParameterCubicBezier_IntersectionWithPlane(
    double& t,
    const T& org,
    const T& nrm,
    const T& p1,
    const T& p2,
    const T& p3,
    const T& p4)
{
  double h1 = (p1-org).dot(nrm);
  double h2 = (p2-org).dot(nrm);
  double h3 = (p3-org).dot(nrm);
  double h4 = (p4-org).dot(nrm);
  double ref = fabs(h1) + fabs(h2) + fabs(h3) + fabs(h4);
  double eps = 1.0e-5;
  if( fabs(h1)<ref*eps && fabs(h4)<ref*eps ) return false;
  if( h1*h4>0 ) return false;
  if( fabs(h1)<ref*eps ){ t = 0.0; return true; }
  if( fabs(h4)<ref*eps ){ t = 1.0; return true; }
  t = 0.5;
  for(int itr=0;itr<10;++itr){
    double tp = 1.0-t;
    double h = t*t*t*h4 + 3*t*t*tp*h3 + 3*t*tp*tp*h2 + tp*tp*tp*h1;
    if( fabs(h)<1.0e-6*ref  ) return true;
    double dh = 3*t*t*h4 + 3*t*(2-3*t)*h3 + 3*tp*(1-3*t)*h2 - 3*tp*tp*h1;
    t -= (h/dh);
  }
  return false;
}

// Bezier
// p00: u=0 v=0
// p03: u=0 v=1
// p30: u=1 v=0
// p33: u=1 v=1
template <typename T>
T getPointSurfaceBezierCubic(
    double u, double v,
    const T& p00, const T& p01, const T& p02, const T& p03,
    const T& p10, const T& p11, const T& p12, const T& p13,
    const T& p20, const T& p21, const T& p22, const T& p23,
    const T& p30, const T& p31, const T& p32, const T& p33)
{
  double up = 1.0-u;
  double u3 = u*u*u;
  double u2 = 3*u*u*up;
  double u1 = 3*u*up*up;
  double u0 = up*up*up;
  ////
  double vp = 1.0-v;
  double v3 = v*v*v;
  double v2 = 3*v*v*vp;
  double v1 = 3*v*vp*vp;
  double v0 = vp*vp*vp;
  ////
  return
  + (u0*v0)*p00 + (u0*v1)*p01 + (u0*v2)*p02 + (u0*v3)*p03
  + (u1*v0)*p10 + (u1*v1)*p11 + (u1*v2)*p12 + (u1*v3)*p13
  + (u2*v0)*p20 + (u2*v1)*p21 + (u2*v2)*p22 + (u2*v3)*p23
  + (u3*v0)*p30 + (u3*v1)*p31 + (u3*v2)*p32 + (u3*v3)*p33;
}

template <typename T>
void getCubicBezierCurve(
    const int n,
    std::vector<T>& aP,
    const std::vector<T>& aCP)
{
  int ns = (int)(aCP.size()/3);
  aP.resize(ns*n+1);
  for (int is = 0; is<ns; is++){
    for (int i = 0; i<n; i++){
      double t = (double)i/n;
      aP[is*n+i] = getPointCubicBezierCurve(t, aCP[is*3+0], aCP[is*3+1], aCP[is*3+2], aCP[is*3+3]);
    }
  }
  aP[ns*n] = aCP[ns*3];
}


template <typename T>
T DeBoorBSpline(
    double u,
    int ndegree,
    const std::vector<T>& aCP,
    const std::vector<double>& aKnot)
{
  assert( ndegree>0 );
  assert( aKnot.size() == aCP.size()+ndegree+1 );
  const double eps = 1.0e-10;
  ////
  unsigned int iks;
  {
    for(iks=ndegree;iks<aKnot.size()-ndegree;++iks){
      double u0 = aKnot[iks];
      double u1 = aKnot[iks+1];
      if( u >= u0-eps && u <= u1+eps ){ break; }
    }
  }
  std::vector<T> aTmp;
  {
    aTmp.reserve(ndegree+1);
    for(unsigned int ik=iks-ndegree;ik<iks+1;++ik){
      assert(ik>=0);
      aTmp.push_back(aCP[ik]);
    }
  }
  for(int r=0;r<ndegree;++r){
    for(int j=0;j<ndegree-r;++j){
      double u0 = aKnot[j+iks-ndegree+1+r];
      double u1 = aKnot[j+iks+1];
//      assert(u>=u0-eps && u<u1+eps);
      double a = (u-u0)/(u1-u0);
      aTmp[j] = (1-a)*aTmp[j]+a*aTmp[j+1];
    }
  }
  return aTmp[0];
}

template <typename T>
void SampleBSpline(
    std::vector<T>& polyline0,
    const int nsmpl,
    const int ndegree,
    const std::vector<double>& aKnotFlat,
    const std::vector<T>& aCtrlPoint)
{
  polyline0.clear();
  double u0 = aKnotFlat[0];
  double u1 = aKnotFlat[aKnotFlat.size()-1];
  for(int i=0;i<nsmpl+1;++i){
    double u = (double)i*(u1-u0)/nsmpl+u0;
    T p = DeBoorBSpline(u, ndegree,aCtrlPoint,aKnotFlat);
    polyline0.push_back(p);
  }
}
  
}

#endif /* paramgeo_vec23_hpp */
