/*
 * Copyright (c) 2019 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */

#include "delfem2/vec3.h"
#include "delfem2/paramgeo_v23.h"

CVector3 getPointCoonsQuad_CubicBezier
(double u, double v,
 CVector3 aP[12])
{
  CVector3 p01u = getPointCubicBezierCurve(u, aP[0], aP[1], aP[2], aP[3]);
  CVector3 p32u = getPointCubicBezierCurve(u, aP[9], aP[8], aP[7], aP[6]);
  CVector3 p = (1-v)*p01u + v*p32u;
  ///
  CVector3 q03v = getPointCubicBezierCurve(v, aP[0], aP[11], aP[10], aP[9]);
  CVector3 q12v = getPointCubicBezierCurve(v, aP[3], aP[ 4], aP[ 5], aP[6]);
  CVector3 q = (1-u)*q03v + u*q12v;
  
  CVector3 r = (1-u)*(1-v)*aP[0] + u*(1-v)*aP[3] + u*v*aP[6] + (1-u)*v*aP[9];
  
  return p+q-r;
}

CVector3 getPointCoonsTri_CubicBezierEdge
(double u, double v, double w,
 CVector3 aP[9])
{
  CVector3 peu = getPointCubicBezierCurve(w/(1-u), aP[3], aP[4], aP[5], aP[6]);
  CVector3 pev = getPointCubicBezierCurve(u/(1-v), aP[6], aP[7], aP[8], aP[0]);
  CVector3 pew = getPointCubicBezierCurve(v/(1-w), aP[0], aP[1], aP[2], aP[3]);
  CVector3 pu = (1-u)*peu + u*aP[0];
  CVector3 pv = (1-v)*pev + v*aP[3];
  CVector3 pw = (1-w)*pew + w*aP[6];
  CVector3 pl = u*aP[0] + v*aP[3] + w*aP[6];
  
  //  double tmp = 1.0/(1.0/u + 1.0/v + 1.0/w);
  //  return (1.0/u)*tmp*pu + (1.0/v)*tmp*pv + (1.0/w)*tmp*pw;
  //  return
  return pu+pv+pw-2*pl;
  
}

CVector3 getPointHermetianQuad
(double u, double v,
 CVector3 aP[12])
{
  double u0 = +2*u*u*u-3*u*u+1;
  double u1 = -2*u*u*u+3*u*u;
  double du0 = +1*u*u*u-2*u*u+u;
  double du1 = +1*u*u*u-1*u*u;
  /////
  double v0 = +2*v*v*v-3*v*v+1;
  double v1 = -2*v*v*v+3*v*v;
  double dv0 = +1*v*v*v-2*v*v+v;
  double dv1 = +1*v*v*v-1*v*v;
  ///
  CVector3 p = aP[0]*u0*v0 + aP[3]*u1*v0 + aP[6]*u1*v1 + aP[9]*u0*v1;
  CVector3 q = 3*(aP[ 1]-aP[0])*du0*v0 + 3*(aP[3]-aP[2])*du1*v0 + 3*(aP[6]-aP[7])*du1*v1 + 3*(aP[8]-aP[ 9])*du0*v1;
  CVector3 r = 3*(aP[11]-aP[0])*u0*dv0 + 3*(aP[4]-aP[3])*u1*dv0 + 3*(aP[6]-aP[5])*u1*dv1 + 3*(aP[9]-aP[10])*u0*dv1;
  return p+q+r;
}

bool getParameterCubicBezier_IntersectionWithPlane
(double& t,
 const CVector3& org, const CVector3& nrm,
 const CVector3& p1, const CVector3& p2, const CVector3& p3, const CVector3& p4)
{
  double h1 = (p1-org)*nrm;
  double h2 = (p2-org)*nrm;
  double h3 = (p3-org)*nrm;
  double h4 = (p4-org)*nrm;
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

CVector3 getPointCubicBezierCurve
(double t,
 const CVector3& p1, const CVector3& p2, const CVector3& p3, const CVector3& p4)
{
  double tp = 1.0-t;
  return t*t*t*p4
  + 3*t*t*tp*p3
  + 3*t*tp*tp*p2
  + tp*tp*tp*p1;
}

CVector3 getTangentCubicBezierCurve
(double t,
 const CVector3& p1, const CVector3& p2, const CVector3& p3, const CVector3& p4)
{
  double tp = 1.0-t;
  return 3*t*t*p4
  + 3*t*(2-3*t)*p3
  + 3*tp*(1-3*t)*p2
  - 3*tp*tp*p1;
}

CVector3 QuadBilinear
(int iq, double r0, double r1,
 std::vector<int>& aQuad,
 std::vector<CVector3>& aPoint)
{
  int i0 = aQuad[iq*4+0];
  int i1 = aQuad[iq*4+1];
  int i2 = aQuad[iq*4+2];
  int i3 = aQuad[iq*4+3];
  return (1-r0)*(1-r1)*aPoint[i0] + r0*(1-r1)*aPoint[i1] + r0*r1*aPoint[i2] + (1-r0)*r1*aPoint[i3];
}


void getCubicBezierCurve
(const int n,
 std::vector<CVector3>& aP,
 const std::vector<CVector3>& aCP)
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

// p00: u=0 v=0
// p03: u=0 v=1
// p30: u=1 v=0
// p33: u=1 v=1
CVector3 getPointSurfaceBezierCubic
(double u, double v,
 const CVector3& p00, const CVector3& p01, const CVector3& p02, const CVector3& p03,
 const CVector3& p10, const CVector3& p11, const CVector3& p12, const CVector3& p13,
 const CVector3& p20, const CVector3& p21, const CVector3& p22, const CVector3& p23,
 const CVector3& p30, const CVector3& p31, const CVector3& p32, const CVector3& p33)
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

void getCubicBezierSurface
(const int n, // number of segment
 std::vector<CVector3>& aP,
 const std::vector<CVector3>& aCP)
{
  aP.resize((n+1)*(n+1));
  for (int i = 0; i<(n+1); ++i){
    for (int j = 0; j<(n+1); ++j){
      double u = (double)i/n;
      double v = (double)j/n;
      aP[i*(n+1)+j] = getPointSurfaceBezierCubic(u,v,
                                                 aCP[ 0], aCP[ 1], aCP[ 2], aCP[ 3],
                                                 aCP[ 4], aCP[ 5], aCP[ 6], aCP[ 7],
                                                 aCP[ 8], aCP[ 9], aCP[10], aCP[11],
                                                 aCP[12], aCP[13], aCP[14], aCP[15]);
    }
  }
}


void FlatKnot
(std::vector<double>& aKnotFlat,
 const std::vector<int>& aKnotMulti,
 const std::vector<double>& aKnot)
{
  assert(aKnot.size()==aKnotMulti.size());
  aKnotFlat.clear();
  for(int ik=0;ik<aKnot.size();++ik){
    for(int iik=0;iik<aKnotMulti[ik];++iik){
      aKnotFlat.push_back(aKnot[ik]);
    }
  }
}
