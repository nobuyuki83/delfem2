#ifndef PARAMGEO_V23_H
#define PARAMGEO_V23_H

#include <stdio.h>

#include "delfem2/vec3.h"

CVector3 getPointCoonsQuad_CubicBezierEdge
(double u, double v,
 CVector3 aP[12]);

CVector3 getPointCoonsTri_CubicBezierEdge
(double u, double v, double w,
 CVector3 aP[9]);

CVector3 getPointHermetianQuad
(double u, double v,
 CVector3 aP[12]);

CVector3 getPointCubicBezierCurve
(double t,
 const CVector3& p1, const CVector3& p2, const CVector3& p3, const CVector3& p4);

CVector3 getTangentCubicBezierCurve
(double t,
 const CVector3& p1, const CVector3& p2, const CVector3& p3, const CVector3& p4);

bool getParameterCubicBezier_IntersectionWithPlane
(double& t,
 const CVector3& org, const CVector3& nrm,
 const CVector3& p1, const CVector3& p2, const CVector3& p3, const CVector3& p4);

// Bezier
CVector3 getPointSurfaceBezierCubic
(double u, double v,
 const CVector3& p00, const CVector3& p01, const CVector3& p02, const CVector3& p03,
 const CVector3& p10, const CVector3& p11, const CVector3& p12, const CVector3& p13,
 const CVector3& p20, const CVector3& p21, const CVector3& p22, const CVector3& p23,
 const CVector3& p30, const CVector3& p31, const CVector3& p32, const CVector3& p33);

void getCubicBezierCurve
(const int n,
 std::vector<CVector3>& aP,
 const std::vector<CVector3>& aCP);


void FlatKnot
(std::vector<double>& aKnotFlat,
 const std::vector<int>& aKnotMulti,
 const std::vector<double>& aKnot);

template <typename T>
T DeBoorBSpline
(double u,
 int ndegree,
 const std::vector<T>& aCP,
 const std::vector<double>& aKnot)
{
  assert( ndegree>0 );
  assert( aKnot.size() == aCP.size()+ndegree+1 );
  const double eps = 1.0e-10;
  ////
  int iks;
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
    for(int ik=iks-ndegree;ik<iks+1;++ik){
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
void SampleBSpline
(std::vector<T>& polyline0,
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

#endif /* paramgeo_vec23_hpp */
