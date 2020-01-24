/*
 * Copyright (c) 2019 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */

#ifndef DFM2_PARAMGEO_V23_H
#define DFM2_PARAMGEO_V23_H

#include <stdio.h>
#include "delfem2/vec3.h"

namespace delfem2 {
  
CVec3d QuadBilinear(int iq, double r0, double r1,
                      std::vector<int>& aQuad,
                      std::vector<CVec3d>& aPoint);

  
CVec3d getPointCoonsQuad_CubicBezier(double u, double v,
                                        CVec3d aP[12]);
  
void getCubicBezierSurface(const int n, // number of segment
                           std::vector<CVec3d>& aP,
                           const std::vector<CVec3d>& aCP);

CVec3d getPointCoonsQuad_CubicBezierEdge(double u, double v,
                                           CVec3d aP[12]);

CVec3d getPointCoonsTri_CubicBezierEdge(double u, double v, double w,
                                          CVec3d aP[9]);

CVec3d getPointHermetianQuad(double u, double v,
                               CVec3d aP[12]);

CVec3d getPointCubicBezierCurve
(double t,
 const CVec3d& p1, const CVec3d& p2, const CVec3d& p3, const CVec3d& p4);

CVec3d getTangentCubicBezierCurve
(double t,
 const CVec3d& p1, const CVec3d& p2, const CVec3d& p3, const CVec3d& p4);

bool getParameterCubicBezier_IntersectionWithPlane
(double& t,
 const CVec3d& org, const CVec3d& nrm,
 const CVec3d& p1, const CVec3d& p2, const CVec3d& p3, const CVec3d& p4);

// Bezier
CVec3d getPointSurfaceBezierCubic
(double u, double v,
 const CVec3d& p00, const CVec3d& p01, const CVec3d& p02, const CVec3d& p03,
 const CVec3d& p10, const CVec3d& p11, const CVec3d& p12, const CVec3d& p13,
 const CVec3d& p20, const CVec3d& p21, const CVec3d& p22, const CVec3d& p23,
 const CVec3d& p30, const CVec3d& p31, const CVec3d& p32, const CVec3d& p33);

void getCubicBezierCurve(const int n,
                         std::vector<CVec3d>& aP,
                         const std::vector<CVec3d>& aCP);

void FlatKnot(std::vector<double>& aKnotFlat,
              const std::vector<int>& aKnotMulti,
              const std::vector<double>& aKnot);

template <typename T>
T DeBoorBSpline(double u,
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
void SampleBSpline(std::vector<T>& polyline0,
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
