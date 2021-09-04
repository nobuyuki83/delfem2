/*
 * Copyright (c) 2019 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */

#include <cstdlib>

#include "delfem2/geo_polyline2.h"

// ----------------------------------------------------------------------------------
// std::vector starts from here

template<typename T>
void delfem2::Polyline_CubicBezierCurve(
    std::vector< CVec2<T> >& aP,
    const int n,
    const std::vector< CVec2<T> >& aCP) {
  int ns = (int) (aCP.size() / 3);
  aP.resize(ns * n + 1);
  for (int is = 0; is < ns; is++) {
    for (int i = 0; i < n; i++) {
      double t = (double) i / n;
      aP[is * n + i] = PointOnCubicBezierCurve(t,
                                               aCP[is * 3 + 0],
                                               aCP[is * 3 + 1],
                                               aCP[is * 3 + 2],
                                               aCP[is * 3 + 3]);
    }
  }
  aP[ns * n] = aCP[ns * 3];
}
#ifdef DFM2_STATIC_LIBRARY
template void delfem2::Polyline_CubicBezierCurve(std::vector<CVec2d>& aP,
                                                 const int n,
                                                 const std::vector<CVec2d>& aCP);
#endif


// ------------------------------

template<typename T>
void delfem2::Polyline_BezierCubic(
    std::vector< CVec2<T> >& aP,
    const unsigned int n,
    const CVec2<T> &p1,
    const CVec2<T> &p2,
    const CVec2<T> &p3,
    const CVec2<T> &p4) {
  aP.resize(n);
  for (unsigned int i = 0;i < n; ++i) {
    const double t = (double) i / (n - 1.0);
    aP[i] = PointOnCubicBezierCurve(
        t,
        p1, p2, p3, p4);
  }
}
template void delfem2::Polyline_BezierCubic(
    std::vector<CVec2d> &aP,
    const unsigned int n,
    const CVec2d &p1,
    const CVec2d &p2,
    const CVec2d &p3,
    const CVec2d &p4);

// --------------

template<typename T>
void delfem2::Polyline_BezierQuadratic(
    std::vector<CVec2 < T>>& aP,
    const unsigned int n,
    const CVec2<T> &p1,
    const CVec2<T> &p2,
    const CVec2<T> &p3) {
  aP.resize(n);
  for (unsigned int i = 0; i < n; ++i) {
    const double t = (double) i / (n - 1.0);
    aP[i] = PointOnQuadraticBezierCurve(
        t,
        p1, p2, p3 );
  }
}
#ifdef DFM2_STATIC_LIBRARY
template void delfem2::Polyline_BezierQuadratic(
    std::vector<CVec2d>& aP,
    const unsigned int n,
    const CVec2d& p1,
    const CVec2d& p2,
    const CVec2d& p3);
#endif

// ---------------------------------------

template<typename T>
void delfem2::Translate(
    std::vector<CVec2<T> >& aP,
    double dx,
    double dy) {
  for ( auto &ip : aP) {
    ip.p[0] += dx;
    ip.p[1] += dy;
  }
}
#ifdef DFM2_STATIC_LIBRARY
template void delfem2::Translate(std::vector<CVec2d>& aP, double dx, double dy);
#endif

// ----------------------------

namespace delfem2 {

template<>
DFM2_INLINE void Rotate(std::vector<CVec2d> &aP, double dt) {
  for (auto &ip : aP) {
    double x0 = ip.p[0];
    double y0 = ip.p[1];
    ip.p[0] = cos(dt) * x0 - sin(dt) * y0;
    ip.p[1] = sin(dt) * x0 + cos(dt) * y0;
  }
}

}

// ----------------------------

template<typename T>
std::vector< delfem2::CVec2<T> >
delfem2::Polyline_Resample_Polyline(
    const std::vector< CVec2<T> >& stroke0,
    double l)
{
  if (stroke0.empty()) {
    std::vector<CVec2<T>> a;
    return a;
  }
  std::vector<CVec2 < T>> stroke;
  stroke.push_back( stroke0[0] );
  int jcur = 0;
  double rcur = 0;
  double lcur = l;
  for(;;){
    if( jcur >= (int)stroke0.size()-1 ) break;
    double lenj = (stroke0[jcur + 1] - stroke0[jcur]).norm();
    double lenjr = lenj * (1.0 - rcur);
    if( lenjr > lcur ){ // put point in this segment
      rcur += lcur / lenj;
      stroke.push_back((1-rcur)*stroke0[jcur] +rcur *stroke0[jcur + 1]);
      lcur = l;
    }
    else{ // next segment
      lcur -= lenjr;
      rcur = 0;
      jcur++;
    }
  }
//  stroke.push_back( stroke0.back() );
  return stroke;
}
#ifdef DFM2_STATIC_LIBRARY
template std::vector<delfem2::CVec2d> delfem2::Polyline_Resample_Polyline(
    const std::vector<CVec2d>& stroke0, double l);
#endif

// -------------------------------


