/*
 * Copyright (c) 2019 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */

#include "delfem2/geo_polyline.h"

#include "delfem2/vec2.h"  // instantiation for vec2

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

template<class VEC>
std::vector<VEC> delfem2::Polyline_Resample_Polyline(
    const std::vector<VEC>& stroke0,
    typename VEC::Scalar l)
{
  if (stroke0.empty()) {
    return std::vector<VEC>{};
  }
  std::vector<VEC> stroke;
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
  return stroke;
}
#ifdef DFM2_STATIC_LIBRARY
template std::vector<delfem2::CVec2d> delfem2::Polyline_Resample_Polyline(
    const std::vector<CVec2d>& stroke0, double l);
template std::vector<delfem2::CVec2f> delfem2::Polyline_Resample_Polyline(
    const std::vector<CVec2f>& stroke0, float l);
#endif

// -------------------------------


