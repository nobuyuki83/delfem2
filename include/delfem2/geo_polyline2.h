/*
 * Copyright (c) 2019 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */

#ifndef DFM2_GEO_POLYLINE2_H
#define DFM2_GEO_POLYLINE2_H

#include <cstdlib>
#include <cassert>
#include <iostream>
#include <cmath>
#include <vector>
#include <climits>

#include "delfem2/vec2.h"
#include "delfem2/pgeo.h"
#include "delfem2/dfm2_inline.h"

// -----------------------------------------------------

namespace delfem2 {

//! @brief translate all the points
template<typename T>
void Translate(
    std::vector<CVec2<T> > &aP,
    double dx,
    double dy);

template<typename T>
DFM2_INLINE void Rotate(
    std::vector<CVec2<T> > &aP,
    double dt);

template<typename T>
void Polyline_CubicBezierCurve(
    std::vector<CVec2<T> > &aP,
    int n,
    const std::vector<CVec2<T> > &aCP);

template<typename T>
void Polyline_BezierCubic(
    std::vector<CVec2<T> > &aP,
    unsigned int n,
    const CVec2<T> &p1,
    const CVec2<T> &p2,
    const CVec2<T> &p3,
    const CVec2<T> &p4);

template<typename T>
void Polyline_BezierQuadratic(
    std::vector<CVec2<T>> &aP,
    unsigned int n,
    const CVec2<T> &p1,
    const CVec2<T> &p2,
    const CVec2<T> &p3);

template<typename T>
std::vector<CVec2<T> > Polyline_Resample_Polyline(
    const std::vector<CVec2<T> > &stroke0,
    double l);

/**
 *
 * @tparam VEC delfem2::CVecX or Eigen::VectorX
 * @param polyline
 * @param scr
 * @return
 */
template<typename VEC>
unsigned int FindNearestPointInPoints(
    const std::vector<VEC>& polyline,
    const VEC& scr){
  unsigned int idx_point_min_dist = UINT_MAX;
  float dist_min = -1;
  for(unsigned int ip=0;ip<polyline.size();++ip){
    auto dist = (scr-polyline[ip]).norm();
    if( dist_min < 0 || dist < dist_min ){
      idx_point_min_dist = ip;
      dist_min = dist;
    }
  }
  return idx_point_min_dist;
}

/**
 *
 * @tparam VEC delfem2::CVecX or Eigen::VectorX
 * @param polyline
 * @param scr
 * @return
 */
template<typename VEC>
VEC FindNearestPointInPolyline(
    const std::vector<VEC>& polyline,
    const VEC& scr){
  float dist_min = -1;
  VEC p_min;
  for(unsigned int ip=0;ip<polyline.size()-1;++ip){
    unsigned int jp = ip+1;
    VEC p_near = GetNearest_LineSeg_Point(scr, polyline[ip], polyline[jp]);
    float dist = (p_near-scr).norm();
    if( dist_min < 0 || dist < dist_min ){
      dist_min = dist;
      p_min = p_near;
    }
  }
  return p_min;
}

} // namespace delfem2

#ifndef DFM2_STATIC_LIBRARY
#  include "delfem2/geo_polyline2.cpp"
#endif

#endif // DFM2_GEO_POLYLINE2_H


