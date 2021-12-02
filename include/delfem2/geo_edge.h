/*
* Copyright (c) 2019 Nobuyuki Umetani
*
* This source code is licensed under the MIT license found in the
* LICENSE file in the root directory of this source tree.
*/


#ifndef DFM2_GEO_EDGE_H
#define DFM2_GEO_EDGE_H

#include <cassert>
#include <cmath>
#include <iostream>
#include <vector>

#include "delfem2/vec3.h"
#include "delfem2/vec2.h"
#include "delfem2/dfm2_inline.h"

#define NEARLY_ZERO 1.e-16

namespace delfem2 {

// r0==0 -> p0==org
// r0==1 -> p1==org
template<class VEC>
VEC Nearest_Origin_Edge(
  typename VEC::Scalar &r0,
  const VEC &p0, // start
  const VEC &p1); // end

template<class VEC>
VEC Nearest_Origin_Edge(
  const VEC &s, // start
  const VEC &e); // end

//! @details  get parameter 't' of the line against point. t=0 is po_s, t=1 is po_e
template<class VEC>
VEC Nearest_Edge_Point(
  const VEC &po_c,
  const VEC &po_s,
  const VEC &po_e);

template<class VEC>
VEC Nearest_Edge_Point(
  typename VEC::Scalar &t,
  const VEC &p, // point
  const VEC &s,
  const VEC &e); // direction

//! @details  get parameter 't' of the line against point. t=0 is po_s, t=1 is po_e
template<class VEC>
typename VEC::Scalar Distance_Edge_Point(
  const VEC &po_c,
  const VEC &po_s,
  const VEC &po_e);

// 2D and 3D
// --------------------------
// 2D

template<typename T>
bool IsCross_LineSeg_LineSeg(
  const CVec2<T> &po_s0,
  const CVec2<T> &po_e0,
  const CVec2<T> &po_s1,
  const CVec2<T> &po_e1);

template<typename T>
double Distance_Edge2_Edge2(
  const CVec2<T> &po_s0,
  const CVec2<T> &po_e0,
  const CVec2<T> &po_s1,
  const CVec2<T> &po_e1);

// above: 2D
// ---------------------------
// below: 3D

DFM2_INLINE void Nearest_Edge3_Point3(
  double nearest_position[3],
  const double point_pos[3], // point
  const double edge_pos0[3], // source
  const double edge_pos1[3]); // end

template<class VEC>
void Nearest_Edge3_Line3(
  VEC &nearest_lineseg,
  VEC &nearest_line,
  const VEC &lineseg_start,
  const VEC &lineseg_end,
  const VEC &line_origin,
  const VEC &line_direction);

template<typename T>
double Distance_Edge3_Edge3(
  const CVec3<T> &p0,
  const CVec3<T> &p1,
  const CVec3<T> &q0,
  const CVec3<T> &q1,
  double &ratio_p,
  double &ratio_q);

template<typename T>
bool IsContact_Edge3_Edge3_Proximity(
  int ino0, int ino1, int jno0, int jno1,
  const CVec3<T> &p0,
  const CVec3<T> &p1,
  const CVec3<T> &q0,
  const CVec3<T> &q1,
  const double delta);

} // end namespace delfem2

#ifndef DFM2_STATIC_LIBRARY
#  include "delfem2/geo_edge.cpp"
#endif

#endif // DFM2_GEO_EDGE_H
