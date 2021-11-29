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
#include "delfem2/dfm2_inline.h"

#define NEARLY_ZERO 1.e-16

namespace delfem2 {


DFM2_INLINE void Nearest_Edge3_Point3(
    double nearest_position[3],
    const double point_pos[3], // point
    const double edge_pos0[3], // source
    const double edge_pos1[3]); // end

template<typename T>
CVec3<T> Nearest_Edge3_Point3(
  const CVec3<T> &p, // point
  const CVec3<T> &s, // source
  const CVec3<T> &e); // end

template<typename T>
CVec3<T> Nearest_Edge3_Point3(
  T &t,
  const CVec3<T> &p, // point
  const CVec3<T> &s,
  const CVec3<T> &e); // direction

template<class VEC>
void Nearest_Edge3_Line3(
  VEC &nearest_lineseg,
  VEC &nearest_line,
  const VEC &lineseg_start,
  const VEC &lineseg_end,
  const VEC &line_origin,
  const VEC &line_direction);

template<typename T>
CVec3<T> Nearest_Origin3_Edge3(
  const CVec3<T> &s, // start
  const CVec3<T> &e); // end

// r0==0 -> p0==org
// r0==1 -> p1==org
template<typename T>
CVec3<T> Nearest_Origin3_Edge3(
  double &r0,
  const CVec3<T> &p0, // start
  const CVec3<T> &p1); // end

template<typename T>
double Distance_Edge3_Edge3(
    const CVec3<T> &p0, const CVec3<T> &p1,
    const CVec3<T> &q0, const CVec3<T> &q1,
    double &ratio_p, double &ratio_q);

template<typename T>
bool IsContact_Edge3_Edge3_Proximity(
  int ino0, int ino1, int jno0, int jno1,
  const CVec3<T> &p0, const CVec3<T> &p1, const CVec3<T> &q0, const CVec3<T> &q1,
  const double delta);

} // end namespace delfem2

#ifndef DFM2_STATIC_LIBRARY
#  include "delfem2/geo_edge.cpp"
#endif

#endif // DFM2_GEO_EDGE_H
