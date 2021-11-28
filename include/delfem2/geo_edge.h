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


DFM2_INLINE void Nearest_LineSeg3_Point3(
    double nearest_position[3],
    const double point_position[3], // point
    const double lineseg_s[3], // source
    const double lineseg_e[3]); // end

template<typename T>
CVec3<T> Nearest_LineSeg3_Point3(
  const CVec3<T> &p, // point
  const CVec3<T> &s, // source
  const CVec3<T> &e); // end

template<typename T>
CVec3<T> Nearest_LineSeg3_Point3(
  T &t,
  const CVec3<T> &p, // point
  const CVec3<T> &s,
  const CVec3<T> &e); // direction

template<class VEC, typename T>
void Nearest_LineSeg3_Line3(
  VEC &nearest_lineseg,
  VEC &nearest_line,
  const VEC &lineseg_start,
  const VEC &lineseg_end,
  const VEC &line_origin,
  const VEC &line_direction);

template<typename T>
CVec3<T> Nearest_Origin3_LineSeg3(
  const CVec3<T> &s, // start
  const CVec3<T> &e); // end

// r0==0 -> p0==org
// r0==1 -> p1==org
template<typename T>
CVec3<T> Nearest_Origin3_LineSeg3(
  double &r0,
  const CVec3<T> &p0, // start
  const CVec3<T> &p1); // end

template<typename T>
double DistanceEdgeEdge(
    const CVec3<T> &p0, const CVec3<T> &p1,
    const CVec3<T> &q0, const CVec3<T> &q1,
    double &ratio_p, double &ratio_q);

template<typename T>
bool IsContact_EE_Proximity(
    int ino0, int ino1, int jno0, int jno1,
    const CVec3<T> &p0, const CVec3<T> &p1, const CVec3<T> &q0, const CVec3<T> &q1,
    double delta);

} // end namespace delfem2

#ifndef DFM2_STATIC_LIBRARY
#  include "delfem2/geo_edge.cpp"
#endif

#endif // DFM2_GEO_NEAREST3_H
