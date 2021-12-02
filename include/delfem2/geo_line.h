/*
* Copyright (c) 2019 Nobuyuki Umetani
*
* This source code is licensed under the MIT license found in the
* LICENSE file in the root directory of this source tree.
*/


#ifndef DFM2_GEO_LINE_H
#define DFM2_GEO_LINE_H

#include <cassert>
#include <cmath>
#include <iostream>
#include <vector>

#include "delfem2/vec3.h"
#include "delfem2/dfm2_inline.h"

#define NEARLY_ZERO 1.e-16

namespace delfem2 {


// get parameter 't' of the line against point. t=0 is po_s, t=1 is po_e
// this one has implementation in header because GetDist_LineSeg_Point below refers this
template <class VEC>
double FindNearestPointParameter_Line_Point(
    const VEC &po_c,
    const VEC &po_s,
    const VEC &po_e) {
  using SCALAR = typename VEC::Scalar;
  const VEC &es = po_e - po_s;
  const SCALAR a = es.squaredNorm();
  const SCALAR b = es.dot(po_s - po_c);
  return -b / a;
}


template<class VEC>
VEC Nearest_Line3_Point3(
  const VEC &point,
  const VEC &line_src,
  const VEC &line_dir);

template<typename T>
CVec3<T> Nearest_Line3_Point3(
  double &t,
  const CVec3<T> &p, // point
  const CVec3<T> &s, // source
  const CVec3<T> &d); // direction


/**
 *  @param D (out) scaling factor
 *  @param Da (out) nearest point scaled by D on line A
 *  @param Db (out) nearest point scaled by D on line B
 *  @param Dta (out) parameter for nearest pont one line A. Da = D*pa_ + Dta*va
 *  @param Dtb (out) parameter for nearest pont one line B. Db = D*pb_ + Dtb*vb
 */
template<class VEC, typename T>
void Nearest_Line3_Line3(
  T &D,
  VEC &Da,
  VEC &Db,
  T &Dta,
  T &Dtb,
  const VEC &pa_,
  const VEC &va,
  const VEC &pb_,
  const VEC &vb);


/**
 *  @param[out] scale scaling factor
 *  @param[out] scaled_nearest_a nearest point scaled by D on line A
 *  @param[out] scaled_nearest_b nearest point scaled by D on line B
 */
template<class VEC, typename T>
void Nearest_Line3_Line3(
  T &scale,
  VEC &scaled_neraest_a,
  VEC &scaled_nearest_b,
  const VEC &line_org_a,
  const VEC &line_dir_a,
  const VEC &line_org_b,
  const VEC &line_dir_b);

// above: line
// --------------------------------------------------------------
// below: plane
// --------------------------------
// below: circle

/**
 * @param p0 (out)  nearest point on line
 * @param q0 (out)  nearest point on circle
 */
template<typename T>
void Nearest_Line_Circle(
    CVec3<T> &p0,
    CVec3<T> &q0,
    const CVec3<T> &src,
    const CVec3<T> &dir,
    const CVec3<T> &org,
    const CVec3<T> &normal,
    T rad);


} // end namespace delfem2

#ifndef DFM2_STATIC_LIBRARY
#  include "delfem2/geo_line.cpp"
#endif

#endif // DFM2_GEO_NEAREST3_H
