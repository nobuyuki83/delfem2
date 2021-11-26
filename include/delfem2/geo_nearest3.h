/*
* Copyright (c) 2019 Nobuyuki Umetani
*
* This source code is licensed under the MIT license found in the
* LICENSE file in the root directory of this source tree.
*/


#ifndef DFM2_GEO_NEAREST3_H
#define DFM2_GEO_NEAREST3_H

#include <cassert>
#include <cmath>
#include <iostream>
#include <vector>

#include "delfem2/vec3.h"
#include "delfem2/dfm2_inline.h"

#define NEARLY_ZERO 1.e-16

namespace delfem2 {

template<class VEC, typename T>
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
// -----------------
// below: lineseg

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

DFM2_INLINE double Nearest_LineSeg_LineSeg_CCD_Iteration(
  double p[3],
  const CVec3d &p0s,
  const CVec3d &p0e,
  const CVec3d &p1s,
  const CVec3d &p1e,
  const CVec3d &q0s,
  const CVec3d &q0e,
  const CVec3d &q1s,
  const CVec3d &q1e,
  unsigned int nitr);

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

// above: lineseg
// --------------------------------------------------------------
// below: plane

template<typename T>
CVec3<T> Nearest_Plane3_Point3(
  const CVec3<T> &p, // point
  const CVec3<T> &o, // origin
  const CVec3<T> &n); // normal

template<typename T>
CVec3<T> Nearest_Orgin3_PlaneTri3(
  T &r0,
  T &r1,
  const CVec3<T> &q0,
  const CVec3<T> &q1,
  const CVec3<T> &q2);

// ---------------------
// below: triangle

DFM2_INLINE void Nearest_Triangle3_Point3(
  double nearest_position[3],
  double &r0,
  double &r1,
  const double ps[3], // origin point
  const double q0[3],
  const double q1[3],
  const double q2[3]);


template<typename T>
CVec3<T> Nearest_Origin3_Tri3(
    T &r0,
    T &r1,
    const CVec3<T> &q0,
    const CVec3<T> &q1,
    const CVec3<T> &q2);

// --------------------------------
// below: quad

template<typename T>
CVec3<T> nearst_Origin_Quad(
    double &s0, double &s1,
    const CVec3<T> &p,
    const CVec3<T> &q0,
    const CVec3<T> &q1,
    const CVec3<T> &q2,
    const CVec3<T> &q3);

template<typename T>
CVec3<T> Nearst_Origin3_Quad3(
  double &s0,
  double &s1,
  const CVec3<T> &q0,
  const CVec3<T> &q1,
  const CVec3<T> &q2,
  const CVec3<T> &q3);

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
#  include "delfem2/geo_nearest3.cpp"
#endif

#endif // DFM2_GEO_NEAREST3_H
