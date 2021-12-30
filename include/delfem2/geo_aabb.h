/*
* Copyright (c) 2019 Nobuyuki Umetani
*
* This source code is licensed under the MIT license found in the
* LICENSE file in the root directory of this source tree.
*/

/**
 * @detail The order of dependency in delfem2:
 * aabb ->
 * line -> ray -> edge -> polyline ->
 * curve_quadratic -> curve_cubic -> curve_ndegree ->
 * plane < tri < quad
 */

#ifndef DFM2_AABB_H
#define DFM2_AABB_H

#include <cassert>
#include <cmath>
#include <iostream>
#include <vector>

#include "delfem2/dfm2_inline.h"

namespace delfem2 {

/**
 * check if the padded AABB of point3 contains the origin
 * @tparam VEC delfem2::CVec2, delfem2::CVec3, Eigen::Vector
 * @param p0
 * @param p1
 * @param p2
 * @param d padding size. How much AABB is padded
 * @return
 */
template<typename VEC>
bool IsContact_Orgin_AabbOfPoint3(
  const VEC &p0,
  const VEC &p1,
  const VEC &p2,
  typename VEC::Scalar d) {
  constexpr int ndim = VEC::SizeAtCompileTime;
  for (unsigned int idim = 0; idim < ndim; ++idim) {
    if (p0[idim] > +d && p1[idim] > +d && p2[idim] > +d) { return false; }
    if (p0[idim] < -d && p1[idim] < -d && p2[idim] < -d) { return false; }
  }
  return true;
}

/**
 * check if the padded AABB of point3 contains the origin
 * @tparam VEC delfem2::CVec2, delfem2::CVec3, Eigen::Vector
 * @param p0
 * @param p1
 * @param p2
 * @param d padding size. How much AABB is padded
 * @return
 */
template<typename VEC>
bool IsContact_Orgin_AabbOfPoint4(
  const VEC &p0,
  const VEC &p1,
  const VEC &p2,
  const VEC &p3,
  typename VEC::Scalar d) {
  constexpr int ndim = VEC::SizeAtCompileTime;
  for (unsigned int idim = 0; idim < ndim; ++idim) {
    if (p0[idim] > +d && p1[idim] > +d && p2[idim] > +d && p3[idim] > +d) { return false; }
    if (p0[idim] < -d && p1[idim] < -d && p2[idim] < -d && p3[idim] < -d) { return false; }
  }
  return true;
}

/**
 * check if the padded AABB of point3 contains the origin
 * @tparam VEC delfem2::CVec2, delfem2::CVec3, Eigen::Vector
 * @param p0
 * @param p1
 * @param p2
 * @param d padding size. How much AABB is padded
 * @return
 */
template<typename VEC>
bool IsContact_Orgin_AabbOfPoint5(
  const VEC &p0,
  const VEC &p1,
  const VEC &p2,
  const VEC &p3,
  const VEC &p4,
  typename VEC::Scalar d) {
  constexpr int ndim = VEC::SizeAtCompileTime;
  for (unsigned int idim = 0; idim < ndim; ++idim) {
    if (p0[idim] > +d && p1[idim] > +d && p2[idim] > +d && p3[idim] > +d && p4[idim] > +d) { return false; }
    if (p0[idim] < -d && p1[idim] < -d && p2[idim] < -d && p3[idim] < -d && p4[idim] < -d) { return false; }
  }
  return true;
}

/**
 * check if the padded AABB of point3 contains the origin
 * @tparam VEC delfem2::CVec2, delfem2::CVec3, Eigen::Vector
 * @param p0
 * @param p1
 * @param p2
 * @param d padding size. How much AABB is padded
 * @return
 */
template<typename VEC>
bool IsContact_Orgin_AabbOfPoint6(
  const VEC &p0,
  const VEC &p1,
  const VEC &p2,
  const VEC &p3,
  const VEC &p4,
  const VEC &p5,
  typename VEC::Scalar d) {
  constexpr int ndim = VEC::SizeAtCompileTime;
  for (unsigned int idim = 0; idim < ndim; ++idim) {
    if (p0[idim] > +d && p1[idim] > +d && p2[idim] > +d && p3[idim] > +d && p4[idim] > +d
      && p5[idim] > +d) { return false; }
    if (p0[idim] < -d && p1[idim] < -d && p2[idim] < -d && p3[idim] < -d && p4[idim] < -d
      && p5[idim] < -d) { return false; }
  }
  return true;
}

} // end namespace delfem2

#endif // DFM2_AABB_H
