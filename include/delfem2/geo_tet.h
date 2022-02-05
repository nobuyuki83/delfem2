/*
* Copyright (c) 2019 Nobuyuki Umetani
*
* This source code is licensed under the MIT license found in the
* LICENSE file in the root directory of this source tree.
*/


#ifndef DFM2_GEO_TET_H
#define DFM2_GEO_TET_H

#include <cassert>
#include <cmath>
#include <iostream>
#include <vector>
#include <array>

#include "delfem2/dfm2_inline.h"
#include "delfem2/geo_meta_funcs.h"

namespace delfem2 {


//! Hight of a tetrahedra
template<typename VEC, typename T = value_type<VEC>>
double Height_Tet(
    const VEC &v1,
    const VEC &v2,
    const VEC &v3,
    const VEC &v4);

// ---------------------------------------------------------

template<typename VEC, typename T = value_type<VEC>>
bool barycentricCoord_Origin_Tet(
    T &r0,
    T &r1,
    T &r2,
    const VEC &p0,
    const VEC &p1,
    const VEC &p2,
    const VEC &p3);

// -----------------------------------------------

//! Volume of a tetrahedra v0 is orgin
template<typename VEC, typename T = value_type<VEC>>
T Volume_OrgTet(
    const VEC &v1,
    const VEC &v2,
    const VEC &v3);

// -----------------------------------------------

/**
 * Volume of a tetrahedra
 * @tparam VEC dfm2::CVec3, Eigen::Vector3, std::array<*,3>, * [3]
 * example: Volume_Tet<double [3], double >(v0,...)
 */
template<typename VEC, typename T = value_type<VEC>>
T Volume_Tet(
    const VEC &v0,
    const VEC &v1,
    const VEC &v2,
    const VEC &v3);

// -----------------------------------------------

/**
 * three basis of a tetrahedra (origin is the first vertex)
 * @tparam VEC dfm2::CVec3, Eigen::Vector3, std::array<*,3>, * [3]
 * example: Mat3_3BasesOfTet<double [3], double >(v0,...)
 */
template <typename VEC, typename REAL = value_type<VEC>>
std::array<REAL,9> Mat3_3BasesOfTet(
    const VEC &p0,
    const VEC &p1,
    const VEC &p2,
    const VEC &p3)
{
  return  {
      p1[0]-p0[0], p2[0]-p0[0], p3[0]-p0[0],
      p1[1]-p0[1], p2[1]-p0[1], p3[1]-p0[1],
      p1[2]-p0[2], p2[2]-p0[2], p3[2]-p0[2]};
}

// -----------------------------------------------

template<typename VEC, typename T = value_type<VEC>>
std::array<T, 9> DeformationGradientOfTet(
    const VEC &P0,
    const VEC &P1,
    const VEC &P2,
    const VEC &P3,
    const VEC &p0,
    const VEC &p1,
    const VEC &p2,
    const VEC &p3);

// -----------------------------------------------

template<typename VEC, typename T = value_type<VEC>>
void DiffDeformationGradientOfTet(
    T dF[4][3],
    const VEC &P0,
    const VEC &P1,
    const VEC &P2,
    const VEC &P3);

// ---------------------------------------------

template<typename VEC, typename T = value_type<VEC>>
T LongestSquaredEdgeLengthOfTet(
    const VEC &ipo0,
    const VEC &ipo1,
    const VEC &ipo2,
    const VEC &ipo3);

// -----------------------------------------------

template<typename VEC, typename T = value_type<VEC>>
T ShortestSquaredEdgeLengthOfTet(
    const VEC &ipo0,
    const VEC &ipo1,
    const VEC &ipo2,
    const VEC &ipo3);

// ----------------------------------------------------------
// here starts std::vector<CVector3>

template<typename VEC, typename T = value_type<VEC>>
T Volume_Tet(
    int iv1,
    int iv2,
    int iv3,
    int iv4,
    const std::vector<VEC> &aPoint) {
  return Volume_Tet(aPoint[iv1], aPoint[iv2], aPoint[iv3], aPoint[iv4]);
}

} // end namespace delfem2

#ifndef DFM2_STATIC_LIBRARY
#  include "delfem2/geo_tet.cpp"
#endif

#endif // DFM2_GEO_TET_H
