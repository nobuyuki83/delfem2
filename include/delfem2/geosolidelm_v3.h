/*
* Copyright (c) 2019 Nobuyuki Umetani
*
* This source code is licensed under the MIT license found in the
* LICENSE file in the root directory of this source tree.
*/


#ifndef DFM2_GEOSOLIDELM_V3_H
#define DFM2_GEOSOLIDELM_V3_H

#include <cassert>
#include <cmath>
#include <iostream>
#include <vector>

#include "delfem2/vec3.h"
#include "delfem2/dfm2_inline.h"

#define NEARLY_ZERO 1.e-16

namespace delfem2 {

template<typename T>
double Height(
    const CVec3<T> &v1,
    const CVec3<T> &v2,
    const CVec3<T> &v3,
    const CVec3<T> &v4);

// ---------------------------------------------------------

template<typename T>
bool barycentricCoord_Origin_Tet(
    double &r0,
    double &r1,
    double &r2,
    const CVec3<T> &p0,
    const CVec3<T> &p1,
    const CVec3<T> &p2,
    const CVec3<T> &p3);

template<typename T>
bool barycentricCoord_Origin_Pyramid(
    double &r0,
    double &r1,
    double &r2,
    const CVec3<T> &p0,
    const CVec3<T> &p1,
    const CVec3<T> &p2,
    const CVec3<T> &p3,
    const CVec3<T> &p4);

template<typename T>
bool barycentricCoord_Origin_Wedge(
    double &r0,
    double &r1,
    double &r2,
    const CVec3<T> &p0,
    const CVec3<T> &p1,
    const CVec3<T> &p2,
    const CVec3<T> &p3,
    const CVec3<T> &p4,
    const CVec3<T> &p5);

template<typename T>
CVec3<T> positionBarycentricCoord_Pyramid(
    double r0,
    double r1,
    double r2,
    const CVec3<T> &p0,
    const CVec3<T> &p1,
    const CVec3<T> &p2,
    const CVec3<T> &p3,
    const CVec3<T> &p4);

template<typename T>
CVec3<T> positionBarycentricCoord_Wedge(
    double r0,
    double r1,
    double r2,
    const CVec3<T> &p0,
    const CVec3<T> &p1,
    const CVec3<T> &p2,
    const CVec3<T> &p3,
    const CVec3<T> &p4,
    const CVec3<T> &p5);

template<typename T>
void iteration_barycentricCoord_Origin_Solid(
    double &r0,
    double &r1,
    double &r2,
    const CVec3<T> &q, // q=positionBarycentricCoord_Wedge
    const CVec3<T> &dpdr0,
    const CVec3<T> &dpdr1,
    const CVec3<T> &dpdr2,
    double damp = 1.0);


// -----------------------------------------------

template<typename T>
T Volume_OrgTet(
    const CVec3<T> &v1,
    const CVec3<T> &v2,
    const CVec3<T> &v3);

/**
 * Volume of a tetrahedra
 * @tparam VEC dfm2::CVec3, Eigen::Vector3, std::array<*,3>
 */
template<typename VEC, typename T = typename VEC::Scalar>
T Volume_Tet(
    const VEC &v0,
    const VEC &v1,
    const VEC &v2,
    const VEC &v3) {
  T v = (v1[0] - v0[0]) * ((v2[1] - v0[1]) * (v3[2] - v0[2]) - (v3[1] - v0[1]) * (v2[2] - v0[2]))
      + (v1[1] - v0[1]) * ((v2[2] - v0[2]) * (v3[0] - v0[0]) - (v3[2] - v0[2]) * (v2[0] - v0[0]))
      + (v1[2] - v0[2]) * ((v2[0] - v0[0]) * (v3[1] - v0[1]) - (v3[0] - v0[0]) * (v2[1] - v0[1]));
  return v * 0.16666666666666666666666666666667;
}

template<typename VEC>
std::array<VEC, 4> DiffDeformationGradient(
    const VEC &P0,
    const VEC &P1,
    const VEC &P2,
    const VEC &P3);

template<typename T>
double Volume_Pyramid(
    const CVec3<T> &p0,
    const CVec3<T> &p1,
    const CVec3<T> &p2,
    const CVec3<T> &p3,
    const CVec3<T> &p4);

template<typename T>
double Volume_Wedge(
    const CVec3<T> &p0,
    const CVec3<T> &p1,
    const CVec3<T> &p2,
    const CVec3<T> &p3,
    const CVec3<T> &p4,
    const CVec3<T> &p5);

// ---------------------------------------------

template<typename T>
double SqareLongestEdgeLength(
    const CVec3<T> &ipo0,
    const CVec3<T> &ipo1,
    const CVec3<T> &ipo2,
    const CVec3<T> &ipo3);

template<typename T>
double LongestEdgeLength(
    const CVec3<T> &ipo0,
    const CVec3<T> &ipo1,
    const CVec3<T> &ipo2,
    const CVec3<T> &ipo3);

template<typename T>
double SqareShortestEdgeLength(
    const CVec3<T> &ipo0,
    const CVec3<T> &ipo1,
    const CVec3<T> &ipo2,
    const CVec3<T> &ipo3);

template<typename T>
double ShortestEdgeLength(
    const CVec3<T> &ipo0,
    const CVec3<T> &ipo1,
    const CVec3<T> &ipo2,
    const CVec3<T> &ipo3);

// ----------------------------------------------------------
// here starts std::vector<CVector3>

template<typename T>
double Volume_Tet(
    int iv1,
    int iv2,
    int iv3,
    int iv4,
    const std::vector<CVec3<T> > &aPoint);

template<typename T>
double SolidAngleTri(
    const CVec3<T> &v1,
    const CVec3<T> &v2,
    const CVec3<T> &v3);

} // end namespace delfem2

#ifndef DFM2_STATIC_LIBRARY
#  include "delfem2/geosolidelm_v3.cpp"
#endif

// ----------------------------
// local function interface

namespace delfem2::solidelm {
  template<typename REAL>
  DFM2_INLINE void MyInverse_Mat3(
      REAL Ainv[9],
      const REAL A[9]);
}

// ---------------------------
// implementation of exposed functions

template <typename VEC>
std::array<VEC, 4> delfem2::DiffDeformationGradient(
    const VEC &P0,
    const VEC &P1,
    const VEC &P2,
    const VEC &P3) {
  using SCALAR = typename VEC::Scalar;
  const SCALAR Basis0[9] = {
      P1[0] - P0[0], P2[0] - P0[0], P3[0] - P0[0],
      P1[1] - P0[1], P2[1] - P0[1], P3[1] - P0[1],
      P1[2] - P0[2], P2[2] - P0[2], P3[2] - P0[2] };
  SCALAR Bi0[9];
  ::delfem2::solidelm::MyInverse_Mat3(Bi0, Basis0);
  return {
      VEC{-Bi0[0] - Bi0[3] - Bi0[6],
          -Bi0[1] - Bi0[4] - Bi0[7],
          -Bi0[2] - Bi0[5] - Bi0[8]},
      VEC{Bi0[0], Bi0[1], Bi0[2]},
      VEC{Bi0[3], Bi0[4], Bi0[5]},
      VEC{Bi0[6], Bi0[7], Bi0[8]}};
}

#endif // DFM2_GEOSOLIDELM_V3_H
