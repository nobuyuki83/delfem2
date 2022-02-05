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

#include "delfem2/dfm2_inline.h"
#include "delfem2/geo_meta_funcs.h"

#define NEARLY_ZERO 1.e-16

namespace delfem2 {


//! Hight of a tetrahedra
template<typename VEC, typename T = value_type<VEC>>
double Height_Tet(
    const VEC &v1,
    const VEC &v2,
    const VEC &v3,
    const VEC &v4) {
  T n[3];
  Normal_Tri3(
      n,
      v1, v2, v3);
  Normalize3(n);
  return (v4[0] - v1[0]) * n[0] + (v4[1] - v1[1]) * n[1] + (v4[2] - v1[2]) * n[2];
}

// ---------------------------------------------------------

template<typename T>
bool barycentricCoord_Origin_Tet(
    double &r0,
    double &r1,
    double &r2,
    const CVec3<T> &p0,
    const CVec3<T> &p1,
    const CVec3<T> &p2,
    const CVec3<T> &p3) {
  double v0 = Volume_OrgTet(p1, p2, p3);
  double v1 = Volume_OrgTet(p2, p0, p3);
  double v2 = Volume_OrgTet(p1, p3, p0);
  double v3 = Volume_OrgTet(p1, p0, p2);
  double vt_inv = 1.0 / (v0 + v1 + v2 + v3);
  r0 = v0 * vt_inv;
  r1 = v1 * vt_inv;
  r2 = v2 * vt_inv;
  return true;
}

// -----------------------------------------------

//! Volume of a tetrahedra v0 is orgin
template<typename VEC, typename T = value_type<VEC>>
T Volume_OrgTet(
    const VEC &v1,
    const VEC &v2,
    const VEC &v3) {
  T v = v1[0] * (v2[1] * v3[2] - v3[1] * v2[2])
      + v1[1] * (v2[2] * v3[0] - v3[2] * v2[0])
      + v1[2] * (v2[0] * v3[1] - v3[0] * v2[1]);
  return v * static_cast<T>(1.0 / 6.0);
}

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
    const VEC &v3) {
  T v = (v1[0] - v0[0]) * ((v2[1] - v0[1]) * (v3[2] - v0[2]) - (v3[1] - v0[1]) * (v2[2] - v0[2]))
      + (v1[1] - v0[1]) * ((v2[2] - v0[2]) * (v3[0] - v0[0]) - (v3[2] - v0[2]) * (v2[0] - v0[0]))
      + (v1[2] - v0[2]) * ((v2[0] - v0[0]) * (v3[1] - v0[1]) - (v3[0] - v0[0]) * (v2[1] - v0[1]));
  return v * 0.16666666666666666666666666666667;
}

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
    const VEC &ipo3) {
  T edge1, edge2;
  edge1 = SquareDistance3(ipo0, ipo1);
  edge2 = SquareDistance3(ipo0, ipo2);
  if (edge2 > edge1) edge1 = edge2;
  edge2 = SquareDistance3(ipo0, ipo3);
  if (edge2 > edge1) edge1 = edge2;
  edge2 = SquareDistance3(ipo1, ipo2);
  if (edge2 > edge1) edge1 = edge2;
  edge2 = SquareDistance3(ipo1, ipo3);
  if (edge2 > edge1) edge1 = edge2;
  edge2 = SquareDistance3(ipo2, ipo3);
  if (edge2 > edge1) edge1 = edge2;
  return edge1;
}

template<typename VEC, typename T = value_type<VEC>>
double ShortestSquaredEdgeLengthOfTet(
    const VEC &ipo0,
    const VEC &ipo1,
    const VEC &ipo2,
    const VEC &ipo3) {
  double edge1, edge2;
  edge1 = SquareDistance3(ipo0, ipo1);
  edge2 = SquareDistance3(ipo0, ipo2);
  if (edge2 < edge1) edge1 = edge2;
  edge2 = SquareDistance3(ipo0, ipo3);
  if (edge2 < edge1) edge1 = edge2;
  edge2 = SquareDistance3(ipo1, ipo2);
  if (edge2 < edge1) edge1 = edge2;
  edge2 = SquareDistance3(ipo1, ipo3);
  if (edge2 < edge1) edge1 = edge2;
  edge2 = SquareDistance3(ipo2, ipo3);
  if (edge2 < edge1) edge1 = edge2;
  return edge1;
}

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

// =====================================
// below: implementation of exposed functions

namespace delfem2::geo_tet{

template<typename REAL>
void Inverse_Mat3(
    REAL A[9]) {
  const REAL B[9] = {
      A[0], A[1], A[2],
      A[3], A[4], A[5],
      A[6], A[7], A[8]};
  const REAL det =
      +B[0] * B[4] * B[8] + B[3] * B[7] * B[2] + B[6] * B[1] * B[5]
          - B[0] * B[7] * B[5] - B[6] * B[4] * B[2] - B[3] * B[1] * B[8];
  const REAL inv_det = 1 / det;
  A[0] = inv_det * (B[4] * B[8] - B[5] * B[7]);
  A[1] = inv_det * (B[2] * B[7] - B[1] * B[8]);
  A[2] = inv_det * (B[1] * B[5] - B[2] * B[4]);
  A[3] = inv_det * (B[5] * B[6] - B[3] * B[8]);
  A[4] = inv_det * (B[0] * B[8] - B[2] * B[6]);
  A[5] = inv_det * (B[2] * B[3] - B[0] * B[5]);
  A[6] = inv_det * (B[3] * B[7] - B[4] * B[6]);
  A[7] = inv_det * (B[1] * B[6] - B[0] * B[7]);
  A[8] = inv_det * (B[0] * B[4] - B[1] * B[3]);
}

}

template <typename VEC, typename SCALAR>
void delfem2::DiffDeformationGradientOfTet(
    SCALAR dF[4][3],
    const VEC &P0,
    const VEC &P1,
    const VEC &P2,
    const VEC &P3) {
  SCALAR Bi0[9] = {
      P1[0] - P0[0], P2[0] - P0[0], P3[0] - P0[0],
      P1[1] - P0[1], P2[1] - P0[1], P3[1] - P0[1],
      P1[2] - P0[2], P2[2] - P0[2], P3[2] - P0[2] };
  ::delfem2::geo_tet::Inverse_Mat3(Bi0);
  dF[0][0] = -Bi0[0] - Bi0[3] - Bi0[6];
  dF[0][1] = -Bi0[1] - Bi0[4] - Bi0[7];
  dF[0][2] = -Bi0[2] - Bi0[5] - Bi0[8];
  dF[1][0] = Bi0[0];
  dF[1][1] = Bi0[1];
  dF[1][2] = Bi0[2];
  dF[2][0] = Bi0[3];
  dF[2][1] = Bi0[4];
  dF[2][2] = Bi0[5];
  dF[3][0] = Bi0[6];
  dF[3][1] = Bi0[7];
  dF[3][2] = Bi0[8];
}

#endif // DFM2_GEOSOLIDELM_V3_H
