//
// Created by Nobuyuki Umetani on 2022/02/05.
//


#include "delfem2/geo_tet.h"

#include "delfem2/mat3_funcs.h"
#include "delfem2/geo_tri.h"
#include "delfem2/vec3_funcs.h"

template<typename VEC, typename T>
T delfem2::Volume_Tet(
    const VEC &v0,
    const VEC &v1,
    const VEC &v2,
    const VEC &v3) {
  T v = (v1[0] - v0[0]) * ((v2[1] - v0[1]) * (v3[2] - v0[2]) - (v3[1] - v0[1]) * (v2[2] - v0[2]))
      + (v1[1] - v0[1]) * ((v2[2] - v0[2]) * (v3[0] - v0[0]) - (v3[2] - v0[2]) * (v2[0] - v0[0]))
      + (v1[2] - v0[2]) * ((v2[0] - v0[0]) * (v3[1] - v0[1]) - (v3[0] - v0[0]) * (v2[1] - v0[1]));
  return v * 0.16666666666666666666666666666667;
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
  ::delfem2::Inverse_Mat3(Bi0);
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

template<typename VEC, typename T>
std::array<T, 9> delfem2::DeformationGradientOfTet(
    const VEC &P0,
    const VEC &P1,
    const VEC &P2,
    const VEC &P3,
    const VEC &p0,
    const VEC &p1,
    const VEC &p2,
    const VEC &p3) {
  namespace dfm2 = delfem2;
  std::array<T, 9> Basis0 = delfem2::Mat3_3BasesOfTet<VEC,T>(P0,P1,P2,P3);
  const std::array<T, 9> basis0 = delfem2::Mat3_3BasesOfTet<VEC,T>(p0,p1,p2,p3);
  delfem2::Inverse_Mat3(Basis0.data());
  std::array<T, 9> r{};
  delfem2::MatMat3(r.data(), basis0.data(), Basis0.data());
  return r;
}

template<typename VEC, typename T>
double delfem2::Height_Tet(
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

template<typename VEC, typename T>
T delfem2::Volume_OrgTet(
    const VEC &v1,
    const VEC &v2,
    const VEC &v3) {
  T v = v1[0] * (v2[1] * v3[2] - v3[1] * v2[2])
      + v1[1] * (v2[2] * v3[0] - v3[2] * v2[0])
      + v1[2] * (v2[0] * v3[1] - v3[0] * v2[1]);
  return v * static_cast<T>(1.0 / 6.0);
}

template<typename VEC, typename T>
bool delfem2::barycentricCoord_Origin_Tet(
    T &r0,
    T &r1,
    T &r2,
    const VEC &p0,
    const VEC &p1,
    const VEC &p2,
    const VEC &p3) {
  const T v0 = Volume_OrgTet(p1, p2, p3);
  const T v1 = Volume_OrgTet(p2, p0, p3);
  const T v2 = Volume_OrgTet(p1, p3, p0);
  const T v3 = Volume_OrgTet(p1, p0, p2);
  const T vt_inv = 1 / (v0 + v1 + v2 + v3);
  r0 = v0 * vt_inv;
  r1 = v1 * vt_inv;
  r2 = v2 * vt_inv;
  return true;
}

template<typename VEC, typename T>
T delfem2::ShortestSquaredEdgeLengthOfTet(
    const VEC &ipo0,
    const VEC &ipo1,
    const VEC &ipo2,
    const VEC &ipo3) {
  T edge1, edge2;
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

template<typename VEC, typename T>
T delfem2::LongestSquaredEdgeLengthOfTet(
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

// ============================================

#ifdef DFM2_STATIC_LIBRARY

#include "delfem2/vec3.h"

namespace delfem2 {

using f0 = float [3];
using d0 = double [3];
using f1 = float *;
using d1 = double *;
using f2 = std::array<float, 3>;
using d2 = std::array<double, 3>;
using f3 = CVec3f;
using d3 = CVec3d;
//
template double Height_Tet(const d0 &, const d0 &, const d0 &, const d0 &);
template double Height_Tet(const d1 &, const d1 &, const d1 &, const d1 &);
template double Height_Tet(const d2 &, const d2 &, const d2 &, const d2 &);
template double Height_Tet(const d3 &, const d3 &, const d3 &, const d3 &);
template double Height_Tet(const f0 &, const f0 &, const f0 &, const f0 &);
template double Height_Tet(const f1 &, const f1 &, const f1 &, const f1 &);
template double Height_Tet(const f2 &, const f2 &, const f2 &, const f2 &);
template double Height_Tet(const f3 &, const f3 &, const f3 &, const f3 &);
//
template float Volume_OrgTet(const f0&, const f0&, const f0&);
template float Volume_OrgTet(const f1&, const f1&, const f1&);
template float Volume_OrgTet(const f2&, const f2&, const f2&);
template float Volume_OrgTet(const f3&, const f3&, const f3&);
template double Volume_OrgTet(const d0&, const d0&, const d0&);
template double Volume_OrgTet(const d1&, const d1&, const d1&);
template double Volume_OrgTet(const d2&, const d2&, const d2&);
template double Volume_OrgTet(const d3&, const d3&, const d3&);
//
template float Volume_Tet(const f0&, const f0&, const f0&, const f0&);
template float Volume_Tet(const f1&, const f1&, const f1&, const f1&);
template float Volume_Tet(const f2&, const f2&, const f2&, const f2&);
template float Volume_Tet(const f3&, const f3&, const f3&, const f3&);
template double Volume_Tet(const d0&, const d0&, const d0&, const d0&);
template double Volume_Tet(const d1&, const d1&, const d1&, const d1&);
template double Volume_Tet(const d2&, const d2&, const d2&, const d2&);
template double Volume_Tet(const d3&, const d3&, const d3&, const d3&);
//
template std::array<float, 9> DeformationGradientOfTet(const f0&, const f0&, const f0&, const f0&,
                                                       const f0&, const f0&, const f0&, const f0&);
template std::array<float, 9> DeformationGradientOfTet(const f1&, const f1&, const f1&, const f1&,
                                                       const f1&, const f1&, const f1&, const f1&);
template std::array<float, 9> DeformationGradientOfTet(const f2&, const f2&, const f2&, const f2&,
                                                       const f2&, const f2&, const f2&, const f2&);
template std::array<float, 9> DeformationGradientOfTet(const f3&, const f3&, const f3&, const f3&,
                                                       const f3&, const f3&, const f3&, const f3&);
template std::array<double, 9> DeformationGradientOfTet(const d0&, const d0&, const d0&, const d0&,
                                                        const d0&, const d0&, const d0&, const d0&);
template std::array<double, 9> DeformationGradientOfTet(const d1&, const d1&, const d1&, const d1&,
                                                        const d1&, const d1&, const d1&, const d1&);
template std::array<double, 9> DeformationGradientOfTet(const d2&, const d2&, const d2&, const d2&,
                                                        const d2&, const d2&, const d2&, const d2&);
template std::array<double, 9> DeformationGradientOfTet(const d3&, const d3&, const d3&, const d3&,
                                                        const d3&, const d3&, const d3&, const d3&);
//
template void DiffDeformationGradientOfTet(double [4][3], const d0 &, const d0 &, const d0 &, const d0 &);
template void DiffDeformationGradientOfTet(double [4][3], const d1 &, const d1 &, const d1 &, const d1 &);
template void DiffDeformationGradientOfTet(double [4][3], const d2 &, const d2 &, const d2 &, const d2 &);
template void DiffDeformationGradientOfTet(double [4][3], const d3 &, const d3 &, const d3 &, const d3 &);
}

#endif
