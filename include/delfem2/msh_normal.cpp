/*
 * Copyright (c) 2019 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */

#include "delfem2/msh_normal.h"

#include <cassert>
#include <cmath>
#include <vector>
#include <functional>
#include <array>

#ifndef M_PI
#  define M_PI 3.14159265358979323846
#endif

// ------------------------------------------------

namespace delfem2::msh_normal {

template <typename T>
DFM2_INLINE T Length3(const T p[3]) {
  return std::sqrt(p[0] * p[0] + p[1] * p[1] + p[2] * p[2]);
}

template<typename T>
DFM2_INLINE void UnitNormalAreaTri3(
    T n[3],
    T &a,
    const T v1[3],
    const T v2[3],
    const T v3[3]) {
  n[0] = (v2[1] - v1[1]) * (v3[2] - v1[2]) - (v3[1] - v1[1]) * (v2[2] - v1[2]);
  n[1] = (v2[2] - v1[2]) * (v3[0] - v1[0]) - (v3[2] - v1[2]) * (v2[0] - v1[0]);
  n[2] = (v2[0] - v1[0]) * (v3[1] - v1[1]) - (v3[0] - v1[0]) * (v2[1] - v1[1]);
  a = std::sqrt(n[0] * n[0] + n[1] * n[1] + n[2] * n[2]) / 2;
  const T invlen = 1 / (a * 2);
  n[0] *= invlen;
  n[1] *= invlen;
  n[2] *= invlen;
}

template<typename T>
std::array<T,3> Normal_Tri3(
    const T v1[3],
    const T v2[3],
    const T v3[3]) {
  return {
      (v2[1] - v1[1]) * (v3[2] - v1[2]) - (v3[1] - v1[1]) * (v2[2] - v1[2]),
      (v2[2] - v1[2]) * (v3[0] - v1[0]) - (v3[2] - v1[2]) * (v2[0] - v1[0]),
      (v2[0] - v1[0]) * (v3[1] - v1[1]) - (v3[0] - v1[0]) * (v2[1] - v1[1]) };
}

}

// static function above
// ==============================================
// exposed function below

template<typename REAL>
void delfem2::Normal_MeshTri3D(
    REAL *vtx_normal,
    const REAL *vtx_xyz,
    size_t num_vtx,
    const unsigned int *tri_vtx,
    size_t num_tri) {
  for (unsigned int i = 0; i < num_vtx * 3; i++) { vtx_normal[i] = 0; }
  for (unsigned int itri = 0; itri < num_tri; itri++) {
    const unsigned int i0 = tri_vtx[itri * 3 + 0];
    assert(i0 < num_vtx);
    const unsigned int i1 = tri_vtx[itri * 3 + 1];
    assert(i1 < num_vtx);
    const unsigned int i2 = tri_vtx[itri * 3 + 2];
    assert(i2 < num_vtx);
    const REAL *p0 = vtx_xyz + i0 * 3;
    const REAL *p1 = vtx_xyz + i1 * 3;
    const REAL *p2 = vtx_xyz + i2 * 3;
    REAL un[3], area;
    msh_normal::UnitNormalAreaTri3(un, area, p0, p1, p2);
    vtx_normal[i0 * 3 + 0] += un[0];
    vtx_normal[i0 * 3 + 1] += un[1];
    vtx_normal[i0 * 3 + 2] += un[2];
    vtx_normal[i1 * 3 + 0] += un[0];
    vtx_normal[i1 * 3 + 1] += un[1];
    vtx_normal[i1 * 3 + 2] += un[2];
    vtx_normal[i2 * 3 + 0] += un[0];
    vtx_normal[i2 * 3 + 1] += un[1];
    vtx_normal[i2 * 3 + 2] += un[2];
  }
  for (unsigned int ino = 0; ino < num_vtx; ino++) {
    const REAL invlen = 1 / msh_normal::Length3(vtx_normal + ino * 3);
    vtx_normal[ino * 3 + 0] *= invlen;
    vtx_normal[ino * 3 + 1] *= invlen;
    vtx_normal[ino * 3 + 2] *= invlen;
  }
}
#ifdef DFM2_STATIC_LIBRARY
template void delfem2::Normal_MeshTri3D(
    float *,
    const float *,
    size_t,
    const unsigned int *,
    size_t);
template void delfem2::Normal_MeshTri3D(
    double *,
    const double *,
    size_t,
    const unsigned int *,
    size_t);
#endif

// --------------------------------

template<typename REAL>
std::array<REAL,3> delfem2::Normal_TriInMeshTri3(
    unsigned int itri,
    const REAL *vtx_xyz,
    const unsigned int *tri_vtx) {
  return msh_normal::Normal_Tri3(
      vtx_xyz + tri_vtx[itri * 3 + 0] * 3,
      vtx_xyz + tri_vtx[itri * 3 + 1] * 3,
      vtx_xyz + tri_vtx[itri * 3 + 2] * 3);
}

// ---------------------------------------

template<typename REAL>
void delfem2::Normal_MeshQuad3(
    std::vector<REAL> &aNorm,
    const std::vector<REAL> &aXYZ,
    const std::vector<unsigned int> &aQuad) {
  const size_t nXYZ = aXYZ.size() / 3;
  const size_t nQuad = aQuad.size() / 4;
  aNorm.resize(nXYZ * 3);
  // -------
  for (unsigned int i = 0; i < nXYZ * 3; i++) { aNorm[i] = 0; }
  for (unsigned int iquad = 0; iquad < nQuad; ++iquad) {
    const unsigned int i0 = aQuad[iquad * 4 + 0];
    const unsigned int i1 = aQuad[iquad * 4 + 1];
    const unsigned int i2 = aQuad[iquad * 4 + 2];
    const unsigned int i3 = aQuad[iquad * 4 + 3];
    const REAL *p0 = aXYZ.data() + i0 * 3;
    const REAL *p1 = aXYZ.data() + i1 * 3;
    const REAL *p2 = aXYZ.data() + i2 * 3;
    const REAL *p3 = aXYZ.data() + i3 * 3;
    REAL un0[3], a0;
    msh_normal::UnitNormalAreaTri3(un0, a0, p3, p0, p1);
    REAL un1[3], a1;
    msh_normal::UnitNormalAreaTri3(un1, a1, p0, p1, p2);
    REAL un2[3], a2;
    msh_normal::UnitNormalAreaTri3(un2, a2, p1, p2, p3);
    REAL un3[3], a3;
    msh_normal::UnitNormalAreaTri3(un3, a3, p2, p3, p0);
    aNorm[i0 * 3 + 0] += un0[0];
    aNorm[i0 * 3 + 1] += un0[1];
    aNorm[i0 * 3 + 2] += un0[2];
    aNorm[i1 * 3 + 0] += un1[0];
    aNorm[i1 * 3 + 1] += un1[1];
    aNorm[i1 * 3 + 2] += un1[2];
    aNorm[i2 * 3 + 0] += un2[0];
    aNorm[i2 * 3 + 1] += un2[1];
    aNorm[i2 * 3 + 2] += un2[2];
    aNorm[i3 * 3 + 0] += un3[0];
    aNorm[i3 * 3 + 1] += un3[1];
    aNorm[i3 * 3 + 2] += un3[2];
  }
  for (unsigned int ino = 0; ino < nXYZ; ino++) {
    const REAL n[3] = {aNorm[ino * 3 + 0], aNorm[ino * 3 + 1], aNorm[ino * 3 + 2]};
    const REAL invlen = 1 / msh_normal::Length3(n);
    aNorm[ino * 3 + 0] *= invlen;
    aNorm[ino * 3 + 1] *= invlen;
    aNorm[ino * 3 + 2] *= invlen;
  }
}
#ifdef DFM2_STATIC_LIBRARY
template void delfem2::Normal_MeshQuad3(
    std::vector<float> &,
    const std::vector<float> &,
    const std::vector<unsigned int> &);
template void delfem2::Normal_MeshQuad3(
    std::vector<double> &,
    const std::vector<double> &,
    const std::vector<unsigned int> &);
#endif

// ======================

#ifdef DFM2_STATIC_LIBRARY

namespace delfem2 {

template std::array<float,3> Normal_TriInMeshTri3(
    unsigned int itri,
    const float *vtx_xyz,
    const unsigned int *tri_vtx);
template std::array<double,3> Normal_TriInMeshTri3(
    unsigned int itri,
    const double *vtx_xyz,
    const unsigned int *tri_vtx);

}

#endif
