/*
 * Copyright (c) 2019 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */

// DONE(2020/12/5): rename "CPointElemSolid" -> CPtElm3
// DONE(2020/12/5): rename "CPointElemSurf" -> CPtElm2

#ifndef DFM2_SRCHUNI_V3_H
#define DFM2_SRCHUNI_V3_H

#include <cstdio>
#include <climits>
#include <map>
#include <array>

#include "delfem2/vec3.h"
#include "delfem2/dfm2_inline.h"
#include "delfem2/point_on_surface_mesh.h"

// -----------------------------------

namespace delfem2 {

/**
 * @tparam T
 */
template<typename T>
class PointInSolidMesh {
 public:
  PointInSolidMesh() : ielem(-1), r0(0), r1(0), r2(0) {}
  PointInSolidMesh(int itet, double r0, double r1, double r2) : ielem(itet), r0(r0), r1(r1), r2(r2) {}
  [[nodiscard]] bool isInside(double eps) const {
    double r3 = (1 - r0 - r1 - r2);
    if (r0 > eps && r1 > eps && r2 > eps && r3 > eps) { return true; }
    return false;
  }
  [[nodiscard]] int indFace(double eps) const {
    double r3 = (1 - r0 - r1 - r2);
    if (fabs(r0) < eps && fabs(r1) > eps && fabs(r2) > eps && fabs(r3) > eps) { return 0; }
    if (fabs(r0) > eps && fabs(r1) < eps && fabs(r2) > eps && fabs(r3) > eps) { return 1; }
    if (fabs(r0) > eps && fabs(r1) > eps && fabs(r2) < eps && fabs(r3) > eps) { return 2; }
    if (fabs(r0) > eps && fabs(r1) > eps && fabs(r2) > eps && fabs(r3) < eps) { return 3; }
    return -1;
  }
  delfem2::CVec3<T> getPos_Tet(
      const std::vector<double> &aXYZ,
      const std::vector<int> &aTet) const;
  void setPos_Tet(
      int it0, const delfem2::CVec3<T> &q,
      const std::vector<double> &aXYZ,
      const std::vector<int> &aTet);
 public:
  int ielem;
  double r0, r1, r2;
};

// ----------------------------------------------------------

template<typename REAL>
std::vector<PointOnSurfaceMesh<REAL>> IntersectionLine_MeshTri3(
    const delfem2::CVec3<REAL> &org, const delfem2::CVec3<REAL> &dir,
    const std::vector<unsigned int> &aTri,
    const std::vector<REAL> &aXYZ,
    REAL eps);

template<typename REAL>
void IntersectionRay_MeshTri3(
    std::map<REAL, PointOnSurfaceMesh<REAL>> &mapDepthPES,
    const delfem2::CVec3<REAL> &org,
    const delfem2::CVec3<REAL> &dir,
    const std::vector<unsigned int> &aTri,
    const std::vector<REAL> &aXYZ,
    REAL eps);

template<typename REAL>
void IntersectionRay_MeshTri3DPart(
    std::map<REAL, PointOnSurfaceMesh<REAL>> &mapDepthPES,
    const delfem2::CVec3<REAL> &org,
    const delfem2::CVec3<REAL> &dir,
    const std::vector<unsigned int> &aTri,
    const std::vector<REAL> &aXYZ,
    const std::vector<unsigned int> &aIndTri,
    REAL eps);

DFM2_INLINE void IntersectionLine_Hightfield(
    std::vector<PointOnSurfaceMesh<double>> &aPos,
    const double src[3],
    const double dir[3],
    unsigned int nx,
    unsigned int ny,
    const std::vector<float> &aH);

// above functions for ray interesection
// -----------------------------------------------------------
// below functions for nearest

template<typename T>
PointOnSurfaceMesh<T> Nearest_Point_MeshTri3D(
    const delfem2::CVec3<T> &point,
    const std::vector<double> &vtx_xyz,
    const std::vector<unsigned int> &tri_vtx);

template<typename T>
PointOnSurfaceMesh<T> Nearest_Point_MeshTri3DPart(
    const delfem2::CVec3<T> &q,
    const std::vector<double> &aXYZ,
    const std::vector<unsigned int> &aTri,
    const std::vector<int> &aIndTri_Cand);

template<typename T>
PointOnSurfaceMesh<T> Nearest_Point_MeshTetFace3D(
    const delfem2::CVec3<T> &p0,
    const std::vector<double> &aXYZ,
    const std::vector<int> &aTet,
    const std::vector<int> &aTetFaceSrf);

template<typename T>
PointInSolidMesh<T> Nearest_Point_MeshTet3D(
    const delfem2::CVec3<T> &q,
    const std::vector<double> &aXYZ,
    const std::vector<int> &aTet);

template<typename T>
PointInSolidMesh<T> Nearest_Point_MeshTet3D(
    const delfem2::CVec3<T> &p,
    int itet_start, // starting triangle
    const std::vector<double> &aXYZ,
    const std::vector<int> &aTet,
    const std::vector<int> &aTetSurRel);

template<typename T>
double SDFNormal_NearestPoint(
    delfem2::CVec3<T> &n0,
    const delfem2::CVec3<T> &p0,
    const PointOnSurfaceMesh<T> &pes,
    const std::vector<double> &aXYZ,
    const std::vector<unsigned int> &aTri,
    const std::vector<double> &aNorm);

template<typename T>
double SDFNormal_NearestPoint(
    delfem2::CVec3<T> &n0,
    const delfem2::CVec3<T> &p0,
    const PointOnSurfaceMesh<T> &pes,
    const double *aXYZ,
    unsigned int nXYZ,
    const unsigned int *aTri,
    unsigned int nTri,
    const double *aNorm);

template<typename T>
double DistanceToTri(
    PointOnSurfaceMesh<T> &pes,
    const delfem2::CVec3<T> &p,
    unsigned int itri0,
    const std::vector<double> &aXYZ,
    const std::vector<unsigned int> &aTri);

template<typename T>
double DistanceToTri(
    PointOnSurfaceMesh<T> &pes,
    const delfem2::CVec3<T> &p,
    unsigned int itri0,
    const double *aXYZ,
    size_t nXYZ,
    const unsigned int *aTri,
    size_t nTri);

}

#ifndef DFM2_STATIC_LIBRARY
#  include "delfem2/srchuni_v3.cpp"
#endif

#endif /* DFM2_SRCHUNI_V3_H */
