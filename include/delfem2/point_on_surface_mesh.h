/*
 * Copyright (c) 2019 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */

// DONE(2020/12/5): rename "CPointElemSolid" -> CPtElm3
// DONE(2020/12/5): rename "CPointElemSurf" -> CPtElm2

#ifndef DFM2_POINT_ON_SURFACE_MESH_H
#define DFM2_POINT_ON_SURFACE_MESH_H

#include <cstdio>
#include <climits>
#include <map>
#include <array>

#include "delfem2/vec3.h"
#include "delfem2/dfm2_inline.h"

// -----------------------------------

namespace delfem2 {

// class for point in a triangle
/**
 * @tparam T
 */
template<typename T>
class PointOnSurfaceMesh {
 public:
  PointOnSurfaceMesh()
      : itri(UINT_MAX), r0(0), r1(0) {}

  PointOnSurfaceMesh(unsigned int itri, double r0, double r1)
      : itri(itri), r0(r0), r1(r1) {}

  [[nodiscard]] std::array<T,3> PositionOnMeshTri3(
      const std::vector<double> &aXYZ,
      const std::vector<unsigned int> &aTri) const;

  std::array<T,3> PositionOnMeshTri3(
      const double *aXYZ,
      size_t num_vtx,
      const unsigned int *aTri,
      size_t num_tri) const;

  [[nodiscard]] std::array<T,3> PositionOnFaceOfMeshTet(
      const std::vector<double> &aXYZ,
      const std::vector<unsigned int> &aTet,
      const std::vector<unsigned int> &aTetFace) const;

  std::array<T,3> PositionOnGrid2(
      unsigned int nx,
      unsigned int ny,
      double el,
      std::vector<float> &aH) const;

  [[nodiscard]] std::array<T,3> UnitNormalOnMeshTri3(
      const std::vector<double> &aXYZ,
      const std::vector<unsigned int> &aTri,
      const std::vector<double> &aNorm) const;

  std::array<T,3> UnitNormalOnMeshTri3(
      const double *vtx_xyz,
      unsigned int num_vtx,
      const unsigned int *tri_vtx,
      unsigned int num_tri,
      const double *vtx_norm) const;

  [[nodiscard]] bool Check(
      const std::vector<double> &vtx_xyz,
      const std::vector<unsigned int> &tri_vtx,
      double eps) const;
 public:
  unsigned int itri; // can be UINT_MAX
  double r0, r1;
};

using PointOnSurfaceMeshd = PointOnSurfaceMesh<double>;
using PointOnSurfaceMeshf = PointOnSurfaceMesh<float>;

template<typename T>
std::ostream &operator<<(
    std::ostream &output,
    const PointOnSurfaceMesh<T> &v) {
  output.setf(std::ios::scientific);
  output << v.itri << " " << v.r0 << " " << v.r1;
  return output;
}

template<typename T>
std::istream &operator>>(
    std::istream &input,
    PointOnSurfaceMesh<T> &v) {
  input >> v.itri >> v.r0 >> v.r1;
  return input;
}

}

#ifndef DFM2_STATIC_LIBRARY
#  include "delfem2/point_on_surface_mesh.cpp"
#endif

#endif /* search_mesh_hpp */
