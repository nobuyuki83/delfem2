/*
 * Copyright (c) 2019 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */

#include "delfem2/point_on_surface_mesh.h"

#include <map>

#include "delfem2/vec3.h"

template<typename T>
std::array<T,3> delfem2::PointOnSurfaceMesh<T>::PositionOnMeshTri3(
    const std::vector<double> &vtx_xyz,
    const std::vector<unsigned int> &tri_vtx) const {
  assert(itri < tri_vtx.size() / 3);
  const unsigned int i0 = tri_vtx[itri * 3 + 0];
  const unsigned int i1 = tri_vtx[itri * 3 + 1];
  const unsigned int i2 = tri_vtx[itri * 3 + 2];
  const CVec3<T> p0(vtx_xyz.data() + i0 * 3);
  const CVec3<T> p1(vtx_xyz.data() + i1 * 3);
  const CVec3<T> p2(vtx_xyz.data() + i2 * 3);
  const CVec3<T> res = r0 * p0 + r1 * p1 + (1.0 - r0 - r1) * p2;
  return {res.x, res.y, res.z};
}
#ifdef DFM2_STATIC_LIBRARY
template std::array<double,3> delfem2::PointOnSurfaceMesh<double>::PositionOnMeshTri3(
    const std::vector<double> &,
    const std::vector<unsigned int> &) const;
#endif

// --------------------------------------------

template<typename T>
std::array<T,3> delfem2::PointOnSurfaceMesh<T>::PositionOnMeshTri3(
    const double *vtx_xyz,
    [[maybe_unused]] size_t num_vtx,
    const unsigned int *tri_vtx,
    [[maybe_unused]] size_t num_tri) const {
  assert(itri < num_tri);
  const unsigned int i0 = tri_vtx[itri * 3 + 0];
  const unsigned int i1 = tri_vtx[itri * 3 + 1];
  const unsigned int i2 = tri_vtx[itri * 3 + 2];
  const CVec3<T> p0(vtx_xyz+i0 * 3);
  const CVec3<T> p1(vtx_xyz+i1 * 3);
  const CVec3<T> p2(vtx_xyz+i2 * 3);
  const CVec3<T> q = r0 * p0 + r1 * p1 + (1.0 - r0 - r1) * p2;
  return {q.x, q.y, q.z};
}
#ifdef DFM2_STATIC_LIBRARY
template std::array<double,3> delfem2::PointOnSurfaceMesh<double>::PositionOnMeshTri3(
    const double *,
    size_t,
    const unsigned int *,
    size_t) const;
#endif

// -----------------------------------------------


template<typename T>
std::array<T,3> delfem2::PointOnSurfaceMesh<T>::UnitNormalOnMeshTri3(
    const std::vector<double> &vtx_normal,
    const std::vector<unsigned int> &tri_vtx) const {
  assert(itri < tri_vtx.size() / 3);
  const unsigned int i0 = tri_vtx[itri * 3 + 0];
  const unsigned int i1 = tri_vtx[itri * 3 + 1];
  const unsigned int i2 = tri_vtx[itri * 3 + 2];
  const CVec3<T> n0(vtx_normal.data()+i0 * 3);
  const CVec3<T> n1(vtx_normal.data()+i1 * 3);
  const CVec3<T> n2(vtx_normal.data()+i2 * 3);
  const CVec3<T> nr = (r0 * n0 + r1 * n1 + (1.0 - r0 - r1) * n2).normalized();
  return {nr.x, nr.y, nr.z};
}
#ifdef DFM2_STATIC_LIBRARY
template std::array<double,3> delfem2::PointOnSurfaceMesh<double>::UnitNormalOnMeshTri3(
    const std::vector<double> &,
    const std::vector<unsigned int> &) const;
#endif

// ------------------------------------------

template<typename T>
std::array<T,3> delfem2::PointOnSurfaceMesh<T>::UnitNormalOnMeshTri3(
    [[maybe_unused]] const double *vtx_xyz,
    [[maybe_unused]] unsigned int num_vtx,
    const unsigned int *tri_vtx,
    [[maybe_unused]] unsigned int num_tri,
    const double *vtx_norm) const {
  assert(itri < num_tri);
  const unsigned int i0 = tri_vtx[itri * 3 + 0];
  const unsigned int i1 = tri_vtx[itri * 3 + 1];
  const unsigned int i2 = tri_vtx[itri * 3 + 2];
  const CVec3<T> n0(vtx_norm+i0 * 3);
  const CVec3<T> n1(vtx_norm+i1 * 3);
  const CVec3<T> n2(vtx_norm+i2 * 3);
  const CVec3<T> nr = (r0 * n0 + r1 * n1 + (1.0 - r0 - r1) * n2).normalized();
  return {nr.x, nr.y, nr.z};
}
#ifdef DFM2_STATIC_LIBRARY
template std::array<double,3> delfem2::PointOnSurfaceMesh<double>::UnitNormalOnMeshTri3(
    const double *vtx_xyz,
    unsigned int num_vtx,
    const unsigned int *tri_vtx,
    unsigned int num_tri,
    const double *vtx_norm) const;
#endif

// ------------------------------------------

template<typename T>
std::array<T,3> delfem2::PointOnSurfaceMesh<T>::PositionOnFaceOfMeshTet(
    const std::vector<double> &aXYZ,
    const std::vector<unsigned int> &aTet,
    const std::vector<unsigned int> &aTetFace) const {
  const int noelTetFace[4][3] = {
      {1, 2, 3},
      {0, 3, 2},
      {0, 1, 3},
      {0, 2, 1}};
  int itet = aTetFace[itri * 2 + 0];
  int iface = aTetFace[itri * 2 + 1];
  double r2 = (1 - r0 - r1);
  int ielno0 = noelTetFace[iface][0];
  int ielno1 = noelTetFace[iface][1];
  int ielno2 = noelTetFace[iface][2];
  const unsigned int iq0 = aTet[itet * 4 + ielno0];
  const unsigned int iq1 = aTet[itet * 4 + ielno1];
  const unsigned int iq2 = aTet[itet * 4 + ielno2];
  CVec3<T> p;
  p.p[0] = r0 * aXYZ[iq0 * 3 + 0] + r1 * aXYZ[iq1 * 3 + 0] + r2 * aXYZ[iq2 * 3 + 0];
  p.p[1] = r0 * aXYZ[iq0 * 3 + 1] + r1 * aXYZ[iq1 * 3 + 1] + r2 * aXYZ[iq2 * 3 + 1];
  p.p[2] = r0 * aXYZ[iq0 * 3 + 2] + r1 * aXYZ[iq1 * 3 + 2] + r2 * aXYZ[iq2 * 3 + 2];
  return {p.x, p.y, p.z};
}
#ifdef DFM2_STATIC_LIBRARY
template std::array<double,3> delfem2::PointOnSurfaceMesh<double>::PositionOnFaceOfMeshTet(
    const std::vector<double> &aXYZ,
    const std::vector<unsigned int> &aTet,
    const std::vector<unsigned int> &aTetFace) const;
#endif

// ----------------------------------------

template<typename T>
std::array<T,3> delfem2::PointOnSurfaceMesh<T>::PositionOnGrid2(
    unsigned int nx,
    [[maybe_unused]] unsigned int ny,
    double el,
    std::vector<float> &aH) const {
  int iey = (this->itri / 2) / nx;
  int iex = (this->itri / 2) - nx * iey;
  CVec3<T> p00((iex + 0) * el, (iey + 0) * el, aH[(iey + 0) * nx + (iex + 0)]);
  CVec3<T> p10((iex + 1) * el, (iey + 0) * el, aH[(iey + 0) * nx + (iex + 1)]);
  CVec3<T> p01((iex + 0) * el, (iey + 1) * el, aH[(iey + 1) * nx + (iex + 0)]);
  CVec3<T> p11((iex + 1) * el, (iey + 1) * el, aH[(iey + 1) * nx + (iex + 1)]);
  if (this->itri % 2 == 0) {
    CVec3<T> pr = p00 * r0 + p10 * r1 + p11 * (1 - r0 - r1);
    return {pr.x, pr.y, pr.z};
  }
  CVec3<T> pr = p00 * r0 + p11 * r1 + p01 * (1 - r0 - r1);
  return {pr.x, pr.y, pr.z};
}
#ifdef DFM2_STATIC_LIBRARY
template std::array<double,3> delfem2::PointOnSurfaceMesh<double>::PositionOnGrid2(
    unsigned int nx, unsigned int ny,
    double el, std::vector<float> &aH) const;
#endif

// ----------------------------------------

template<typename T>
bool delfem2::PointOnSurfaceMesh<T>::Check(
    [[maybe_unused]] const std::vector<double> &vtx_xyz,
    const std::vector<unsigned int> &tri_vtx,
    double eps) const {
  if (itri >= tri_vtx.size() / 3) { return false; }
  if (r0 < -eps || r0 > 1 + eps) { return false; }
  if (r1 < -eps || r1 > 1 + eps) { return false; }
  double r2 = 1 - r0 - r1;
  if (r2 < -eps || r2 > 1 + eps) { return false; }
  return true;
}
#ifdef DFM2_STATIC_LIBRARY
template bool delfem2::PointOnSurfaceMesh<double>::Check(
    const std::vector<double> &vtx_xyz,
    const std::vector<unsigned int> &tri_vtx,
    double eps) const;
#endif
