/*
 * Copyright (c) 2019 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */

#include "delfem2/srch_bruteforce.h"

#include <map>
#include <set>
#include <stack>

#include "delfem2/geo_plane.h"
#include "delfem2/geo_edge.h"
#include "delfem2/geo_tri.h"
#include "delfem2/vec3.h"
#include "delfem2/vec3_funcs.h"

namespace delfem2::srchuni {

//! Volume of a tetrahedra
template<typename T>
T Volume_Tet(
    const CVec3<T> &v0,
    const CVec3<T> &v1,
    const CVec3<T> &v2,
    const CVec3<T> &v3) {
  const T x0 = (v2.p[1] - v0.p[1]) * (v3.p[2] - v0.p[2]) - (v3.p[1] - v0.p[1]) * (v2.p[2] - v0.p[2]);
  const T y0 = (v2.p[2] - v0.p[2]) * (v3.p[0] - v0.p[0]) - (v3.p[2] - v0.p[2]) * (v2.p[0] - v0.p[0]);
  const T z0 = (v2.p[0] - v0.p[0]) * (v3.p[1] - v0.p[1]) - (v3.p[0] - v0.p[0]) * (v2.p[1] - v0.p[1]);
  const T v = (v1.p[0] - v0.p[0]) * x0 + (v1.p[1] - v0.p[1]) * y0 + (v1.p[2] - v0.p[2]) * z0;
  return v * static_cast<T>(0.16666666666666666666666666666667);
}

}

// ----------------------------------------------

template<typename T>
delfem2::CVec3<T> delfem2::PointInSolidMesh<T>::getPos_Tet(
    const std::vector<double> &aXYZ_,
    const std::vector<int> &aTet_) const {
  assert(ielem >= 0 && ielem < (int) aTet_.size() / 4);
  int ip0 = aTet_[ielem * 4 + 0];
  int ip1 = aTet_[ielem * 4 + 1];
  int ip2 = aTet_[ielem * 4 + 2];
  int ip3 = aTet_[ielem * 4 + 3];
  const CVec3<T> p0(aXYZ_[ip0 * 3 + 0], aXYZ_[ip0 * 3 + 1], aXYZ_[ip0 * 3 + 2]);
  const CVec3<T> p1(aXYZ_[ip1 * 3 + 0], aXYZ_[ip1 * 3 + 1], aXYZ_[ip1 * 3 + 2]);
  const CVec3<T> p2(aXYZ_[ip2 * 3 + 0], aXYZ_[ip2 * 3 + 1], aXYZ_[ip2 * 3 + 2]);
  const CVec3<T> p3(aXYZ_[ip3 * 3 + 0], aXYZ_[ip3 * 3 + 1], aXYZ_[ip3 * 3 + 2]);
  return r0 * p0 + r1 * p1 + r2 * p2 + (1.0 - r0 - r1 - r2) * p3;
}
#ifdef DFM2_STATIC_LIBRARY
template delfem2::CVec3d delfem2::PointInSolidMesh<double>::getPos_Tet
    (const std::vector<double> &aXYZ,
     const std::vector<int> &aTet) const;
#endif

// ----------------------------------------

template<typename T>
void delfem2::PointInSolidMesh<T>::setPos_Tet(
    int it0,
    const CVec3<T> &q,
    const std::vector<double> &aXYZ_,
    const std::vector<int> &aTet_) {
  namespace lcl = delfem2::srchuni;
  assert(it0 >= 0 && it0 < (int) aTet_.size() / 4);
  int ip0 = aTet_[it0 * 4 + 0];
  int ip1 = aTet_[it0 * 4 + 1];
  int ip2 = aTet_[it0 * 4 + 2];
  int ip3 = aTet_[it0 * 4 + 3];
  const CVec3<T> p0(aXYZ_[ip0 * 3 + 0], aXYZ_[ip0 * 3 + 1], aXYZ_[ip0 * 3 + 2]);
  const CVec3<T> p1(aXYZ_[ip1 * 3 + 0], aXYZ_[ip1 * 3 + 1], aXYZ_[ip1 * 3 + 2]);
  const CVec3<T> p2(aXYZ_[ip2 * 3 + 0], aXYZ_[ip2 * 3 + 1], aXYZ_[ip2 * 3 + 2]);
  const CVec3<T> p3(aXYZ_[ip3 * 3 + 0], aXYZ_[ip3 * 3 + 1], aXYZ_[ip3 * 3 + 2]);
  double v0 = lcl::Volume_Tet(q, p1, p2, p3);
  double v1 = lcl::Volume_Tet(p0, q, p2, p3);
  double v2 = lcl::Volume_Tet(p0, p1, q, p3);
  //    double v3 = volume_Tet(p0,p1,p2, q);
  double vt = lcl::Volume_Tet(p0, p1, p2, p3);
  this->ielem = it0;
  this->r0 = v0 / vt;
  this->r1 = v1 / vt;
  this->r2 = v2 / vt;
}
#ifdef DFM2_STATIC_LIBRARY
template void delfem2::PointInSolidMesh<double>::setPos_Tet
    (int, const CVec3d &,
     const std::vector<double> &,
     const std::vector<int> &);
#endif

// ---------------------------------------------------------------------------------

template<typename T>
std::vector<delfem2::PointOnSurfaceMesh<T>>
delfem2::IntersectionLine_MeshTri3(
    const CVec3<T> &org,
    const CVec3<T> &dir,
    const std::vector<unsigned int> &aTri_,
    const std::vector<T> &aXYZ_,
    T eps) {
  std::vector<PointOnSurfaceMesh<T>> aPES;
  for (unsigned int itri = 0; itri < aTri_.size() / 3; ++itri) {
    const unsigned int ip0 = aTri_[itri * 3 + 0];
    const unsigned int ip1 = aTri_[itri * 3 + 1];
    const unsigned int ip2 = aTri_[itri * 3 + 2];
    assert(ip0 < aXYZ_.size() / 3);
    assert(ip1 < aXYZ_.size() / 3);
    assert(ip2 < aXYZ_.size() / 3);
    const CVec3<T> p0(aXYZ_[ip0 * 3 + 0], aXYZ_[ip0 * 3 + 1], aXYZ_[ip0 * 3 + 2]);
    const CVec3<T> p1(aXYZ_[ip1 * 3 + 0], aXYZ_[ip1 * 3 + 1], aXYZ_[ip1 * 3 + 2]);
    const CVec3<T> p2(aXYZ_[ip2 * 3 + 0], aXYZ_[ip2 * 3 + 1], aXYZ_[ip2 * 3 + 2]);
    double r0, r1;
    bool res = delfem2::IntersectRay_Tri3(r0, r1,
                                          org, dir, p0, p1, p2,
                                          eps);
    if (!res) { continue; }
    aPES.emplace_back(itri, r0, r1);
  }
  return aPES;
}
#ifdef DFM2_STATIC_LIBRARY
template
std::vector<delfem2::PointOnSurfaceMesh<double>> delfem2::IntersectionLine_MeshTri3(
    const CVec3<double> &org,
    const CVec3<double> &dir,
    const std::vector<unsigned int> &aTri,
    const std::vector<double> &aXYZ,
    double eps);
#endif

// -------------------------------------

template<typename T>
void delfem2::IntersectionRay_MeshTri3(
    std::map<T, PointOnSurfaceMesh<T>> &mapDepthPES,
    const CVec3<T> &org,
    const CVec3<T> &dir,
    const std::vector<unsigned int> &aTri_,
    const std::vector<T> &aXYZ_,
    T eps) {
  const std::vector<PointOnSurfaceMesh<T>> aPES = IntersectionLine_MeshTri3(
      org, dir,
      aTri_, aXYZ_,
      eps);
  mapDepthPES.clear();
  for (auto pes: aPES) {
    CVec3<T> p0 = pes.PositionOnMeshTri3(aXYZ_, aTri_);
    double depth = (p0 - org).dot(dir);
    if (depth < 0) continue;
    mapDepthPES.insert(std::make_pair(depth, pes));
  }
}
#ifdef DFM2_STATIC_LIBRARY
template
void delfem2::IntersectionRay_MeshTri3(
    std::map<double, PointOnSurfaceMesh<double>> &mapDepthPES,
    const CVec3<double> &org, const CVec3<double> &dir,
    const std::vector<unsigned int> &aTri,
    const std::vector<double> &aXYZ,
    double eps);
#endif

// ----------------------

template<typename T>
void delfem2::IntersectionRay_MeshTri3DPart(
    std::map<T, PointOnSurfaceMesh<T>> &mapDepthPES,
    const CVec3<T> &org,
    const CVec3<T> &dir,
    const std::vector<unsigned int> &aTri_,
    const std::vector<T> &aXYZ_,
    const std::vector<unsigned int> &aIndTri,
    T eps) {
  mapDepthPES.clear();
  for (unsigned int itri: aIndTri) {
    const unsigned int ip0 = aTri_[itri * 3 + 0];
    const unsigned int ip1 = aTri_[itri * 3 + 1];
    const unsigned int ip2 = aTri_[itri * 3 + 2];
    assert(ip0 < aXYZ_.size() / 3);
    assert(ip1 < aXYZ_.size() / 3);
    assert(ip2 < aXYZ_.size() / 3);
    const CVec3<T> p0(aXYZ_.data()+ip0 * 3);
    const CVec3<T> p1(aXYZ_.data()+ip1 * 3);
    const CVec3<T> p2(aXYZ_.data()+ip2 * 3);
    double r0, r1;
    bool res = IntersectRay_Tri3(
        r0, r1,
        org, dir, p0, p1, p2, eps);
    if (!res) { continue; }
    double r2 = 1 - r0 - r1;
    CVec3<T> q0 = p0 * r0 + p1 * r1 + p2 * r2;
    double depth = (q0 - org).dot(dir) / dir.squaredNorm();
    if (depth < 0) continue;
    mapDepthPES.insert(std::make_pair(depth, PointOnSurfaceMesh<T>(itri, r0, r1)));
  }
}
#ifdef DFM2_STATIC_LIBRARY
template
void delfem2::IntersectionRay_MeshTri3DPart(
    std::map<double, PointOnSurfaceMesh<double>> &,
    const CVec3<double> &,
    const CVec3<double> &,
    const std::vector<unsigned int> &,
    const std::vector<double> &,
    const std::vector<unsigned int> &,
    double);
#endif

// -----------------------------------------------

DFM2_INLINE void delfem2::IntersectionLine_Hightfield(
    std::vector<PointOnSurfaceMesh<double>> &aPes,
    const double src[3],
    const double dir[3],
    unsigned int nx,
    unsigned int ny,
    const std::vector<float> &aH) {
  assert(aH.size() == nx * ny);
  for (unsigned int iey = 0; iey < ny - 1; ++iey) {
    for (unsigned int iex = 0; iex < nx - 1; ++iex) {
      const double h00 = aH[(iey + 0) * nx + (iex + 0)];
      const double h10 = aH[(iey + 0) * nx + (iex + 1)];
      const double h01 = aH[(iey + 1) * nx + (iex + 0)];
      const double h11 = aH[(iey + 1) * nx + (iex + 1)];
      const double p00[3] = {double(iex + 0), double(iey + 0), h00};
      const double p10[3] = {double(iex + 1), double(iey + 0), h10};
      const double p01[3] = {double(iex + 0), double(iey + 1), h01};
      const double p11[3] = {double(iex + 1), double(iey + 1), h11};
      double r0 = 1.0, r1 = 1.0;
      if (IntersectRay_Tri3(
          r0, r1,
          CVec3d(src), CVec3d(dir),
          CVec3d(p00), CVec3d(p10), CVec3d(p11), 1.0e-3)) {
        aPes.emplace_back((iey * nx + iex) * 2 + 0, r0, r1);
      }
      // ---------------------
      if (IntersectRay_Tri3(
          r0, r1,
          CVec3d(src), CVec3d(dir),
          CVec3d(p00), CVec3d(p11), CVec3d(p01), 1.0e-3)) {
        aPes.emplace_back((iey * nx + iex) * 2 + 1, r0, r1);
      }
    } // iex
  } // iey
}

// ----------------------------------------------------------

template<typename T>
delfem2::PointOnSurfaceMesh<T> delfem2::Nearest_Point_MeshTri3D(
    const CVec3<T> &point,
    const std::vector<double> &vtx_xyz,
    const std::vector<unsigned int> &tri_vtx) {
  PointOnSurfaceMesh<T> pes;
  double min_dist = -1;
  const size_t num_tri = tri_vtx.size() / 3;
  for (unsigned int it = 0; it < num_tri; ++it) {
    const unsigned int i0 = tri_vtx[it * 3 + 0];
    const unsigned int i1 = tri_vtx[it * 3 + 1];
    const unsigned int i2 = tri_vtx[it * 3 + 2];
    const CVec3<T> p0 = CVec3<T>(vtx_xyz.data() + i0 * 3) - point;
    const CVec3<T> p1 = CVec3<T>(vtx_xyz.data() + i1 * 3) - point;
    const CVec3<T> p2 = CVec3<T>(vtx_xyz.data() + i2 * 3) - point;
    double r0, r1;
    CVec3<T> p_min = Nearest_Origin3_Tri3(r0, r1, p0, p1, p2);
    double dist = p_min.squaredNorm();
    if (min_dist < 0 || dist < min_dist) {
      min_dist = dist;
      pes = PointOnSurfaceMesh<T>(it, r0, r1);
    }
  }
  assert(pes.itri != UINT_MAX);
  return pes;
}
#ifdef DFM2_STATIC_LIBRARY
template delfem2::PointOnSurfaceMesh<double> delfem2::Nearest_Point_MeshTri3D(
    const CVec3d &,
    const std::vector<double> &,
    const std::vector<unsigned int> &);
#endif

// ----------------------------------------------


template<typename T>
delfem2::PointOnSurfaceMesh<T> delfem2::Nearest_Point_MeshTri3DPart(
    const CVec3<T> &q,
    const std::vector<double> &aXYZ,
    const std::vector<unsigned int> &aTri,
    const std::vector<int> &aIndTri_Cand) {
  double min_dist = -1;
  PointOnSurfaceMesh<T> pes;
  for (int itri0: aIndTri_Cand) {
    const unsigned int i0 = aTri[itri0 * 3 + 0];
    const unsigned int i1 = aTri[itri0 * 3 + 1];
    const unsigned int i2 = aTri[itri0 * 3 + 2];
    const CVec3<T> p0(aXYZ[i0 * 3 + 0] - q.x, aXYZ[i0 * 3 + 1] - q.y, aXYZ[i0 * 3 + 2] - q.z);
    const CVec3<T> p1(aXYZ[i1 * 3 + 0] - q.x, aXYZ[i1 * 3 + 1] - q.y, aXYZ[i1 * 3 + 2] - q.z);
    const CVec3<T> p2(aXYZ[i2 * 3 + 0] - q.x, aXYZ[i2 * 3 + 1] - q.y, aXYZ[i2 * 3 + 2] - q.z);
    double r0, r1;
    CVec3<T> p_min = Nearest_Origin3_Tri3(r0, r1, p0, p1, p2);
    assert(r0 > -1.0e-10 && r1 > -1.0e-10 && (1 - r0 - r1) > -1.0e-10);
    double dist = p_min.squaredNorm();
    if (min_dist < 0 || dist < min_dist) {
      min_dist = dist;
      pes = PointOnSurfaceMesh<T>(itri0, r0, r1);
    }
  }
  return pes;
}
#ifdef DFM2_STATIC_LIBRARY
template delfem2::PointOnSurfaceMesh<double> delfem2::Nearest_Point_MeshTri3DPart(
    const CVec3d &,
    const std::vector<double> &,
    const std::vector<unsigned int> &,
    const std::vector<int> &);
#endif

// ----------------------------------------------------------------------------

template<typename T>
delfem2::PointInSolidMesh<T> delfem2::Nearest_Point_MeshTet3D(
    const CVec3<T> &q,
    const std::vector<double> &aXYZ,
    const std::vector<int> &aTet) {
  const double eps = 1.0e-4;
  const unsigned int ntet = aTet.size() / 4;
  for (unsigned int itet = 0; itet < ntet; ++itet) {
    PointInSolidMesh<T> pt;
    pt.setPos_Tet(itet, q, aXYZ, aTet);
    if (pt.isInside(-eps)) { return pt; }
  }
  return PointInSolidMesh<T>();
}

// ------------------------------

template<typename T>
delfem2::PointInSolidMesh<T> delfem2::Nearest_Point_MeshTet3D(
    const CVec3<T> &p,
    int itet_start, // starting triangle
    const std::vector<double> &aXYZ,
    const std::vector<int> &aTet,
    const std::vector<int> &aTetSurRel) {
  const double eps = 1.0e-4;
  int itet1 = itet_start;
  if (itet1 < 0 || itet1 >= (int) aTet.size() / 4) { return PointInSolidMesh<T>(); }
  for (int itr = 0; itr < 50; itr++) {
    if (itet1 == -1) return PointInSolidMesh<T>();
    int ip0 = aTet[itet1 * 4 + 0];
    int ip1 = aTet[itet1 * 4 + 1];
    int ip2 = aTet[itet1 * 4 + 2];
    int ip3 = aTet[itet1 * 4 + 3];
    const CVec3<T> p0(aXYZ[ip0 * 3 + 0], aXYZ[ip0 * 3 + 1], aXYZ[ip0 * 3 + 2]);
    const CVec3<T> p1(aXYZ[ip1 * 3 + 0], aXYZ[ip1 * 3 + 1], aXYZ[ip1 * 3 + 2]);
    const CVec3<T> p2(aXYZ[ip2 * 3 + 0], aXYZ[ip2 * 3 + 1], aXYZ[ip2 * 3 + 2]);
    const CVec3<T> p3(aXYZ[ip3 * 3 + 0], aXYZ[ip3 * 3 + 1], aXYZ[ip3 * 3 + 2]);
    double v0 = Volume_Tet(p, p1, p2, p3);
    double v1 = Volume_Tet(p0, p, p2, p3);
    double v2 = Volume_Tet(p0, p1, p, p3);
    double v3 = Volume_Tet(p0, p1, p2, p);
    double vt = (v0 + v1 + v2 + v3);
    if (v0 > -eps * vt && v1 > -eps * vt && v2 > -eps * vt && v3 > -eps * vt) {
      double r0 = v0 / (v0 + v1 + v2 + v3);
      double r1 = v1 / (v0 + v1 + v2 + v3);
      double r2 = v2 / (v0 + v1 + v2 + v3);
      PointInSolidMesh<T> pt(itet1, r0, r1, r2);
      return pt;
    }
    if (v0 < v1 && v0 < v2 && v0 < v3) { itet1 = aTetSurRel[itet1 * 8 + 0 * 2 + 0]; }
    else if (v1 < v0 && v1 < v2 && v1 < v3) { itet1 = aTetSurRel[itet1 * 8 + 1 * 2 + 0]; }
    else if (v2 < v0 && v2 < v1 && v2 < v3) { itet1 = aTetSurRel[itet1 * 8 + 2 * 2 + 0]; }
    else { itet1 = aTetSurRel[itet1 * 8 + 3 * 2 + 0]; }
  }
  return PointInSolidMesh<T>();
}

// ---------------------------------------------------------------

template<typename T>
delfem2::PointOnSurfaceMesh<T> delfem2::Nearest_Point_MeshTetFace3D(
    const CVec3<T> &p0,
    const std::vector<double> &aXYZ,
    const std::vector<int> &aTet,
    const std::vector<int> &aTetFaceSrf) {
  const int noelTetFace[4][3] = {
      {1, 2, 3},
      {0, 3, 2},
      {0, 1, 3},
      {0, 2, 1}};
  ////
  double dist_min = -1.0;
  int itf_min = -1;
  CVec3<T> p_min;
  for (size_t itf = 0; itf < aTetFaceSrf.size() / 2; ++itf) {
    int itet = aTetFaceSrf[itf * 2 + 0];
    int iface = aTetFaceSrf[itf * 2 + 1];
    const int i0 = aTet[itet * 4 + noelTetFace[iface][0]];
    const int i1 = aTet[itet * 4 + noelTetFace[iface][1]];
    const int i2 = aTet[itet * 4 + noelTetFace[iface][2]];
    CVec3<T> q0 = CVec3<T>(aXYZ[i0 * 3 + 0], aXYZ[i0 * 3 + 1], aXYZ[i0 * 3 + 2]) - p0;
    CVec3<T> q1 = CVec3<T>(aXYZ[i1 * 3 + 0], aXYZ[i1 * 3 + 1], aXYZ[i1 * 3 + 2]) - p0;
    CVec3<T> q2 = CVec3<T>(aXYZ[i2 * 3 + 0], aXYZ[i2 * 3 + 1], aXYZ[i2 * 3 + 2]) - p0;
    double r0, r1;
    CVec3<T> p2 = Nearest_Origin3_Tri3(r0, r1, q0, q1, q2);
    double dist = p2.norm();
    if (itf_min == -1 || dist < dist_min) {
      dist_min = dist;
      itf_min = itf;
      p_min = p2;
    }
  }
  assert(itf_min != -1);
  {
    int itet = aTetFaceSrf[itf_min * 2 + 0];
    int iface = aTetFaceSrf[itf_min * 2 + 1];
    const int i0 = aTet[itet * 4 + noelTetFace[iface][0]];
    const int i1 = aTet[itet * 4 + noelTetFace[iface][1]];
    const int i2 = aTet[itet * 4 + noelTetFace[iface][2]];
    CVec3<T> q0(aXYZ[i0 * 3 + 0], aXYZ[i0 * 3 + 1], aXYZ[i0 * 3 + 2]);
    CVec3<T> q1(aXYZ[i1 * 3 + 0], aXYZ[i1 * 3 + 1], aXYZ[i1 * 3 + 2]);
    CVec3<T> q2(aXYZ[i2 * 3 + 0], aXYZ[i2 * 3 + 1], aXYZ[i2 * 3 + 2]);
    double a0 = Area_Tri(p_min, q1, q2);
    double a1 = Area_Tri(p_min, q2, q0);
    double a2 = Area_Tri(p_min, q0, q1);
    double inva = 1.0 / (a0 + a1 + a2);
    a0 *= inva;
    a1 *= inva;
    a2 *= inva;
    PointOnSurfaceMesh<T> ptf;
    ptf.itri = itf_min;
    ptf.r0 = a0;
    ptf.r1 = a1;
    return ptf;
  }
}

// ----------------------------------------

template<typename T>
double delfem2::SDFNormal_NearestPoint(
    CVec3<T> &n0,
    const CVec3<T> &p0,
    const PointOnSurfaceMesh<T> &pes,
    const double *aXYZ, unsigned int nXYZ,
    const unsigned int *aTri, unsigned int nTri,
    const double *aNorm) {
  CVec3<T> q1 = pes.PositionOnMeshTri3(aXYZ, nXYZ, aTri, nTri);
  double dist = (q1 - p0).norm();
  CVec3<T> n1 = pes.UnitNormalOnMeshTri3(aXYZ, nXYZ, aTri, nTri, aNorm);
  if ((q1 - p0).dot(n1) > 0) {  //inside
    if (dist < 1.0e-6) { n0 = n1; }
    else { n0 = (q1 - p0).normalized(); }
    return dist;
  } else { // outside
    if (dist < 1.0e-6) { n0 = n1; }
    else { n0 = (p0 - q1).normalized(); }
    return -dist;
  }
}
#ifdef DFM2_STATIC_LIBRARY
template double delfem2::SDFNormal_NearestPoint
    (CVec3d &n0,
     const CVec3d &p0,
     const PointOnSurfaceMesh<double> &pes,
     const double *aXYZ, unsigned int nXYZ,
     const unsigned int *aTri, unsigned int nTri,
     const double *aNorm);
#endif

// ------------------------------------------------------

template<typename T>
double delfem2::SDFNormal_NearestPoint(
    CVec3<T> &n0,
    const CVec3<T> &p0,
    const PointOnSurfaceMesh<T> &pes,
    const std::vector<double> &aXYZ,
    const std::vector<unsigned int> &aTri,
    const std::vector<double> &aNorm) {
  CVec3<T> q1 = pes.PositionOnMeshTri3(aXYZ, aTri);
  double dist = (q1 - p0).norm();
  CVec3<T> n1 = pes.UnitNormalOnMeshTri3(aNorm, aTri);
  if ((q1 - p0).dot(n1) > 0) {  //inside
    if (dist < 1.0e-6) { n0 = n1; }
    else { n0 = (q1 - p0).normalized(); }
    return dist;
  } else { // outside
    if (dist < 1.0e-6) { n0 = n1; }
    else { n0 = (p0 - q1).normalized(); }
    return -dist;
  }
}
#ifdef DFM2_STATIC_LIBRARY
template double delfem2::SDFNormal_NearestPoint
    (CVec3d &n0,
     const CVec3d &p0,
     const PointOnSurfaceMesh<double> &pes,
     const std::vector<double> &aXYZ,
     const std::vector<unsigned int> &aTri,
     const std::vector<double> &aNorm);
#endif

// ---------------------------------------------

template<typename T>
double delfem2::DistanceToTri(
    PointOnSurfaceMesh<T> &pes,
    const CVec3<T> &p,
    unsigned int itri0,
    const std::vector<double> &aXYZ_,
    const std::vector<unsigned int> &aTri_) {
  const unsigned int i0 = aTri_[itri0 * 3 + 0];
  const unsigned int i1 = aTri_[itri0 * 3 + 1];
  const unsigned int i2 = aTri_[itri0 * 3 + 2];
  const CVec3<T> p0(aXYZ_[i0 * 3 + 0] - p.p[0], aXYZ_[i0 * 3 + 1] - p.p[1], aXYZ_[i0 * 3 + 2] - p.p[2]);
  const CVec3<T> p1(aXYZ_[i1 * 3 + 0] - p.p[0], aXYZ_[i1 * 3 + 1] - p.p[1], aXYZ_[i1 * 3 + 2] - p.p[2]);
  const CVec3<T> p2(aXYZ_[i2 * 3 + 0] - p.p[0], aXYZ_[i2 * 3 + 1] - p.p[1], aXYZ_[i2 * 3 + 2] - p.p[2]);
  double r0, r1;
  CVec3<T> p_min = Nearest_Origin3_Tri3(r0, r1, p0, p1, p2);
  assert(r0 > -1.0e-10 && r1 > -1.0e-10 && (1 - r0 - r1) > -1.0e-10);
  pes.itri = itri0;
  pes.r0 = r0;
  pes.r1 = r1;
  return p_min.norm();
}
#ifdef DFM2_STATIC_LIBRARY
template double delfem2::DistanceToTri(
     PointOnSurfaceMesh<double> &pes,
     const CVec3<double> &p,
     unsigned int itri0,
     const std::vector<double> &aXYZ,
     const std::vector<unsigned int> &aTri);
#endif

// -------------------------------------------

template<typename T>
double delfem2::DistanceToTri(
    PointOnSurfaceMesh<T> &pes,
    const CVec3<T> &p,
    unsigned int itri0,
    const double *aXYZ_,
    [[maybe_unused]] size_t nXYZ,
    const unsigned int *aTri_,
    [[maybe_unused]] size_t nTri) {
  const unsigned int i0 = aTri_[itri0 * 3 + 0];
  const unsigned int i1 = aTri_[itri0 * 3 + 1];
  const unsigned int i2 = aTri_[itri0 * 3 + 2];
  const CVec3<T> p0(aXYZ_[i0 * 3 + 0] - p.p[0], aXYZ_[i0 * 3 + 1] - p.p[1], aXYZ_[i0 * 3 + 2] - p.p[2]);
  const CVec3<T> p1(aXYZ_[i1 * 3 + 0] - p.p[0], aXYZ_[i1 * 3 + 1] - p.p[1], aXYZ_[i1 * 3 + 2] - p.p[2]);
  const CVec3<T> p2(aXYZ_[i2 * 3 + 0] - p.p[0], aXYZ_[i2 * 3 + 1] - p.p[1], aXYZ_[i2 * 3 + 2] - p.p[2]);
  double r0, r1;
  CVec3<T> p_min = Nearest_Origin3_Tri3(r0, r1, p0, p1, p2);
  assert(r0 > -1.0e-10 && r1 > -1.0e-10 && (1 - r0 - r1) > -1.0e-10);
  pes.itri = itri0;
  pes.r0 = r0;
  pes.r1 = r1;
  return p_min.norm();
}
#ifdef DFM2_STATIC_LIBRARY
template double delfem2::DistanceToTri(
    PointOnSurfaceMesh<double> &pes,
    const CVec3<double> &p,
    unsigned int itri0,
    const double *aXYZ,
    size_t nXYZ,
    const unsigned int *aTri,
    size_t nTri);
#endif

// ---------------------------------

std::vector<unsigned int> delfem2::IndexesOfConnectedTriangleInSphere(
  const std::array<double, 3> pos,
  double rad,
  unsigned int itri0,
  const std::vector<double> &vtx_xyz,
  const std::vector<unsigned int> &tri_vtx,
  const std::vector<unsigned int> &tri_adjtri) {
  std::vector<unsigned int> res;
  std::set<unsigned int> searched;
  std::stack<unsigned int> next0;
  next0.push(itri0);
  while (!next0.empty()) {
    unsigned int iel0 = next0.top();
    next0.pop();
    if (searched.find(iel0) != searched.end()) { continue; } // already studied
    searched.insert(iel0);
    double dist_min = -1;
    {
      double pn[3], r0, r1;
      delfem2::Nearest_Triangle3_Point3(
        pn, r0, r1,
        pos.data(),
        vtx_xyz.data() + tri_vtx[iel0 * 3 + 0] * 3,
        vtx_xyz.data() + tri_vtx[iel0 * 3 + 1] * 3,
        vtx_xyz.data() + tri_vtx[iel0 * 3 + 2] * 3);
      dist_min = delfem2::Distance3(pn,pos.data());
    }
    if (dist_min > rad) { continue; }
    res.push_back(iel0);
    for (unsigned int ie = 0; ie < 3; ++ie) {
      unsigned int iel1 = tri_adjtri[iel0 * 3 + ie];
      if (iel1 == UINT_MAX) { continue; }
      next0.push(iel1);
    }
  }
  return res;
}

bool delfem2::IsTherePointOnMeshInsideSphere(
  const std::tuple<unsigned int, double, double> &smpli,
  double rad,
  const std::vector<std::tuple<unsigned int, double, double> > &samples,
  const std::multimap<unsigned int, unsigned int> &el_smpl,
  const std::vector<double> &vtx_xyz,
  const std::vector<unsigned int> &tri_vtx,
  const std::vector<unsigned int> &tri_adjtri){
  namespace dfm2 = delfem2;
  const auto posi = delfem2::PointOnSurfaceMesh<double>(smpli).PositionOnMeshTri3(vtx_xyz, tri_vtx);
  std::vector<unsigned int> aIE = dfm2::IndexesOfConnectedTriangleInSphere(
    posi, rad,
    std::get<0>(smpli), vtx_xyz, tri_vtx, tri_adjtri);
  for (auto ie: aIE) {
    const auto[il, iu] = el_smpl.equal_range(ie);
    for (auto it = il; it != iu; ++it) {
      const unsigned int jsmpl = it->second;
      const auto smplj = samples[jsmpl];
      const auto posj = delfem2::PointOnSurfaceMesh<double>(smplj).PositionOnMeshTri3(vtx_xyz, tri_vtx);
      const double dist = dfm2::Distance3(posi.data(), posj.data());
      if (dist < rad) { return true; }
    }
  }
  return false;
}
