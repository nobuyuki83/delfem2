/*
 * Copyright (c) 2019 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */

// TOOD: split the code w.r.t the dependencies.
// there are too many dependencies here. code is too long.
// name the file properly

#ifndef DFM2_SRCH_V3BVHMSHTOPO_H
#define DFM2_SRCH_V3BVHMSHTOPO_H

#include <cmath>
#include <vector>

#include "delfem2/srchuni_v3.h" // CPointElemSurf
#include "delfem2/srchbvh.h"
#include "delfem2/mshuni.h" // sourrounding relationship
#include "delfem2/vec3.h"
#include "delfem2/mat4.h"
#include "delfem2/mshmisc.h"
#include "delfem2/points.h"

namespace delfem2 {

template <class BV>
void ConstructBVHTriangleMeshMortonCode(
    std::vector<delfem2::CNodeBVH2> &aNodeBVH,
    std::vector<BV> &aAABB,
    const std::vector<double> &aXYZ,
    const std::vector<unsigned int> &aTri)
{
  namespace dfm2 = delfem2;
  std::vector<double> vec_center_of_tri;
  dfm2::CentsMaxRad_MeshTri3(
      vec_center_of_tri,
      aXYZ, aTri);
  double min_xyz[3], max_xyz[3];
  delfem2::BoundingBox3_Points3(
      min_xyz, max_xyz,
      vec_center_of_tri.data(),
      static_cast<unsigned int>(vec_center_of_tri.size() / 3));
  std::vector<unsigned int> aSortedId;
  std::vector<std::uint32_t> aSortedMc;
  dfm2::SortedMortenCode_Points3(
      aSortedId, aSortedMc,
      vec_center_of_tri, min_xyz, max_xyz);
  dfm2::BVHTopology_Morton(
      aNodeBVH,
      aSortedId, aSortedMc);
#ifndef NDEBUG
  dfm2::Check_MortonCode_Sort(
      aSortedId, aSortedMc, vec_center_of_tri,
      min_xyz, max_xyz);
  dfm2::Check_MortonCode_RangeSplit(
      aSortedMc);
#endif
  dfm2::CLeafVolumeMaker_Mesh<BV, double> lvm(
      1.0e-10,
      aXYZ.data(), aXYZ.size() / 3,
      aTri.data(), aTri.size() / 3, 3);
  dfm2::BVH_BuildBVHGeometry(
      aAABB,
      0, aNodeBVH,
      lvm);
#ifndef NDEBUG
  dfm2::Check_BVH(aNodeBVH, vec_center_of_tri.size() / 3);
#endif
}

/**
 * @brief potential maximum distance of the nearest point
 */
template<typename BV, typename REAL>
void BVH_NearestPoint_MeshTri3D(
    double &dist_min,
    CPtElm2<REAL> &pes,
    //
    double px, double py, double pz,
    const std::vector<double> &aXYZ,
    const std::vector<unsigned int> &aTri,
    int ibvh,
    const std::vector<delfem2::CNodeBVH2> &aBVH,
    const std::vector<BV> &aBB) {
  double min0, max0;
  aBB[ibvh].Range_DistToPoint(min0, max0, px, py, pz);
  assert(min0 >= 0 && max0 >= min0);
  //
  if (dist_min >= 0 && min0 > dist_min) { return; }
  const unsigned int ichild0 = aBVH[ibvh].ichild[0];
  const unsigned int ichild1 = aBVH[ibvh].ichild[1];
  if (ichild1 == UINT_MAX) { // leaf
    const unsigned int itri0 = ichild0;
    CPtElm2<REAL> pes_tmp;
    double dist = DistanceToTri(
        pes_tmp,
        CVec3<REAL>(px, py, pz),
        itri0, aXYZ, aTri);
    if (dist_min < 0 || dist < dist_min) {
      dist_min = dist;
      pes = pes_tmp;
    }
    return;
  }
  //
  BVH_NearestPoint_MeshTri3D(dist_min, pes, px, py, pz, aXYZ, aTri, ichild0, aBVH, aBB);
  BVH_NearestPoint_MeshTri3D(dist_min, pes, px, py, pz, aXYZ, aTri, ichild1, aBVH, aBB);
}

// potential maximum distance of the nearest point
template<typename BV, typename REAL>
void BVH_NearestPoint_IncludedInBVH_MeshTri3D(
    double &dist_tri, // minimum distance to triangle
    double &dist_bv, // minimum distance to leaf bounding volume
    CPtElm2<REAL> &pes,
    //
    double px,
    double py,
    double pz,
    double rad_exp, // exploring distance
    const double *aXYZ,
    size_t nXYZ,
    const unsigned int *aTri,
    size_t nTri,
    unsigned int ibvh,
    const std::vector<delfem2::CNodeBVH2> &aBVH,
    const std::vector<BV> &aBB) {
  assert(ibvh < aBB.size());
  double min0, max0;
  aBB[ibvh].Range_DistToPoint(min0, max0,
                              px, py, pz);
  if (min0 > rad_exp) { return; }
  //
  const unsigned int ichild0 = aBVH[ibvh].ichild[0];
  const unsigned int ichild1 = aBVH[ibvh].ichild[1];
  if (ichild1 == UINT_MAX) { // leaf
    if (min0 < dist_bv) { dist_bv = min0; }
    if (min0 == 0.0) {
      dist_bv = 0.0;
      const unsigned int itri0 = ichild0;
      CPtElm2<REAL> pes_tmp;
      const double dist0 = DistanceToTri(
          pes_tmp,
          //
          CVec3<REAL>(px, py, pz),
          itri0, aXYZ, nXYZ, aTri, nTri);
      if (dist_tri < 0 || dist0 < dist_tri) {
        dist_tri = dist0;
        pes = pes_tmp;
      }
    }
    return;
  }
  BVH_NearestPoint_IncludedInBVH_MeshTri3D(
      dist_tri, dist_bv, pes,
      //
      px, py, pz, rad_exp,
      aXYZ, nXYZ, aTri, nTri,
      ichild0, aBVH, aBB);
  BVH_NearestPoint_IncludedInBVH_MeshTri3D(
      dist_tri, dist_bv, pes,
      //
      px, py, pz, rad_exp,
      aXYZ, nXYZ, aTri, nTri,
      ichild1, aBVH, aBB);
}

template<typename BV, typename REAL>
class CBVH_MeshTri3D {
 public:
  CBVH_MeshTri3D() :
      iroot_bvh(0) {}
  void Init(
      const double *pXYZ,
      size_t nXYZ,
      const unsigned int *pTri,
      size_t nTri,
      double margin) {
    assert(margin >= 0);
    { // make BVH topology
      std::vector<double> aElemCenter(nTri * 3);
      for (unsigned int itri = 0; itri < nTri; ++itri) {
        const unsigned int i0 = pTri[itri * 3 + 0];
        const unsigned int i1 = pTri[itri * 3 + 1];
        const unsigned int i2 = pTri[itri * 3 + 2];
        double x0 = (pXYZ[i0 * 3 + 0] + pXYZ[i1 * 3 + 0] + pXYZ[i2 * 3 + 0]) / 3.0;
        double y0 = (pXYZ[i0 * 3 + 1] + pXYZ[i1 * 3 + 1] + pXYZ[i2 * 3 + 1]) / 3.0;
        double z0 = (pXYZ[i0 * 3 + 2] + pXYZ[i1 * 3 + 2] + pXYZ[i2 * 3 + 2]) / 3.0;
        aElemCenter[itri * 3 + 0] = x0;
        aElemCenter[itri * 3 + 1] = y0;
        aElemCenter[itri * 3 + 2] = z0;
      }
      {
        std::vector<unsigned int> aTriSuTri;
        ElSuEl_MeshElem(aTriSuTri,
                        pTri, nTri,
                        delfem2::MESHELEM_TRI, nXYZ);
        iroot_bvh = BVHTopology_TopDown_MeshElem(aNodeBVH,
                                                 3, aTriSuTri,
                                                 aElemCenter);
      }
    }
    this->UpdateGeometry(pXYZ, nXYZ,
                         pTri, nTri,
                         margin);
    assert(aBB_BVH.size() == aNodeBVH.size());
  }
  void UpdateGeometry(
      const double *pXYZ,
      size_t nXYZ,
      const unsigned int *pTri,
      size_t nTri,
      double margin) {
    assert(margin >= 0);
    CLeafVolumeMaker_Mesh<BV, REAL> lvm(
        margin,
        pXYZ, nXYZ,
        pTri, nTri, 3);
    BVH_BuildBVHGeometry(
        aBB_BVH,
        //
        iroot_bvh, aNodeBVH,
        lvm);
    assert(aBB_BVH.size() == aNodeBVH.size());
  }
  double Nearest_Point_IncludedInBVH(
      CPtElm2<REAL> &pes,
      const CVec3<REAL> &p0,
      double rad_exp, // look leaf inside this radius
      const double *aXYZ0,
      size_t nXYZ,
      const unsigned int *aTri,
      size_t nTri) const {
    assert(aBB_BVH.size() == aNodeBVH.size());
    double dist = -1, dist_min = rad_exp;
    pes.itri = UINT_MAX;
    delfem2::BVH_NearestPoint_IncludedInBVH_MeshTri3D(
        dist, dist_min, pes,
        p0.x, p0.y, p0.z, rad_exp,
        aXYZ0, nXYZ, aTri, nTri,
        iroot_bvh, aNodeBVH, aBB_BVH);
    if (pes.itri == UINT_MAX) { return dist_min; }
    return dist;
  }
  CPtElm2<REAL> NearestPoint_Global(
      const CVec3<REAL> &p0,
      const std::vector<double> &aXYZ,
      const std::vector<unsigned int> &aTri) const {
    assert(aBB_BVH.size() == aNodeBVH.size());
    CPtElm2<REAL> pes;
    double dist_min = -1;
    BVH_NearestPoint_MeshTri3D(dist_min, pes,
                               p0.x, p0.y, p0.z,
                               aXYZ, aTri, iroot_bvh, aNodeBVH, aBB_BVH);
    return pes;
  }
  // inside positive
  double SignedDistanceFunction(
      CVec3<REAL> &n0,
      //
      const CVec3<REAL> &p0,
      const std::vector<double> &aXYZ0,
      const std::vector<unsigned int> &aTri0,
      const std::vector<double> &aNorm) const {
    assert(aBB_BVH.size() == aNodeBVH.size());
    CPtElm2<REAL> pes;
    {
      double dist_min = -1;
      delfem2::BVH_NearestPoint_MeshTri3D(dist_min, pes,
                                          p0.x, p0.y, p0.z,
                                          aXYZ0, aTri0,
                                          iroot_bvh, aNodeBVH, aBB_BVH);
    }
    const CVec3<REAL> q0 = pes.Pos_Tri(aXYZ0, aTri0);
    double dist = (q0 - p0).norm();
    if (!aBB_BVH[iroot_bvh].isInclude_Point(p0.x, p0.y, p0.z)) { // outside
      n0 = (p0 - q0).normalized();
      return -dist;
    }
    const CVec3<REAL> n1 = pes.UNorm_Tri(aXYZ0, aTri0, aNorm);
    if (dist < 1.0e-6) {
      n0 = n1;
      if ((q0 - p0).dot(n1) > 0) { return dist; } //inside
      return -dist; // outside
    }
    CVec3<REAL> dir = (CG_Tri3(pes.itri, aTri0, aXYZ0) - p0).normalized();
    if ((q0 - p0).dot(n1) < 0) { dir = -dir; } // probaby outside so shoot ray outside
    std::vector<unsigned int> aIndElem;
    BVH_GetIndElem_Predicate(aIndElem,
                             CIsBV_IntersectRay<BV>(p0.p, dir.p),
                             iroot_bvh, aNodeBVH, aBB_BVH);
    std::map<double, CPtElm2<REAL>> mapDepthPES1;
    IntersectionRay_MeshTri3DPart(mapDepthPES1,
                                  p0, dir,
                                  aTri0, aXYZ0, aIndElem,
                                  0.0);
    if (mapDepthPES1.size() % 2 == 0) { // outside
      n0 = (p0 - q0).normalized();
      return -dist;
    }
    n0 = (q0 - p0).normalized();
    return +dist;
  }
 public:
  int iroot_bvh;
  std::vector<delfem2::CNodeBVH2> aNodeBVH; // array of BVH node
  std::vector<BV> aBB_BVH;
};

template<typename T, typename REAL>
void Project_PointsIncludedInBVH_Outside(
    std::vector<double> &aXYZt,
    double cc,
    const CBVH_MeshTri3D<T, REAL> &bvh,
    const std::vector<double> &aXYZ0,
    const std::vector<unsigned int> &aTri0,
    const std::vector<double> &aNorm0) {
  for (unsigned int ip = 0; ip < aXYZt.size() / 3; ++ip) {
    CVec3<REAL> p0(aXYZt[ip * 3 + 0], aXYZt[ip * 3 + 1], aXYZt[ip * 3 + 2]);
    CPtElm2<REAL> pes;
    bvh.Nearest_Point_IncludedInBVH(pes,
                                    p0, 0.0,
                                    aXYZ0.data(), aXYZ0.size() / 3,
                                    aTri0.data(), aTri0.size() / 3);
    if (pes.itri == UINT_MAX) { continue; }
    CVec3<REAL> n0;
    double sdf = SDFNormal_NearestPoint(n0,
                                        p0, pes, aXYZ0, aTri0, aNorm0);
    if (sdf + cc < 0) continue;
    aXYZt[ip * 3 + 0] += (sdf + cc) * n0.x();
    aXYZt[ip * 3 + 1] += (sdf + cc) * n0.y();
    aXYZt[ip * 3 + 2] += (sdf + cc) * n0.z();
  }
}

template<typename REAL>
class CInfoNearest {
 public:
  CInfoNearest() :
      sdf(0.0), is_active(false) {}

 public:
  CPtElm2<REAL> pes;
  CVec3<REAL> pos;
  double sdf;
  bool is_active;
};

/**
 *
 * @tparam BV
 * @tparam REAL
 * @param aXYZt
 * @param aInfoNearest
 * @param nXYZt
 * @param cc (in) contact clearance
 * @param bvh
 * @param pXYZ0
 * @param nXYZ0
 * @param pTri0
 * @param nTri0
 * @param pNorm0
 * @param rad_explore
 */
template<typename BV, typename REAL>
void Project_PointsIncludedInBVH_Outside_Cache(
    double *aXYZt,
    std::vector<CInfoNearest<REAL>> &aInfoNearest,
    unsigned int nXYZt,
    double cc,
    const CBVH_MeshTri3D<BV, REAL> &bvh,
    const double *pXYZ0, unsigned int nXYZ0,
    const unsigned int *pTri0, unsigned int nTri0,
    const double *pNorm0,
    double rad_explore) {
  const unsigned int np = nXYZt;
  aInfoNearest.resize(np);
  for (unsigned int ip = 0; ip < np; ++ip) {
    CVec3<REAL> p0(aXYZt[ip * 3 + 0], aXYZt[ip * 3 + 1], aXYZt[ip * 3 + 2]);
    if (aInfoNearest[ip].is_active) {
      double dp = Distance(p0, aInfoNearest[ip].pos);
      if (aInfoNearest[ip].sdf + dp + cc < 0) {
        continue;
      }
    }
    aInfoNearest[ip].pos = p0;
    double dist0 = bvh.Nearest_Point_IncludedInBVH(aInfoNearest[ip].pes,
                                                   aInfoNearest[ip].pos, rad_explore,
                                                   pXYZ0, nXYZ0,
                                                   pTri0, nTri0);
    if (aInfoNearest[ip].pes.itri == UINT_MAX) {
      if (aInfoNearest[ip].is_active) {
        if (aInfoNearest[ip].sdf < 0) { aInfoNearest[ip].sdf = -dist0; }
        else { aInfoNearest[ip].sdf = +dist0; }
      } else {
        aInfoNearest[ip].sdf = -dist0;
        aInfoNearest[ip].is_active = true;
      }
      continue;
    }
    CVec3<REAL> n0;
    double sdf = SDFNormal_NearestPoint(
        n0,
        aInfoNearest[ip].pos,
        aInfoNearest[ip].pes,
        pXYZ0, nXYZ0,
        pTri0, nTri0,
        pNorm0);
    aInfoNearest[ip].sdf = sdf;
    aInfoNearest[ip].is_active = true;
    if (sdf + cc < 0) continue;
    aXYZt[ip * 3 + 0] += (sdf + cc) * n0.x;
    aXYZt[ip * 3 + 1] += (sdf + cc) * n0.y;
    aXYZt[ip * 3 + 2] += (sdf + cc) * n0.z;
  }
}

template<typename BV>
void Intersection_ImageRay_TriMesh3(
    std::vector<delfem2::CPtElm2<double> > &aPointElemSurf,
    unsigned int nheight,
    unsigned int nwidth,
    const float mMVPf[16],
    const std::vector<CNodeBVH2> &aNodeBVH,
    const std::vector<BV> &aAABB,
    const std::vector<double> &aXYZ, // 3d points
    const std::vector<unsigned int> &aTri) {
  aPointElemSurf.resize(nheight * nwidth);
  double mMVPd[16];
  for (int i = 0; i < 16; ++i) { mMVPd[i] = mMVPf[i]; }
  double mMVPd_inv[16];
  Inverse_Mat4(mMVPd_inv, mMVPd);
  std::vector<unsigned int> aIndElem;
  for (unsigned int ih = 0; ih < nheight; ++ih) {
    for (unsigned int iw = 0; iw < nwidth; ++iw) {
      //
      const double ps[4] = {-1. + (2. / nwidth) * (iw + 0.5), -1. + (2. / nheight) * (ih + 0.5), -1., 1.};
      const double pe[4] = {-1. + (2. / nwidth) * (iw + 0.5), -1. + (2. / nheight) * (ih + 0.5), +1., 1.};
      double qs[3];
      Vec3_Vec3Mat4_AffineProjection(qs, ps, mMVPd_inv);
      double qe[3];
      Vec3_Vec3Mat4_AffineProjection(qe, pe, mMVPd_inv);
      const CVec3d src1(qs);
      const CVec3d dir1 = CVec3d(qe) - src1;
      //
      aIndElem.resize(0);
      BVH_GetIndElem_Predicate(
          aIndElem,
          CIsBV_IntersectLine<BV, double>(src1.p, dir1.p),
          0, aNodeBVH, aAABB);
      if (aIndElem.empty()) { continue; } // no bv hit the ray
      std::map<double, CPtElm2<double>> mapDepthPES;
      IntersectionRay_MeshTri3DPart(
          mapDepthPES,
          src1, dir1,
          aTri, aXYZ, aIndElem, 1.0e-10);
      if (mapDepthPES.empty()) { continue; }
      aPointElemSurf[ih * nwidth + iw] = mapDepthPES.begin()->second;
    }
  }
}

template<typename BV>
void BuildBVH_MeshTri3D_Morton(
    std::vector<CNodeBVH2> &aNodeBVH,
    std::vector<BV> &aAABB,
    const std::vector<double> &aXYZ, // 3d points
    const std::vector<unsigned int> &aTri) {
  std::vector<double> aCent;
  double rad = CentsMaxRad_MeshTri3(
      aCent,
      aXYZ, aTri);
  double min_xyz[3], max_xyz[3];
  BoundingBox3_Points3(min_xyz, max_xyz,
                       aCent.data(), aCent.size() / 3);
  min_xyz[0] -= 1.0e-3 * rad;
  min_xyz[1] -= 1.0e-3 * rad;
  min_xyz[2] -= 1.0e-3 * rad;
  max_xyz[0] += 1.0e-3 * rad;
  max_xyz[1] += 1.0e-3 * rad;
  max_xyz[2] += 1.0e-3 * rad;
  std::vector<unsigned int> aSortedId;
  std::vector<std::uint32_t> aSortedMc;
  SortedMortenCode_Points3(aSortedId, aSortedMc,
                           aCent, min_xyz, max_xyz);
#ifndef NDEBUG
  Check_MortonCode_Sort(aSortedId, aSortedMc, aCent, min_xyz, max_xyz);
  Check_MortonCode_RangeSplit(aSortedMc);
#endif
  BVHTopology_Morton(aNodeBVH,
                     aSortedId, aSortedMc);
#ifndef NDEBUG
  Check_BVH(aNodeBVH, aCent.size() / 3);
#endif
  CLeafVolumeMaker_Mesh<BV, double> lvm(
      1.0e-10,
      aXYZ.data(), aXYZ.size() / 3,
      aTri.data(), aTri.size() / 3, 3);
  BVH_BuildBVHGeometry(aAABB,
                       0, aNodeBVH,
                       lvm);
}

} // namespace delfem2

#endif
