/*
 * Copyright (c) 2019 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */

/**
 * @brief dynamic triangle mesh class & functions
 */

#ifndef DFM2_DTRI2_V2DTRI_H
#define DFM2_DTRI2_V2DTRI_H

#include <map>
#include <algorithm>
#include <stack>

#include "delfem2/dfm2_inline.h"
#include "delfem2/vec2.h"
#include "delfem2/dtri.h"
#include "delfem2/geo_polygon2.h"

// -------------------------------

namespace delfem2 {

/**
 * @brief assert whether the dynamic triangulation data is broken or not
 * @brief this function do nothing if the code is compiled with NDEBUG
 */
DFM2_INLINE void CheckTri(
    const std::vector<CDynPntSur> &aPo3D,
    const std::vector<CDynTri> &aSTri,
    const std::vector<CVec2d> &aXYZ);

/**
 * @brief compute delaunay trangulation around a point
 * @param ipo0 point index
 */
DFM2_INLINE void DelaunayAroundPoint(
    unsigned int ipo0,
    std::vector<CDynPntSur> &aPo,
    std::vector<CDynTri> &aTri,
    const std::vector<CVec2d> &aVec2);

DFM2_INLINE void MeshTri2D_Export(
    std::vector<double> &aXY_out,
    std::vector<unsigned int> &aTri_out,
    const std::vector<CVec2d> &aVec2,
    const std::vector<CDynTri> &aTri_in);

DFM2_INLINE void Meshing_Initialize(
    std::vector<CDynPntSur> &aPo2D,
    std::vector<CDynTri> &aTri,
    std::vector<CVec2d> &aVec2);

DFM2_INLINE void FlagConnected(
    std::vector<int> &inout_flg,
    const std::vector<CDynTri> &aTri_in,
    unsigned int itri0_ker,
    int iflag);


/**
 * @brief delete triangle where flag is raised
 * @details aFlg.size() == aTri_in.size() needs to hold.
 */
DFM2_INLINE void DeleteTriFlag(
    std::vector<CDynTri> &aTri_in,
    std::vector<int> &aFlg,
    int flag);

DFM2_INLINE void EnforceEdge(
    std::vector<CDynPntSur> &aPo2D,
    std::vector<CDynTri> &aTri,
    unsigned int ip0,
    unsigned int ip1,
    const std::vector<CVec2d> &aVec2);

DFM2_INLINE void Meshing_SingleConnectedShape2D(
    std::vector<CDynPntSur> &aPo2D,
    std::vector<CVec2d> &aVec2,
    std::vector<CDynTri> &aETri,
    const std::vector<int> &loopIP_ind,
    const std::vector<int> &loopIP);

DFM2_INLINE void DeleteUnrefPoints(
    std::vector<CVec2d> &aVec2,
    std::vector<CDynPntSur> &aPo2D,
    std::vector<CDynTri> &aTri_in,
    const std::vector<unsigned int> &aPoDel);

DFM2_INLINE void DeletePointsFlag(
    std::vector<CVec2d> &aVec1,
    std::vector<CDynPntSur> &aPo1,
    std::vector<CDynTri> &aTri,
    std::vector<int> &aFlgPnt1,
    int iflg);

DFM2_INLINE void MakeSuperTriangle(
    std::vector<CVec2d> &aVec2,
    std::vector<CDynPntSur> &aPo2D,
    std::vector<CDynTri> &aTri,
    const double bound_2d[4]);

DFM2_INLINE void AddPointsMesh(
    const std::vector<CVec2d> &aVec2,
    std::vector<CDynPntSur> &aPo2D,
    std::vector<CDynTri> &aTri,
    int ipoin,
    double MIN_TRI_AREA);

DFM2_INLINE void MakeInvMassLumped_Tri(
    std::vector<double> &aInvMassLumped,
    double rho,
    const std::vector<CVec2d> &aVec2,
    const std::vector<CDynTri> &aETri);

DFM2_INLINE void MinMaxTriArea(
    double &min_area,
    double &max_area,
    const std::vector<CVec2d> &aVec2,
    const std::vector<CDynTri> &aETri);

/**
 * Convert dynamic triangle mesh data into array of vertex coordinate and their connectivities
 * @param[out] aXY
 * @param[out] aTri
 * @param[in] aVec2
 * @param[in] aETri
 */
DFM2_INLINE void CMeshTri2D(
    std::vector<double> &aXY,
    std::vector<unsigned int> &aTri,
    const std::vector<CVec2d> &aVec2,
    const std::vector<CDynTri> &aETri);

/**
 * Generate mesh inside a region specified by polygons
 * First polygon specify the outside of the region, subsequent polygons specifies the holes
 * @param[out] aPo2D
 * @param[out] aETri
 * @param[out] aVec2
 * @param[in] aaXY input array of polygons
 * @param[in] resolution_edge
 * @param[in] resolution_face
 */
DFM2_INLINE void GenMesh(
    std::vector<CDynPntSur> &aPo2D,
    std::vector<CDynTri> &aETri,
    std::vector<CVec2d> &aVec2,
    const std::vector<std::vector<double> > &aaXY,
    double resolution_edge,
    double resolution_face);

class CInputTriangulation {
 public:
  virtual double edgeLengthRatio(double px, double py) const = 0;
};

class CInputTriangulation_Uniform : public CInputTriangulation {
 public:
  explicit CInputTriangulation_Uniform(double elen) : elen(elen) {}
  double edgeLengthRatio(
      [[maybe_unused]] double px,
      [[maybe_unused]] double py) const override {
    return 1.0;
  }
 public:
  double elen;
};

/**
 * @param aFlagTri a map from triangle index to cad face indeix
 */
void MeshingInside(
    std::vector<CDynPntSur> &aPo2D,
    std::vector<CDynTri> &aTri,
    std::vector<CVec2d> &aVec2,
    std::vector<int> &aFlagPnt,
    std::vector<unsigned int> &aFlagTri,
    //
    size_t nPointFix,
    unsigned int nflgpnt_offset,
    double len,
    const CInputTriangulation &mesh_density);

class CCmdRefineMesh {
 public:
  class CCmdEdge {
   public:
    CCmdEdge(int i0, int i1, double s0) {
      if (i0 < i1) {
        ipo0 = i0;
        ipo1 = i1;
        r0 = s0;
      }
      else {
        ipo0 = i1;
        ipo1 = i0;
        r0 = 1 - s0;
      }
    }
    bool operator<(const CCmdEdge &rhs) const {
      if (ipo0 != rhs.ipo0) { return ipo0 < rhs.ipo0; }
      return ipo1 < rhs.ipo1;
    }
   public:
    int ipo_new{}, ipo0, ipo1;
    double r0;
  };
 public:
  void Interpolate(double *pVal, int np, int ndim) const {
    for (const auto &cmd : aCmdEdge) {
      const int i0 = cmd.ipo0;
      assert(i0 >= 0 && i0 < np);
      const int i1 = cmd.ipo1;
      assert(i1 >= 0 && i1 < np);
      const int i2 = cmd.ipo_new;
      if (i2 >= np || i2 < 0) { continue; }
      double r0 = cmd.r0;
      for (int idim = 0; idim < ndim; ++idim) {
        pVal[i2 * ndim + idim] = r0 * pVal[i0 * ndim + idim] + (1 - r0) * pVal[i1 * ndim + idim];
      }
    }
  }
 public:
  std::vector<CCmdEdge> aCmdEdge;
};

void RefinementPlan_EdgeLongerThan_InsideCircle(
    CCmdRefineMesh &aCmd,
    double elen,
    double px, double py, double rad,
    const std::vector<CDynPntSur> &aPo2D,
    const std::vector<CVec2d> &aVec2,
    const std::vector<CDynTri> &aETri);

void RefineMesh(
    std::vector<CDynPntSur> &aPo3D,
    std::vector<CDynTri> &aSTri,
    std::vector<CVec2d> &aVec2,
    CCmdRefineMesh &aCmd);

class CMeshDynTri2D {
 public:
  void Initialize(const double *aXY, int nPo,
                  const unsigned int *aTri, int nTri) {
    aVec2.resize(nPo);
    for (int ipo = 0; ipo < nPo; ipo++) {
      aVec2[ipo].p[0] = aXY[ipo * 2 + 0];
      aVec2[ipo].p[1] = aXY[ipo * 2 + 1];
    }
    InitializeMesh(aEPo, aETri,
                   aTri, nTri, nPo);
  }
  void setXY(const double *aXY, int nPo) {
    assert((int) aVec2.size() == nPo);
    for (int ipo = 0; ipo < nPo; ipo++) {
      aVec2[ipo].p[0] = aXY[ipo * 2 + 0];
      aVec2[ipo].p[1] = aXY[ipo * 2 + 1];
    }
  }
  void Check() const {
    AssertDTri(aETri);
    AssertMeshDTri(aEPo, aETri);
    CheckTri(aEPo, aETri, aVec2);
  }
  std::vector<double> MinMax_XYZ() const {
    double x_min, x_max, y_min, y_max;
    x_min = x_max = aVec2[0].x;
    y_min = y_max = aVec2[0].y;
    for (unsigned int ipo = 0; ipo < aEPo.size(); ipo++) {
      const double x = aVec2[ipo].x;
      const double y = aVec2[ipo].y;
      x_min = (x_min < x) ? x_min : x;
      x_max = (x_max > x) ? x_max : x;
      y_min = (y_min < y) ? y_min : y;
      y_max = (y_max > y) ? y_max : y;
    }
    std::vector<double> bb;
    bb.push_back(x_min);
    bb.push_back(x_max);
    bb.push_back(y_min);
    bb.push_back(y_max);
    bb.push_back(0.0);
    bb.push_back(0.0);
    return bb;
  }
  int insertPointElem(int itri0, double r0, double r1) {
    CVec2d v2;
    {
      unsigned int i0 = aETri[itri0].v[0];
      unsigned int i1 = aETri[itri0].v[1];
      unsigned int i2 = aETri[itri0].v[2];
      v2 = r0 * aVec2[i0] + r1 * aVec2[i1] + (1 - r0 - r1) * aVec2[i2];
    }
    const auto ipo0 = static_cast<unsigned int>(aEPo.size());
    aVec2.push_back(v2);
    aEPo.emplace_back();
    InsertPoint_Elem(ipo0, itri0, aEPo, aETri);
    return ipo0;
  }
  void DelaunayAroundPoint(int ipo) {
    delfem2::DelaunayAroundPoint(ipo, aEPo, aETri, aVec2);
  }
  void meshing_loops(const std::vector<std::vector<double> > &aaXY,
                     double edge_length) {
    std::vector<int> loopIP_ind, loopIP;
    {
      JArray_FromVecVec_XY(loopIP_ind, loopIP, aVec2,
                           aaXY);
      if (!CheckInputBoundaryForTriangulation(loopIP_ind, aVec2)) {
        return;
      }
      FixLoopOrientation(loopIP,
                         loopIP_ind, aVec2);
      if (edge_length > 10e-10) {
        ResamplingLoop(
            loopIP_ind, loopIP, aVec2,
            edge_length);
      }
    }
    ////
    Meshing_SingleConnectedShape2D(aEPo, aVec2, aETri,
                                   loopIP_ind, loopIP);
    if (edge_length > 1.0e-10) {
      CInputTriangulation_Uniform param(1.0);
      std::vector<int> aFlgPnt(aVec2.size(), 0);
      std::vector<unsigned int> aFlgTri(aETri.size(), 0);
      MeshingInside(aEPo, aETri, aVec2, aFlgPnt, aFlgTri,
                    aVec2.size(), 0, edge_length, param);
    }
  }
  void RefinementPlan_EdgeLongerThan_InsideCircle(CCmdRefineMesh &aCmd,
                                                  double elen,
                                                  double px, double py, double rad) {
    delfem2::RefinementPlan_EdgeLongerThan_InsideCircle(aCmd,
                                                        elen, px, py, rad,
                                                        aEPo, aVec2, aETri);
    RefineMesh(aEPo, aETri, aVec2, aCmd);
    assert(aEPo.size() == aVec2.size());
  }
  void Export_StlVectors(std::vector<double> &aXY, std::vector<unsigned int> &aTri) const {
    MeshTri2D_Export(aXY, aTri, aVec2, aETri);
  }
  void Clear() {
    aEPo.clear();
    aETri.clear();
    aVec2.clear();
  }
  size_t nTri() const { return aETri.size(); }
  size_t nPoint() const { return aEPo.size(); }
  void DeleteTriEdge(int itri, int iedge) { CollapseEdge_MeshDTri(itri, iedge, aEPo, aETri); }
 public:
  std::vector<CDynPntSur> aEPo;
  std::vector<CDynTri> aETri;
  std::vector<CVec2d> aVec2;
};

} // namespace delfem2


#ifndef DFM2_STATIC_LIBRARY
#  include "delfem2/dtri2_v2dtri.cpp"
#endif

#endif
