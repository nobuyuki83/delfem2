/*
 * Copyright (c) 2019 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */

#ifndef DFM2_DTRI_H
#define DFM2_DTRI_H

#include <vector>
#include <climits>

#include "delfem2/dfm2_inline.h"

// ------------------------------------------------------------------

namespace delfem2 {

/**
 * @brief dynamic triangle class
 */
class CDynTri {
 public:
  /**
   * @brief index of vertex
   */
  unsigned int v[3];

  /**
   * @brief index of face that is adjacent to ith edge;
   * @detail The i-th edge is the edge with (i+1)%3-th and (i+2)%3-th vertex.
   * if this edge is on the boundary, the value is set to UINT_MAX
   */
  unsigned int s2[3];
};

class CDynPntSur {
 public:
  CDynPntSur() {
    e = UINT_MAX;
    d = 0;
  }
  CDynPntSur(const CDynPntSur &) = default;
  CDynPntSur(int ielem, unsigned int idir)
  : e(ielem), d(idir) {}
 public:
  /**
   * @brief index of elements
   * @detail this point is isolated if it is UINT_MAX
   */
  unsigned int e;
  unsigned int d;
};

DFM2_INLINE unsigned int FindAdjEdgeIndex(
    const CDynTri &t0,
    unsigned int ied0,
    const std::vector<delfem2::CDynTri> &aTri);

DFM2_INLINE bool JArray_MakeElSuP(
    std::vector<int> &elsup_ind,
    std::vector<int> &elsup,
    const std::vector<CDynTri> &aTri,
    const unsigned int npoin);

DFM2_INLINE void JArray_PSuP(
    std::vector<int> &psup_ind,
    std::vector<int> &psup,
    const std::vector<CDynTri> &aTri,
    const unsigned int npoin,
    const std::vector<int> &elsup_ind,
    const std::vector<int> &elsup);

DFM2_INLINE void MakeInnerRelationTri(
    std::vector<CDynTri> &aTri,
    const unsigned int npoin,
    const std::vector<int> &elsup_ind,
    const std::vector<int> &elsup);

DFM2_INLINE void JArray_PSuP(
    unsigned int *const edge_ind,
    unsigned int &nedge,
    unsigned int *&edge,
    const std::vector<CDynTri> &aTri,
    const unsigned int npoin,
    const std::vector<int> &elsup_ind,
    const std::vector<int> &elsup);

DFM2_INLINE void AssertDTri(
    const std::vector<CDynTri> &aETri);

DFM2_INLINE void AssertMeshDTri(
    const std::vector<CDynPntSur> &aEPo2,
    const std::vector<CDynTri> &aETri);

DFM2_INLINE void InitializeMesh(
    std::vector<CDynPntSur> &aEPo2,
    std::vector<CDynTri> &aETri,
    const unsigned int *aTri,
    const size_t nTri,
    const size_t nXYZ);

DFM2_INLINE bool FindEdge_LookAroundPoint(
    unsigned int &itri0,
    unsigned int &inotri0,
    unsigned int &inotri1,
    //
    const unsigned int ipo0,
    const unsigned int ipo1,
    const std::vector<CDynPntSur> &aPo,
    const std::vector<CDynTri> &aTri);

DFM2_INLINE bool FindPointAroundPoint(
    std::vector<unsigned int> &aIP,
    const unsigned int ipo0,
    const std::vector<CDynPntSur> &aPo,
    const std::vector<CDynTri> &aTri);

/**
 * @details itir0 stays the same if the edge of the traingle is not found
 */
DFM2_INLINE void FindEdge_LookAllTriangles(
    unsigned int &itri0,
    unsigned int &iedtri0,
    const unsigned int ipo0,
    const unsigned int ipo1,
    const std::vector<CDynTri> &aETri);

DFM2_INLINE void GetTriArrayAroundPoint(
    std::vector<std::pair<unsigned int, unsigned int> > &aTriSurPo,
    unsigned int ipoin,
    const std::vector<CDynPntSur> &aEPo2,
    const std::vector<CDynTri> &aETri);

DFM2_INLINE bool MoveCCW(
    unsigned int &itri_cur,
    unsigned int &inotri_cur,
    unsigned int itri_adj,
    const std::vector<CDynTri> &aTri);

DFM2_INLINE bool MoveCW(
    unsigned int &itri_cur,
    unsigned int &inotri_cur,
    unsigned int itri_adj,
    const std::vector<CDynTri> &aTri);

// ---------------
// topology edit

/**
 * @brief flip edge of a triangle mesh
 * @details after the flip, itri0 will have the point that were previously aTri[itri0].v[(ied0+1)%3]
 * aDPo[aDTri[itri0].v[ied0]].e stays itri0
 * @param itri0 trinangle index 0<=itri0<aTri.size()
 * @param ied0 edge index. 0<=ied0<3
 * @return return false if the edge is on the boundary
 */
DFM2_INLINE bool FlipEdge(
    unsigned int itri0,
    unsigned int ied0,
    std::vector<CDynPntSur> &aDPo,
    std::vector<CDynTri> &aDTri);


// ----------------------
// insert point

DFM2_INLINE bool InsertPoint_ElemEdge(
    const unsigned int ipo_ins,  //!< the index of the new point
    const unsigned int itri_ins, //!< triangle index
    const unsigned int ied_ins,  //!< edge index
    std::vector<CDynPntSur> &aEPo2,
    std::vector<CDynTri> &aETri);

DFM2_INLINE bool InsertPoint_Elem(
    const unsigned int ipo_ins,
    const unsigned int itri_ins,
    std::vector<CDynPntSur> &aEPo2,
    std::vector<CDynTri> &aETri);

// -----------------
// delete point

DFM2_INLINE bool DeleteTri(
    unsigned int itri_to,
    std::vector<CDynPntSur> &aEPo2,
    std::vector<CDynTri> &aETri);

/**
 * @brief delete a point aETri[itri_del].v[(ied_del+2)%3]
 */
DFM2_INLINE bool CollapseEdge_MeshDTri(
    const unsigned int itri_del,
    const unsigned int ied_del,
    std::vector<CDynPntSur> &aEPo2,
    std::vector<CDynTri> &aETri);


// ------------------------------------------

DFM2_INLINE void extractHoles(
    std::vector<std::vector<int> > &aIndP_Hole,
    const int npo,
    const std::vector<CDynTri> &aETri);

}

#ifndef DFM2_STATIC_LIBRARY
#  include "delfem2/dtri.cpp"
#endif

#endif // DFM2_DTRI_H
