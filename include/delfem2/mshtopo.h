/*
 * Copyright (c) 2019 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */

/**
 * @file functions to analyze mesh topology for static meshes
 * @details the functions only care about the topology. Geometry (coordinate) information is not handled in this file
 */


// TODO: separate mixed elem

#ifndef DFM2_MSHTOPO_H
#define DFM2_MSHTOPO_H

#include "delfem2/dfm2_inline.h"
#include <cstdio>
#include <vector>


namespace delfem2 {

enum MESHELEM_TYPE
{
  MESHELEM_TRI = 0,
  MESHELEM_TET = 1,
  MESHELEM_QUAD = 2,
  MESHELEM_HEX = 3,
  MESHELEM_PYRAMID = 4,
  MESHELEM_WEDGE = 5,
  MESHELEM_VOX = 6,
  MESHELEM_LINE = 7,
};

// 5 : VTK_TRIANGLE
// 9 : VTK_QUAD
// 10: VTK_TETRA
// 12: VTK_HEXAHEDRON
// 13: VTK_WEDGE
// 14: VTK_PYRAMD

////
const int mapMeshElemType2NNodeElem[8] = {
  3, // TRI 0
  4, // TET 1
  4, // QUAD 2
  8, // HEX 3
  5, // PYRAM 4
  6, // WEDGE 5
  8,  // VOX 6
  2, // LINE 7
};
////
const int mapMeshElemType2NFaceElem[7] = {
  3, // TRI
  4, // TET
  4, // QUAD
  6, // HEX
  5, // Pyramid
  5, // Wedge
  6  // VOX
};
////
const int mapMeshElemTypeIndFace2NNoelElemFace[7][6] =
{
  {2,2,2,-1,-1,-1}, // TRI
  {3,3,3, 3,-1,-1}, // TET
  {2,2,2, 2,-1,-1}, // QUAD
  {4,4,4, 4, 4, 4}, // HEX
  {4,3,3, 3, 3,-1}, // Pyramid
  {4,4,4, 3, 3,-1}, // Wedge
  {4,4,4, 4, 4, 4}, // VOX
};
////
const int noelElemFace_Tri[3][4] = {
  { 1, 2,-1,-1 }, //edge 0
  { 2, 0,-1,-1 }, //edge 1
  { 0, 1,-1,-1 }, //edge 2
};
const int noelElemFace_Tet[4][4] ={
  { 1, 2, 3,-1 },
  { 0, 3, 2,-1 },
  { 0, 1, 3,-1 },
  { 0, 2, 1,-1 }
};
const int noelElemFace_Quad[4][4] = {
  { 0, 1,-1,-1 },
  { 1, 2,-1,-1 },
  { 2, 3,-1,-1 },
  { 3, 0,-1,-1 },
};
const int noelElemFace_Vox[6][4] = { // this numbering is corresponds to VTK_VOX
  { 0, 4, 6, 2 }, // -x
  { 1, 3, 7, 5 }, // +x
  { 0, 1, 5, 4 }, // -y
  { 2, 6, 7, 3 }, // +y
  { 0, 2, 3, 1 }, // -z
  { 4, 5, 7, 6 }  // +z
};
const int noelElemFace_Hex[6][4] = { // this numbering is corresponds to VTK_HEX
  { 0, 4, 7, 3 }, // -x
  { 1, 2, 6, 5 }, // +x
  { 0, 1, 5, 4 }, // -y
  { 3, 7, 6, 2 }, // +y
  { 0, 3, 2, 1 }, // -z
  { 4, 5, 6, 7 }  // +z
};
const int noelElemFace_Pyram[5][4] = {
  { 0, 3, 2, 1},
  { 0, 1, 4,-1},
  { 1, 2, 4,-1},
  { 2, 3, 4,-1},
  { 3, 0, 4,-1}
};
const int noelElemFace_Wed[5][4] = {
  { 0, 2, 5, 3},
  { 1, 0, 3, 4},
  { 2, 1, 4, 5},
  { 0, 1, 2,-1},
  { 3, 5, 4,-1}
};

const int mapMeshElemType2NEdgeElem[7] = {
  3, // TRI
  6, // TET
  4, // QUAD
  12, // HEX
  8, // PYRAM
  9, // WEDGE
  12  // VOX
};
const int noelElemEdge_Tri[3][2] = {
  { 0, 1},
  { 1, 2},
  { 2, 0},
};
const int noelElemEdge_Tet[6][2] = {
  { 0, 1},
  { 0, 2},
  { 0, 3},
  { 1, 2},
  { 1, 3},
  { 2, 2},
};
const int noelElemEdge_Quad[4][2] = {
  { 0, 1},
  { 1, 2},
  { 2, 3},
  { 3, 0},
};
const int noelElemEdge_Hex[12][2] = {
  {0,1},{1,2},{2,3},{3,0},
  {4,5},{5,6},{6,7},{7,4},
  {0,4},{1,5},{2,6},{3,7}
};
// TODO: this looks wrong...
const int noelElemEdge_Vox[12][2] = {
  {0,1},{3,2},{4,5},{7,6},
  {0,3},{1,2},{4,7},{5,6},
  {0,4},{1,5},{3,7},{2,6} 
};

// -----------------------------
inline int nNodeElem(MESHELEM_TYPE type){
  return mapMeshElemType2NNodeElem[type];
}
inline int nFaceElem(MESHELEM_TYPE type){
  return mapMeshElemType2NFaceElem[type];
}
inline int nNodeElemFace(MESHELEM_TYPE type,int iface){
  return mapMeshElemTypeIndFace2NNoelElemFace[type][iface];
}
inline const int (*noelElemFace(MESHELEM_TYPE type))[4]
{
  const int (*noelElemFace0[7])[4] = {
    noelElemFace_Tri,
    noelElemFace_Tet,
    noelElemFace_Quad,
    noelElemFace_Hex,
    noelElemFace_Pyram,
    noelElemFace_Wed,
    noelElemFace_Vox,
  };
  return noelElemFace0[type];
}
inline const int (*noelElemEdge(MESHELEM_TYPE type))[2]
{
  const int (*noelElemEdge0[4])[2] = {
    noelElemEdge_Tri,
    noelElemEdge_Tet,
    noelElemEdge_Quad,
    noelElemEdge_Hex
  };
  return noelElemEdge0[type];
}
  
// ---------------------------------------------
// function related to jagged array

DFM2_INLINE void JArray_Sort(
    const std::vector<unsigned int>& index,
    std::vector<unsigned int>& array);
DFM2_INLINE void JArray_Sort(
    const unsigned int* index,
    unsigned int size,
    unsigned int* array);

DFM2_INLINE void JArray_AddDiagonal(
    std::vector<unsigned int> &psup_ind1,
    std::vector<unsigned int> &psup1,
    const unsigned int *psup_ind0, int npsup_ind0,
    const unsigned int *psup0, int npsup0);

DFM2_INLINE void JArray_Print(
    const std::vector<int>& index,
    const std::vector<int>& array);

DFM2_INLINE void JArray_AddMasterSlavePattern(
    std::vector<unsigned int> &index,
    std::vector<unsigned int> &array,
    const unsigned int* aMSFlag,
    int ndim,
    const unsigned int *psup_ind0,
    int npsup_ind0,
    const unsigned int *psup0);

/**
 * @details compute 2-ring neighborhood from 1-ring neighborhood
 */
DFM2_INLINE void JArray_Extend(
    std::vector<unsigned int>& psup_ind1,
    std::vector<unsigned int>& psup1,
    const unsigned int *psup_ind0,
    unsigned int npsup_ind0,
    const unsigned int *psup0);

// ---------------------------------------------------

DFM2_INLINE unsigned FindAdjEdgeIndex(
    unsigned int itri,
    unsigned int ied,
    unsigned int jtri,
    const unsigned int* aTri);

DFM2_INLINE void convert2Tri_Quad(
    std::vector<unsigned int>& aTri,
    const std::vector<unsigned int>& aQuad);

DFM2_INLINE void Convert2Tri_MeshMix(
    std::vector<unsigned int>& aTri,
    //
    const std::vector<unsigned int>& aElemInd,
    const std::vector<unsigned int>& aElem,
    const std::vector<delfem2::MESHELEM_TYPE>& aElemType);

/**
 * @brief Make quad mesh from tri mesh by merging adjacent triangle elements
 * @param aQuad (out) element index of quad mesh
 * @param aTri (in) element index of tri mesh
 * @param nTri (in) number of triangle mesh
 * @param np (in) number of points
 */
DFM2_INLINE void ElemQuad_DihedralTri(
    std::vector<unsigned int>& aQuad,
    const unsigned int* aTri,
    unsigned int nTri,
    unsigned int np);

DFM2_INLINE void FlipElement_Tri(
    std::vector<unsigned int>& aTri);

DFM2_INLINE void FlipElement_MeshMix(
    std::vector<int>& aElem_Flip,
    //
    const std::vector<unsigned int>& aElemInd,
    const std::vector<unsigned int>& aElem,
    const std::vector<delfem2::MESHELEM_TYPE>& aElemType);


/**
 * @brief make elem surrounding point
 */
DFM2_INLINE void JArray_ElSuP_MeshElem(
    std::vector<unsigned int> &elsup_ind,
    std::vector<unsigned int> &elsup,
    //
    const unsigned int *pElem,
    unsigned int nElem,
    unsigned int nPoEl,
    unsigned int nPo);

/**
 * @brief make elem surrounding point for triangle mesh
 */
DFM2_INLINE void JArray_ElSuP_MeshTri(
    std::vector<unsigned int> &elsup_ind,
    std::vector<unsigned int> &elsup,
    //
    const std::vector<unsigned int> &aTri,
    int nXYZ);

/**
 * @brief elem surrounding point for mixed element
 */
DFM2_INLINE void JArray_ElSuP_MeshMix(
    std::vector<unsigned int> &elsup_ind,
    std::vector<unsigned int> &elsup,
    //
    const std::vector<unsigned int> &aElemInd,
    const std::vector<unsigned int> &aElem,
    int nPo);


// -----------------
// elem sur elem

/**
 * @brief compute adjacent element index for mesh element
 * @param aElSuEl (out) neighbouring element index (UINT_MAX for boundary)
 * @param aEl (in) array of connectivity
 * @param nEl (in) number of elements
 * @param nNoEl (in) number of nodes in a element
 * @param elsup_ind (in) jagged array index of "elem surrounding point"
 * @param elsup (in) jagged array value of "elem surrounding point"
 * @param nfael (in) number of neibouring elements
 * @param nnofa (in) how many nodes are shared with a nighbouring element
 * @param noelElemFace
 */
DFM2_INLINE void ElSuEl_MeshElem(
    std::vector<unsigned int> &aElSuEl,
    const unsigned int *aEl,
    unsigned int nEl,
    int nNoEl,
    const std::vector<unsigned int> &elsup_ind,
    const std::vector<unsigned int> &elsup,
    int nfael,
    int nnofa,
    const int (*noelElemFace)[4]);

/**
 * @brief compute adjacent element index for mesh element
 * @param aElSuEl (ou) adjacent element index for element edge/face (UINT_MAX if face/edge is on the boundary)
 * @param aElem (in) elemnet index
 * @param nElem (in) number of elements
 * @param type (in) type of element
 * @param nXYZ (in) number of points
 */
DFM2_INLINE void ElSuEl_MeshElem(
    std::vector<unsigned int> &aElSuEl,
    const unsigned int *aElem, size_t nElem,
    delfem2::MESHELEM_TYPE type,
    size_t nXYZ);

DFM2_INLINE void ElSuEl_MeshMix(
    std::vector<int> &aElemFaceInd,
    std::vector<int> &aElemFaceRel,
    const std::vector<unsigned int> &aElemInd,
    const std::vector<unsigned int> &aElem,
    const std::vector<delfem2::MESHELEM_TYPE> &aElemType,
    int nPo);

/**
 * @brief the relationship between neighboring elements for mixed mesh
 */
DFM2_INLINE void ElSuEl_MeshMix(
    std::vector<int> &aElemFaceInd,
    std::vector<int> &aElemFaceRel,
    //
    const std::vector<unsigned int> &aElemInd,
    const std::vector<unsigned int> &aElem,
    const std::vector<delfem2::MESHELEM_TYPE> &aElemType,
    const std::vector<unsigned int> &elsup_ind,
    const std::vector<unsigned int> &elsup);

// -------------------
// make boundary

DFM2_INLINE void Boundary_MeshMix(
    std::vector<unsigned int>& aElemInd_Bound,
    std::vector<unsigned int>& aElem_Bound,
    std::vector<delfem2::MESHELEM_TYPE>& aElemType_Bound,
    //
    const std::vector<unsigned int>& aElemInd,
    const std::vector<unsigned int>& aElem,
    const std::vector<delfem2::MESHELEM_TYPE>& aElemType,
    const std::vector<int>& aElemFaceInd,
    const std::vector<int>& aElemFaceRel);

/**
 * @brief make point surrounding point
 * @details psup -> edge bidirectional
 * edge unidir (ip0<ip1)
 * line (array of 2)
 */
DFM2_INLINE void JArrayPointSurPoint_MeshOneRingNeighborhood(
    std::vector<unsigned int>& psup_ind,
    std::vector<unsigned int>& psup,
    //
    const unsigned int* pElem,
    const std::vector<unsigned int> &elsup_ind,
    const std::vector<unsigned int> &elsup,
    unsigned int nnoel,
    unsigned int nPoint);

/**
 * @brief compute indexes of points surrounding a point as a jagged array
 * @param nPoEl number of nodes in an element 
 */
DFM2_INLINE void JArray_PSuP_MeshElem(
    std::vector<unsigned int>& psup_ind,
    std::vector<unsigned int>& psup,
    //
    const unsigned int* pElem,
    size_t nEl,
    unsigned int nPoEl,
    size_t nPo);

DFM2_INLINE void makeOneRingNeighborhood_TriFan(
    std::vector<int>& psup_ind,
    std::vector<int>& psup,
    //
    const std::vector<int>& aTri,
    const std::vector<int>& aTriSurRel,
    const std::vector<int>& elsup_ind,
    const std::vector<int>& elsup,
    int np);

DFM2_INLINE void JArrayEdgeUnidir_PointSurPoint(
    std::vector<unsigned int> &edge_ind,
    std::vector<unsigned int> &edge,
    //
    const std::vector<unsigned int> &psup_ind,
    const std::vector<unsigned int> &psup);

DFM2_INLINE void JArrayEdge_MeshElem(
    std::vector<unsigned int> &edge_ind,
    std::vector<unsigned int> &edge,
    //
    const unsigned int* aElm0,
    delfem2::MESHELEM_TYPE elem_type,
    const std::vector<unsigned int> &elsup_ind,
    const std::vector<unsigned int> &elsup,
    bool is_bidirectional);
  
DFM2_INLINE void MeshLine_JArrayEdge(
    std::vector<unsigned int>& aLine,
    //
    const std::vector<unsigned int> &psup_ind,
    const std::vector<unsigned int> &psup);

DFM2_INLINE void MeshLine_MeshElem(
    std::vector<unsigned int>& aLine,
    //
    const unsigned int* aElm0,
    unsigned int nElem,
    delfem2::MESHELEM_TYPE elem_type,
    unsigned int nPo);

// ------------------------------------

DFM2_INLINE void MarkConnectedElements(
    std::vector<int>& aIndGroup,
    int itri_ker,
    int igroup,
    const std::vector<int>& aTriSurRel,
    int nfael);

DFM2_INLINE void MarkConnectedElements(
    std::vector<int>& aIndGroup,
    unsigned int itri_ker,
    int igroup,
    const std::vector<int>& aElemFaceInd,
    const std::vector<int>& aElemFaceRel);

DFM2_INLINE void MakeGroupElem(
    int& ngroup,
    std::vector<int>& aIndGroup,
    const std::vector<int>& aElem,
    const std::vector<int>& aElemSurRel,
    int nfael,
    int nnoel);

DFM2_INLINE void MakeGroupElem_Tri(
    int& ngroup,
    std::vector<int>& aIndGroup,
    const std::vector<int>& aTri,
    const std::vector<int>& aTriSurRel);

DFM2_INLINE void MakeGroupElem(
    int& ngroup,
    std::vector<int>& aIndGroup,
    const std::vector<unsigned int>& aElemInd,
    const std::vector<unsigned int>& aElem,
    const std::vector<int>& aElemFaceInd,
    const std::vector<int>& aElemFaceRel);

DFM2_INLINE void MakeGroupElem_MeshMix(
    int& ngroup,
    std::vector<int>& aIndGroup,
    //
    const std::vector<unsigned int>& aElemInd,
    const std::vector<unsigned int>& aElem,
    const std::vector<delfem2::MESHELEM_TYPE>& aElemType,
    int nPo);

// ------------------------------

DFM2_INLINE void ClipGroup_MeshMix(
    std::vector<unsigned int>& aElemInd1,
    std::vector<unsigned int>& aElem1,
    std::vector<delfem2::MESHELEM_TYPE>& aElemType1,
    //
    const std::vector<unsigned int>& aElemInd,
    const std::vector<unsigned int>& aElem,
    const std::vector<delfem2::MESHELEM_TYPE>& aElemType,
    int igroup,
    const std::vector<int>& aIndGroup);

// -----------------------------------------------

/**
 * @brief making topology for subdivision of quad
 * @details new points is in the order of [old points], [edge points], [face points]
 * @param aQuad1 (out) new connectivity
 * @param aEdgeFace0 (out) two end points on a edge and two quads touching the edge
 */
DFM2_INLINE void SubdivTopo_MeshQuad(
    std::vector<unsigned int> &aQuad1,
    std::vector<unsigned int> &psup_ind,
    std::vector<unsigned int> &psup,
    std::vector<int> &aEdgeFace0,
    const unsigned int *aQuad0, unsigned int nQuad0,
    unsigned int nPoint0);

DFM2_INLINE void SubdivTopo_MeshHex(
    std::vector<unsigned int> &aHex1,
    std::vector<unsigned int> &psupIndHex0,
    std::vector<unsigned int> &psupHex0,
    std::vector<unsigned int> &aQuadHex0,
    const unsigned int *aHex0, int nHex0,
    int nhp0);

DFM2_INLINE void SubdivTopo_MeshTet(
    std::vector<unsigned int> &aTet1,
    std::vector<unsigned int> &psup_ind,
    std::vector<unsigned int> &psup,
    const unsigned int *aTet0, int nTet0,
    unsigned int nPoint0);

DFM2_INLINE int findEdge(
    unsigned int ip0, unsigned int ip1,
    const std::vector<unsigned int> &psup_ind,
    const std::vector<unsigned int> &psup);

DFM2_INLINE int findFace(
    unsigned int ip0, unsigned int ip1, unsigned int ip2, unsigned int ip3,
    const std::vector<unsigned int> &aQuad,
    const std::vector<unsigned int> &elsupInd,
    const std::vector<unsigned int> &elsup);

// ---------------------------------------------------

DFM2_INLINE void AddElement(
    const delfem2::MESHELEM_TYPE& femelem_type,
    const std::vector<int>& aElemIn,
    //
    std::vector<unsigned int>& aElemInd,
    std::vector<unsigned int>& aElem,
    std::vector<delfem2::MESHELEM_TYPE>& aElemType);

class CElemMixed{
public:
  CElemMixed(){
    aElemInd.resize(1,0);
  }
  void AddElement(
      const MESHELEM_TYPE& femelem_type,
      const std::vector<int>& aElemIn) {
    ::delfem2::AddElement(femelem_type,aElemIn,
                          aElemInd,aElem,aElemType);
  }
  void MakeElemSurroundingPoint(
      std::vector<unsigned int>& elsup_ind,
      std::vector<unsigned int>& elsup,
      const int nPo) const {
    ::delfem2::JArray_ElSuP_MeshMix(elsup_ind,elsup,
                                        aElemInd,aElem,nPo);
  }
  void MakeSurroundingRelationship(
      std::vector<int>& aElemFaceInd,
      std::vector<int>& aElemFaceRel,
      //
      const std::vector<unsigned int>& elsup_ind,
      const std::vector<unsigned int>& elsup) const {
    ::delfem2::ElSuEl_MeshMix(aElemFaceInd, aElemFaceRel,
                                           aElemInd,aElem,aElemType,
                                           elsup_ind,elsup);
  }
  void MakeSurroundingRelationship(
      std::vector<int>& aElemFaceInd,
      std::vector<int>& aElemFaceRel,
      const int nPo) const {
    ::delfem2::ElSuEl_MeshMix(aElemFaceInd, aElemFaceRel,
                                           aElemInd,aElem,aElemType,nPo);
  }
  int nElem() const {  return (int)aElemInd.size()-1; }
  void makeBoundary(
      CElemMixed& emb,
      const std::vector<int>& aElemFaceInd,
      const std::vector<int>& aElemFaceRel) const {
    ::delfem2::Boundary_MeshMix(emb.aElemInd, emb.aElem, emb.aElemType,
                            aElemInd, aElem, aElemType,
                            aElemFaceInd, aElemFaceRel);
  }
  void makeBoundary(CElemMixed& emb, int nPo ) const{
    std::vector<unsigned int> elsup_ind, elsup;
    this->MakeElemSurroundingPoint(elsup_ind, elsup, nPo);
    std::vector<int> aElemFaceInd, aElemFaceRel;
    this->MakeSurroundingRelationship(aElemFaceInd, aElemFaceRel,
                                   elsup_ind, elsup);
    this->makeBoundary(emb,
                    aElemFaceInd, aElemFaceRel);
  }
  void MakeGroupElem(
      int& ngroup,
      std::vector<int>& aIndGroup,
      const std::vector<int>& aElemFaceInd,
      const std::vector<int>& aElemFaceRel) const {
    ::delfem2::MakeGroupElem(ngroup, aIndGroup,
                             aElemInd, aElem,
                             aElemFaceInd, aElemFaceRel);
  }
  void MakeGroupElem(
      int& ngroup,
      std::vector<int>& aIndGroup,
      int nPo) const{
    ::delfem2::MakeGroupElem_MeshMix(ngroup, aIndGroup,
                           aElemInd, aElem, aElemType, nPo);
  }
  void ClipGroup(
      CElemMixed& em,
      int igroup,
      const std::vector<int>& aIndGroup) const{
    ::delfem2::ClipGroup_MeshMix(em.aElemInd,em.aElem,em.aElemType,
                       aElemInd,aElem,aElemType,
                       igroup,aIndGroup);
  }
  void FlipElement(std::vector<int>& aElem_Flip) const{
    ::delfem2::FlipElement_MeshMix(aElem_Flip,
                           aElemInd,aElem,aElemType);
  }
  void getTriElement(std::vector<unsigned int>& aTri) const{
    ::delfem2::Convert2Tri_MeshMix(aTri,
                           aElemInd,aElem,aElemType);
  }
private:
public:
  std::vector<unsigned int> aElemInd;
  std::vector<unsigned int> aElem;
  std::vector<delfem2::MESHELEM_TYPE> aElemType;
};
  
} // end namespace delfem2

#ifdef DFM2_HEADER_ONLY
#  include "delfem2/mshtopo.cpp"
#endif
 
#endif /* meshtopo_hpp */
