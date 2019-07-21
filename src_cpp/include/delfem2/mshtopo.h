/*
 * Copyright (c) 2019 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */

#ifndef MSHTOPO_H
#define MSHTOPO_H

#include <stdio.h>
#include <vector>

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
const int mapMeshElemType2NNodeElem[7] = {
  3, // TRI
  4, // TET
  4, // QUAD
  8, // HEX
  5, // PYRAM
  6, // WEDGE
  8  // VOX
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
const int noelElemFace_Vox[8][4] = { // this numbering is corresponds to VTK_VOX
  { 0, 4, 6, 2 }, // -x
  { 1, 3, 7, 5 }, // +x
  { 0, 1, 5, 4 }, // -y
  { 2, 6, 7, 3 }, // +y
  { 0, 2, 3, 1 }, // -z
  { 4, 5, 7, 6 }  // +z
};
const int noelElemFace_Hex[8][4] = { // this numbering is corresponds to VTK_HEX
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
const int noelElemEdge_Vox[12][2] = {
  {0,1},{3,2},{4,5},{7,6},
  {0,3},{1,2},{4,7},{5,6},
  {0,4},{1,5},{3,7},{2,6} 
};

/////////////////////////////////////////////
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

void JArray_Sort(const std::vector<int>& index,
                 std::vector<int>& array);
void JArray_Sort(const int* index, const int size,
                 int* array);

void JArray_AddDiagonal(std::vector<int >& psup_ind1,
                      std::vector<int >& psup1,
                      const int* psup_ind0, int npsup_ind0,
                      const int* psup0, int npsup0);

void JArray_Print(const std::vector<int>& index,
                        const std::vector<int>& array);

void JArray_AddMasterSlavePattern(std::vector<int>& index,
                           std::vector<int>& array,
                           const int*aMSFlag,
                           int ndim,
                           const int* psup_ind0,
                           const int npsup_ind0,
                           const int* psup);

//////////////////////////////////////
void AddElement(const MESHELEM_TYPE& femelem_type,
                const std::vector<int>& aElemIn,
                ////
                std::vector<int>& aElemInd,
                std::vector<int>& aElem,
                std::vector<MESHELEM_TYPE>& aElemType);

////
void convert2Tri_Quad(std::vector<unsigned int>& aTri,
                      const std::vector<unsigned int>& aQuad);
void convert2Tri(std::vector<int>& aTri,
                 ////
                 const std::vector<int>& aElemInd,
                 const std::vector<int>& aElem,
                 const std::vector<MESHELEM_TYPE>& aElemType);

void ElemQuad_DihedralTri(std::vector<unsigned int>& aQuad,
                          const unsigned int* aTri, int nTri,
                          int np);
////
void FlipElement_Tri(std::vector<int>& aTri);
void FlipElement(std::vector<int>& aElem_Flip,
          ////
          const std::vector<int>& aElemInd,
          const std::vector<int>& aElem,
          const std::vector<MESHELEM_TYPE>& aElemType);

////////////////////////////////////////
// make elsup
void JArrayElemSurPoint_MeshElem(std::vector<int>& elsup_ind,
                              std::vector<int>& elsup,
                              ////
                              const unsigned int* aElem,
                              int nEl,
                              int nPoEl,
                              int nPo);
void JArrayElemSurPoint_MeshTri(std::vector<int>& elsup_ind,
                                std::vector<int>& elsup,
                                ////
                                const std::vector<unsigned int>& aTri,
                                int nXYZ);
void JArrayElemSurPoint_MeshMix(std::vector<int>& elsup_ind,
                              std::vector<int>& elsup,
                              ////
                              const std::vector<int>& aElemInd,
                              const std::vector<int>& aElem,
                              const int nPo);

////////////////////////////////////////
// elem sur elem
void makeSurroundingRelationship(std::vector<int>& aElSurRel,
                                 const unsigned int* aEl, int nEl, int nNoEl,
                                 const std::vector<int>& elsup_ind,
                                 const std::vector<int>& elsup,
                                 const int nfael,
                                 const int nnofa,
                                 const int noelElemFace[][4]);
void makeSurroundingRelationship(std::vector<int>& aElemSurRel,
                                 const unsigned int* aElem, int nEl,
                                 MESHELEM_TYPE type,
                                 const int nXYZ);
void makeSurroundingRelationship(std::vector<int>& aElemFaceInd,
                                 std::vector<int>& aElemFaceRel,
                                 const std::vector<int>& aElemInd,
                                 const std::vector<int>& aElem,
                                 const std::vector<MESHELEM_TYPE>& aElemType,
                                 const int nPo);
void makeSurroundingRelationship(std::vector<int>& aElemFaceInd,
                                 std::vector<int>& aElemFaceRel,
                                 const std::vector<int>& aElemInd,
                                 const std::vector<int>& aElem,
                                 const std::vector<MESHELEM_TYPE>& aElemType,
                                 const std::vector<int>& elsup_ind,
                                 const std::vector<int>& elsup);

//////////////////////////////////////////
// make boundary
void makeBoundary(std::vector<int>& aElemInd_Bound,
                  std::vector<int>& aElem_Bound,
                  std::vector<MESHELEM_TYPE>& aElemType_Bound,
                  ////
                  const std::vector<int>& aElemInd,
                  const std::vector<int>& aElem,
                  const std::vector<MESHELEM_TYPE>& aElemType,
                  const std::vector<int>& aElemFaceInd,
                  const std::vector<int>& aElemFaceRel);

////////////////////////////////////////
// make psup
// psup -> edge bidirectional
// edge unidir (ip0<ip1)
// line (array of 2)
void JArrayPointSurPoint_MeshOneRingNeighborhood(std::vector<int>& psup_ind,
                                                 std::vector<int>& psup,
                                                 ////
                                                 const unsigned int* pElem,
                                                 const std::vector<int>& elsup_ind,
                                                 const std::vector<int>& elsup,
                                                 int nnoel,
                                                 int nnode);
void JArrayPointSurPoint_MeshOneRingNeighborhood(std::vector<int>& psup_ind,
                                                 std::vector<int>& psup,
                                                 ////
                                                 const unsigned int* pElem,
                                                 int nEl,
                                                 int nPoEl,
                                                 int nPo);
void makeOneRingNeighborhood_TriFan(std::vector<int>& psup_ind,
                                    std::vector<int>& psup,
                                    ////
                                    const std::vector<int>& aTri,
                                    const std::vector<int>& aTriSurRel,
                                    const std::vector<int>& elsup_ind,
                                    const std::vector<int>& elsup,
                                    int np);
void JArrayEdgeUnidir_PointSurPoint(std::vector<int>& edge_ind,
                                std::vector<int>& edge,
                                /////
                                const std::vector<int>& psup_ind,
                                const std::vector<int>& psup);
void JArrayEdge_MeshElem(std::vector<int>& edge_ind,
                         std::vector<int>& edge,
                         ////
                         const unsigned int* aElm0,
                         MESHELEM_TYPE elem_type,
                         const std::vector<int>& elsup_ind,
                         const std::vector<int>& elsup,
                         bool is_bidirectional);
void MeshLine_JArrayEdge(std::vector<unsigned int>& aLine,
                         ////
                         const std::vector<int>& psup_ind,
                         const std::vector<int>& psup);
void MeshLine_MeshElem(std::vector<unsigned int>& aLine,
                       const unsigned int* aElm0,
                       unsigned int nElem,
                       MESHELEM_TYPE elem_type,
                       unsigned int nPo);

//////////////////////////////////////
void MarkConnectedElements(std::vector<int>& aIndGroup,
                            int itri_ker,
                            int igroup,
                            const std::vector<int>& aTriSurRel,
                            const int nfael);
void MakeGroupElem(int& ngroup,
                   std::vector<int>& aIndGroup,
                   const std::vector<int>& aElem,
                   const std::vector<int>& aElemSurRel,
                   const int nfael,
                   const int nnoel);
void MakeGroupElem_Tri(int& ngroup,
                       std::vector<int>& aIndGroup,
                       const std::vector<int>& aTri,
                       const std::vector<int>& aTriSurRel);
void MakeGroupElem(int& ngroup,
                   std::vector<int>& aIndGroup,
                   const std::vector<int>& aElemInd,
                   const std::vector<int>& aElem,
                   const std::vector<int>& aElemFaceInd,
                   const std::vector<int>& aElemFaceRel);
void MakeGroupElem(int& ngroup,
                   std::vector<int>& aIndGroup,
                   /////
                   const std::vector<int>& aElemInd,
                   const std::vector<int>& aElem,
                   const std::vector<MESHELEM_TYPE>& aElemType,
                   int nPo);
void ClipGroup(std::vector<int>& aElemInd1,
               std::vector<int>& aElem1,
               std::vector<MESHELEM_TYPE>& aElemType1,
               ///
               const std::vector<int>& aElemInd,
               const std::vector<int>& aElem,
               const std::vector<MESHELEM_TYPE>& aElemType,
               int igroup,
               const std::vector<int>& aIndGroup);

/////////////////////////////////////////////
void QuadSubdiv(std::vector<unsigned int>& aQuad1,
                std::vector<int>& psup_ind,
                std::vector<int>& psup,
                std::vector<int>& aEdgeFace0,
                const unsigned int* aQuad0, int nQuad0,
                unsigned int nPo0);
void HexSubdiv(std::vector<unsigned int>& aHex1,
               std::vector<int>& psupIndHex0,
               std::vector<int>& psupHex0,
               std::vector<unsigned int>& aQuadHex0,
               const unsigned int* aHex0, int nHex0,
               const int nhp0);
void TetSubdiv(std::vector<unsigned int>& aTet1,
               std::vector<int>& psup_ind,
               std::vector<int>& psup,
               const unsigned int* aTet0, int nTet0,
               unsigned int nPoint0);
int findEdge(int ip0, int ip1,
             const std::vector<int>& psup_ind,
             const std::vector<int>& psup);
int findFace(int ipc0, int ipc1, int ip2, int ip3,
             const std::vector<unsigned int>& aQuad,
             const std::vector<int>& elsupInd,
             const std::vector<int>& elsup);

//////////////////////////////////////////


class CElemMixed{
public:
  CElemMixed(){
    aElemInd.resize(1,0);
  }
  void AddElement(const MESHELEM_TYPE& femelem_type, const std::vector<int>& aElemIn){
    ::AddElement(femelem_type,aElemIn,
                 aElemInd,aElem,aElemType);
  }
  void MakeElemSurroundingPoint(std::vector<int>& elsup_ind, std::vector<int>& elsup,
                                const int nPo){
    ::JArrayElemSurPoint_MeshMix(elsup_ind,elsup,
                               aElemInd,aElem,nPo);
  }
  void MakeSurroundingRelationship(std::vector<int>& aElemFaceInd,
                                   std::vector<int>& aElemFaceRel,
                                   ///
                                   const std::vector<int>& elsup_ind,
                                   const std::vector<int>& elsup) const {
    ::makeSurroundingRelationship(aElemFaceInd, aElemFaceRel,
                                  aElemInd,aElem,aElemType,
                                  elsup_ind,elsup);
  }
  void MakeSurroundingRelationship(std::vector<int>& aElemFaceInd,
                                   std::vector<int>& aElemFaceRel,
                                   const int nPo) const {
    ::makeSurroundingRelationship(aElemFaceInd, aElemFaceRel,
                                  aElemInd,aElem,aElemType,nPo);
  }
  int nElem() const {  return (int)aElemInd.size()-1; }
  void makeBoundary(CElemMixed& emb,
                    const std::vector<int>& aElemFaceInd,
                    const std::vector<int>& aElemFaceRel){
    ::makeBoundary(emb.aElemInd, emb.aElem, emb.aElemType,
                   aElemInd, aElem, aElemType,
                   aElemFaceInd, aElemFaceRel);
  }
  void makeBoundary(CElemMixed& emb, int nPo){
    std::vector<int> elsup_ind, elsup;
    this->MakeElemSurroundingPoint(elsup_ind, elsup, nPo);
    std::vector<int> aElemFaceInd, aElemFaceRel;
    this->MakeSurroundingRelationship(aElemFaceInd, aElemFaceRel,
                                   elsup_ind, elsup);
    this->makeBoundary(emb,
                    aElemFaceInd, aElemFaceRel);
  }
  void MakeGroupElem(int& ngroup,
                     std::vector<int>& aIndGroup,
                     const std::vector<int>& aElemFaceInd,
                     const std::vector<int>& aElemFaceRel){
    ::MakeGroupElem(ngroup, aIndGroup,
                    aElemInd, aElem,
                    aElemFaceInd, aElemFaceRel);
  }
  void MakeGroupElem(int& ngroup,
                     std::vector<int>& aIndGroup,
                     int nPo){
    ::MakeGroupElem(ngroup, aIndGroup,
                    aElemInd, aElem, aElemType, nPo);
  }
  void ClipGroup(CElemMixed& em,
                 int igroup,
                 const std::vector<int>& aIndGroup){
    ::ClipGroup(em.aElemInd,em.aElem,em.aElemType,
                aElemInd,aElem,aElemType,
                igroup,aIndGroup);
  }
  void FlipElement(std::vector<int>& aElem_Flip){
    ::FlipElement(aElem_Flip,
           aElemInd,aElem,aElemType);
  }
  void getTriElement(std::vector<int>& aTri){
    ::convert2Tri(aTri,
                  aElemInd,aElem,aElemType);
  }
private:
public:
  std::vector<int> aElemInd;
  std::vector<int> aElem;
  std::vector<MESHELEM_TYPE> aElemType;
};

#endif /* meshtopo_hpp */
