/*
 * Copyright (c) 2020 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */

/**
 * @file functions to analyze mesh topology for mixed meshes
 * @details the functions only care about the topology. Geometry (coordinate) information is not handled in this file
 */

// (2020/12/00) created

#ifndef DFM2_MSHMIX_H
#define DFM2_MSHMIX_H

#include "delfem2/mshelm.h"
#include "delfem2/dfm2_inline.h"
#include <vector>

namespace delfem2
{


DFM2_INLINE void AddElement(
                            const delfem2::MESHELEM_TYPE& femelem_type,
                            const std::vector<int>& aElemIn,
                            //
                            std::vector<unsigned int>& aElemInd,
                            std::vector<unsigned int>& aElem,
                            std::vector<delfem2::MESHELEM_TYPE>& aElemType);

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

DFM2_INLINE void MakeGroupElem_MeshMix(
                                       int& ngroup,
                                       std::vector<int>& aIndGroup,
                                       //
                                       const std::vector<unsigned int>& aElemInd,
                                       const std::vector<unsigned int>& aElem,
                                       const std::vector<delfem2::MESHELEM_TYPE>& aElemType,
                                       int nPo);


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

DFM2_INLINE void Convert2Tri_MeshMix(
                                     std::vector<unsigned int>& aTri,
                                     //
                                     const std::vector<unsigned int>& aElemInd,
                                     const std::vector<unsigned int>& aElem,
                                     const std::vector<delfem2::MESHELEM_TYPE>& aElemType);

DFM2_INLINE void FlipElement_MeshMix(
                                     std::vector<int>& aElem_Flip,
                                     //
                                     const std::vector<unsigned int>& aElemInd,
                                     const std::vector<unsigned int>& aElem,
                                     const std::vector<delfem2::MESHELEM_TYPE>& aElemType);


DFM2_INLINE void MakeGroupElem(
                               int& ngroup,
                               std::vector<int>& aIndGroup,
                               const std::vector<unsigned int>& aElemInd,
                               const std::vector<unsigned int>& aElem,
                               const std::vector<int>& aElemFaceInd,
                               const std::vector<int>& aElemFaceRel);

DFM2_INLINE void MarkConnectedElements(
                                       std::vector<int>& aIndGroup,
                                       unsigned int itri_ker,
                                       int igroup,
                                       const std::vector<int>& aElemFaceInd,
                                       const std::vector<int>& aElemFaceRel);

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

#ifndef DFM2_STATIC_LIBRARY
#  include "delfem2/mshmix.cpp"
#endif

#endif /* mshmix_h */
