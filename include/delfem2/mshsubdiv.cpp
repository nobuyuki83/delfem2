/*
 * Copyright (c) 2019 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */

#include "delfem2/mshsubdiv.h"

#include <vector>
#include <cassert>
#include <climits>

#include "delfem2/mshuni.h"

// ----------------------------------------------------

DFM2_INLINE unsigned int delfem2::findEdge(
    unsigned int ip0,
    unsigned int ip1,
    const std::vector<unsigned int> &psup_ind,
    const std::vector<unsigned int> &psup) {
  if (ip1 > ip0) {
    for (unsigned int ipsup = psup_ind[ip0]; ipsup < psup_ind[ip0 + 1]; ++ipsup) {
      unsigned int ip2 = psup[ipsup];
      if (ip2 == ip1) { return ipsup; }
    }
  } else {
    for (unsigned int ipsup = psup_ind[ip1]; ipsup < psup_ind[ip1 + 1]; ++ipsup) {
      unsigned int ip2 = psup[ipsup];
      if (ip2 == ip0) { return ipsup; }
    }
  }
  return UINT_MAX;
}

DFM2_INLINE int delfem2::findFace(
    unsigned int ip0,
    unsigned int ip1,
    unsigned int ip2,
    unsigned int ip3,
    const std::vector<unsigned int> &quad_vtx_idx,
    const std::vector<unsigned int> &elsup_ind,
    const std::vector<unsigned int> &elsup) {
  if (ip0 >= elsup_ind.size() - 1) return -1;
  for (unsigned int ielsup = elsup_ind[ip0]; ielsup < elsup_ind[ip0 + 1]; ++ielsup) {
    unsigned int ie0 = elsup[ielsup];
    unsigned int iq0 = quad_vtx_idx[ie0 * 4 + 0];
    unsigned int iq1 = quad_vtx_idx[ie0 * 4 + 1];
    unsigned int iq2 = quad_vtx_idx[ie0 * 4 + 2];
    unsigned int iq3 = quad_vtx_idx[ie0 * 4 + 3];
    if (ip0 != iq0 && ip0 != iq1 && ip0 != iq2 && ip0 != iq3) continue;
    if (ip1 != iq0 && ip1 != iq1 && ip1 != iq2 && ip1 != iq3) continue;
    if (ip2 != iq0 && ip2 != iq1 && ip2 != iq2 && ip2 != iq3) continue;
    if (ip3 != iq0 && ip3 != iq1 && ip3 != iq2 && ip3 != iq3) continue;
    return (int) ie0;
  }
  return -1;
}

// new points is in the order of [old points], [edge points], [face points]
DFM2_INLINE void delfem2::SubdivTopo_MeshQuad(
    std::vector<unsigned int> &quad_vtxidx1,
    std::vector<unsigned int> &psup_ind,
    std::vector<unsigned int> &psup,
    std::vector<unsigned int> &aEdgeFace0, // two points on the edge and two quads touching the edge
    const unsigned int *quad_vtxidx0,
    size_t num_quad0,
    size_t num_vtx) {
  const size_t nq0 = num_quad0;
  const auto np0 = static_cast<unsigned int>(num_vtx);
  std::vector<unsigned int> elsup_ind, elsup;
  JArray_ElSuP_MeshElem(
      elsup_ind, elsup,
      quad_vtxidx0, num_quad0, 4,
      num_vtx);
  JArrayEdge_MeshElem(
      psup_ind, psup,
      quad_vtxidx0, MESHELEM_QUAD,
      elsup_ind, elsup,
      false); // is_bidirectional = false
  const auto ne0 = static_cast<unsigned int>(psup.size());
  aEdgeFace0.resize(0);
  aEdgeFace0.reserve(ne0 * 4);
  for (unsigned int ip = 0; ip < num_vtx; ++ip) {
    for (unsigned int ipsup = psup_ind[ip]; ipsup < psup_ind[ip + 1]; ++ipsup) {
      const unsigned int ip1 = psup[ipsup];
      aEdgeFace0.push_back(ip);
      aEdgeFace0.push_back(ip1);
      unsigned int iq0 = UINT_MAX, iq1 = UINT_MAX;
      for (unsigned int ielsup = elsup_ind[ip]; ielsup < elsup_ind[ip + 1]; ++ielsup) {
        const unsigned int jq0 = elsup[ielsup];
        const unsigned int jp0 = quad_vtxidx0[jq0 * 4 + 0];
        const unsigned int jp1 = quad_vtxidx0[jq0 * 4 + 1];
        const unsigned int jp2 = quad_vtxidx0[jq0 * 4 + 2];
        const unsigned int jp3 = quad_vtxidx0[jq0 * 4 + 3];
        if ((jp0 != ip) && (jp1 != ip) && (jp2 != ip) && (jp3 != ip)) { continue; }
        if ((jp0 != ip1) && (jp1 != ip1) && (jp2 != ip1) && (jp3 != ip1)) { continue; }
        // ----------------------------
        if (iq0 == UINT_MAX) { iq0 = jq0; }
        else {
          assert(iq1 == UINT_MAX);
          iq1 = jq0;
        }
      }
      aEdgeFace0.push_back(iq0);
      aEdgeFace0.push_back(iq1);
    }
  }
  quad_vtxidx1.resize(0);
  quad_vtxidx1.reserve(num_quad0 * 4);
  for (unsigned int iq = 0; iq < nq0; ++iq) {
    const unsigned int ip0 = quad_vtxidx0[iq * 4 + 0];
    const unsigned int ip1 = quad_vtxidx0[iq * 4 + 1];
    const unsigned int ip2 = quad_vtxidx0[iq * 4 + 2];
    const unsigned int ip3 = quad_vtxidx0[iq * 4 + 3];
    const unsigned int ie01 = findEdge(ip0, ip1, psup_ind, psup);
    assert(ie01 != UINT_MAX);
    const unsigned int ie12 = findEdge(ip1, ip2, psup_ind, psup);
    assert(ie12 != UINT_MAX);
    const unsigned int ie23 = findEdge(ip2, ip3, psup_ind, psup);
    assert(ie23 != UINT_MAX);
    const unsigned int ie30 = findEdge(ip3, ip0, psup_ind, psup);
    assert(ie30 != UINT_MAX);
    const unsigned int ip01 = ie01 + np0;
    const unsigned int ip12 = ie12 + np0;
    const unsigned int ip23 = ie23 + np0;
    const unsigned int ip30 = ie30 + np0;
    const unsigned int ip0123 = iq + np0 + ne0;
    quad_vtxidx1.push_back(ip0);
    quad_vtxidx1.push_back(ip01);
    quad_vtxidx1.push_back(ip0123);
    quad_vtxidx1.push_back(ip30);
    quad_vtxidx1.push_back(ip1);
    quad_vtxidx1.push_back(ip12);
    quad_vtxidx1.push_back(ip0123);
    quad_vtxidx1.push_back(ip01);
    quad_vtxidx1.push_back(ip2);
    quad_vtxidx1.push_back(ip23);
    quad_vtxidx1.push_back(ip0123);
    quad_vtxidx1.push_back(ip12);
    quad_vtxidx1.push_back(ip3);
    quad_vtxidx1.push_back(ip30);
    quad_vtxidx1.push_back(ip0123);
    quad_vtxidx1.push_back(ip23);
  }
}


// new points is in the order of [old points], [edge points]
DFM2_INLINE void delfem2::SubdivTopo_MeshTet(
    std::vector<unsigned int> &aTet1,
    std::vector<unsigned int> &psup_ind,
    std::vector<unsigned int> &psup,
    const unsigned int *aTet0,
    int nTet0,
    unsigned int nPoint0) {
  const int nt0 = nTet0;
  std::vector<unsigned int> elsup_ind, elsup;
  JArray_ElSuP_MeshElem(
      elsup_ind, elsup,
      aTet0, nTet0, 4,
      nPoint0);
  JArrayEdge_MeshElem(
      psup_ind, psup,
      aTet0, MESHELEM_TET, elsup_ind, elsup,
      false);
  aTet1.resize(0);
  aTet1.reserve(nTet0 * 4);
  for (int it = 0; it < nt0; ++it) {
    unsigned int ip0 = aTet0[it * 4 + 0];
    unsigned int ip1 = aTet0[it * 4 + 1];
    unsigned int ip2 = aTet0[it * 4 + 2];
    unsigned int ip3 = aTet0[it * 4 + 3];
    const unsigned int ie01 = findEdge(ip0, ip1, psup_ind, psup);
    assert(ie01 != UINT_MAX);
    const unsigned int ie02 = findEdge(ip0, ip2, psup_ind, psup);
    assert(ie02 != UINT_MAX);
    const unsigned int ie03 = findEdge(ip0, ip3, psup_ind, psup);
    assert(ie03 != UINT_MAX);
    const unsigned int ie12 = findEdge(ip1, ip2, psup_ind, psup);
    assert(ie12 != UINT_MAX);
    const unsigned int ie13 = findEdge(ip1, ip3, psup_ind, psup);
    assert(ie13 != UINT_MAX);
    const unsigned int ie23 = findEdge(ip2, ip3, psup_ind, psup);
    assert(ie23 != UINT_MAX);
    unsigned int ip01 = ie01 + nPoint0;
    unsigned int ip02 = ie02 + nPoint0;
    unsigned int ip03 = ie03 + nPoint0;
    unsigned int ip12 = ie12 + nPoint0;
    unsigned int ip13 = ie13 + nPoint0;
    unsigned int ip23 = ie23 + nPoint0;
    aTet1.push_back(ip0);
    aTet1.push_back(ip01);
    aTet1.push_back(ip02);
    aTet1.push_back(ip03);
    aTet1.push_back(ip1);
    aTet1.push_back(ip01);
    aTet1.push_back(ip13);
    aTet1.push_back(ip12);
    aTet1.push_back(ip2);
    aTet1.push_back(ip02);
    aTet1.push_back(ip12);
    aTet1.push_back(ip23);
    aTet1.push_back(ip3);
    aTet1.push_back(ip03);
    aTet1.push_back(ip23);
    aTet1.push_back(ip13);
    aTet1.push_back(ip01);
    aTet1.push_back(ip23);
    aTet1.push_back(ip13);
    aTet1.push_back(ip12);
    aTet1.push_back(ip01);
    aTet1.push_back(ip23);
    aTet1.push_back(ip12);
    aTet1.push_back(ip02);
    aTet1.push_back(ip01);
    aTet1.push_back(ip23);
    aTet1.push_back(ip02);
    aTet1.push_back(ip03);
    aTet1.push_back(ip01);
    aTet1.push_back(ip23);
    aTet1.push_back(ip03);
    aTet1.push_back(ip13);
  }
}

/*
// TODO: This one is imcomplete
void VoxSubdiv
(std::vector<int>& aVox1,
 std::vector<int>& psupIndHex0,
 std::vector<int>& psupHex0,
 std::vector<int>& aQuadHex0,
 ///
 const std::vector<int>& aVox0,
 const int nhp0)
{
  //  int nhp0 = (int)aHexPoint0.size(); // hex point
  std::vector<int> elsupIndHex0, elsupHex0;
  makeElemSurroundingPoint(elsupIndHex0, elsupHex0,
                           aVox0,8,nhp0);
  
  //edge
  makeEdgeVox(psupIndHex0, psupHex0,
              aVox0, elsupIndHex0,elsupHex0, nhp0);
  
  //face
  aQuadHex0.clear();
  {
    std::vector<int> aHexSurRel0;
    makeSurroundingRelationship(aHexSurRel0,
                                aVox0,FEMELEM_VOX,
                                elsupIndHex0,elsupHex0);
    for(int ih=0;ih<aVox0.size()/8;++ih){
      for(int ifh=0;ifh<6;++ifh){
        int jh0 = aHexSurRel0[ih*6*2+ifh*2+0];
        if( jh0!=-1 && ih>jh0 ) continue;
        for(int inofa=0;inofa<4;++inofa){
          int inoel0 = noelElemFace_Hex[ifh][inofa];
          int igp0 = aVox0[ih*8+inoel0];
          aQuadHex0.push_back(igp0);
        }
      }
    }
  }
  std::vector<int> elsupIndQuadHex0, elsupQuadHex0;
  makeElemSurroundingPoint(elsupIndQuadHex0,elsupQuadHex0,
                           aQuadHex0,4,nhp0);
  
  const int neh0 = (int)psupHex0.size();
  const int nfh0 = (int)aQuadHex0.size()/4;
  std::cout << nfh0 << " " << aQuadHex0.size() << std::endl;
  
  // making vox
  aVox1.clear();
  for(int ih=0;ih<aVox0.size()/8;++ih){
    int ihc0 = aVox0[ih*8+0];
    int ihc1 = aVox0[ih*8+1];
    int ihc2 = aVox0[ih*8+2];
    int ihc3 = aVox0[ih*8+3];
    int ihc4 = aVox0[ih*8+4];
    int ihc5 = aVox0[ih*8+5];
    int ihc6 = aVox0[ih*8+6];
    int ihc7 = aVox0[ih*8+7];
    int ihc01 = findEdge(ihc0,ihc1, psupIndHex0,psupHex0)+nhp0; assert(ihc01>=nhp0&&ihc01<nhp0+neh0);
    int ihc32 = findEdge(ihc3,ihc2, psupIndHex0,psupHex0)+nhp0; assert(ihc32>=nhp0&&ihc32<nhp0+neh0);
    int ihc45 = findEdge(ihc4,ihc5, psupIndHex0,psupHex0)+nhp0; assert(ihc45>=nhp0&&ihc45<nhp0+neh0);
    int ihc76 = findEdge(ihc7,ihc6, psupIndHex0,psupHex0)+nhp0; assert(ihc76>=nhp0&&ihc76<nhp0+neh0);
    int ihc03 = findEdge(ihc0,ihc3, psupIndHex0,psupHex0)+nhp0; assert(ihc03>=nhp0&&ihc03<nhp0+neh0);
    int ihc12 = findEdge(ihc1,ihc2, psupIndHex0,psupHex0)+nhp0; assert(ihc12>=nhp0&&ihc12<nhp0+neh0);
    int ihc47 = findEdge(ihc4,ihc7, psupIndHex0,psupHex0)+nhp0; assert(ihc47>=nhp0&&ihc47<nhp0+neh0);
    int ihc56 = findEdge(ihc5,ihc6, psupIndHex0,psupHex0)+nhp0; assert(ihc56>=nhp0&&ihc56<nhp0+neh0);
    int ihc04 = findEdge(ihc0,ihc4, psupIndHex0,psupHex0)+nhp0; assert(ihc04>=nhp0&&ihc04<nhp0+neh0);
    int ihc15 = findEdge(ihc1,ihc5, psupIndHex0,psupHex0)+nhp0; assert(ihc15>=nhp0&&ihc15<nhp0+neh0);
    int ihc37 = findEdge(ihc3,ihc7, psupIndHex0,psupHex0)+nhp0; assert(ihc37>=nhp0&&ihc37<nhp0+neh0);
    int ihc26 = findEdge(ihc2,ihc6, psupIndHex0,psupHex0)+nhp0; assert(ihc26>=nhp0&&ihc26<nhp0+neh0);
    int ihc0462 = findFace(ihc0,ihc4,ihc6,ihc2, aQuadHex0,elsupIndQuadHex0,elsupQuadHex0)+nhp0+neh0; assert(ihc0462>=nhp0+neh0&&ihc0462<nhp0+neh0+nfh0);
    int ihc1375 = findFace(ihc1,ihc3,ihc7,ihc5, aQuadHex0,elsupIndQuadHex0,elsupQuadHex0)+nhp0+neh0; assert(ihc1375>=nhp0+neh0&&ihc1375<nhp0+neh0+nfh0);
    int ihc0154 = findFace(ihc0,ihc1,ihc5,ihc4, aQuadHex0,elsupIndQuadHex0,elsupQuadHex0)+nhp0+neh0; assert(ihc0154>=nhp0+neh0&&ihc0154<nhp0+neh0+nfh0);
    int ihc2673 = findFace(ihc2,ihc6,ihc7,ihc3, aQuadHex0,elsupIndQuadHex0,elsupQuadHex0)+nhp0+neh0; assert(ihc2673>=nhp0+neh0&&ihc2673<nhp0+neh0+nfh0);
    int ihc0231 = findFace(ihc0,ihc2,ihc3,ihc1, aQuadHex0,elsupIndQuadHex0,elsupQuadHex0)+nhp0+neh0; assert(ihc0231>=nhp0+neh0&&ihc0231<nhp0+neh0+nfh0);
    int ihc4576 = findFace(ihc4,ihc5,ihc7,ihc6, aQuadHex0,elsupIndQuadHex0,elsupQuadHex0)+nhp0+neh0; assert(ihc4576>=nhp0+neh0&&ihc4576<nhp0+neh0+nfh0);
    int ihc0473 = findFace(ihc0,ihc4,ihc7,ihc3, aQuadHex0,elsupIndQuadHex0,elsupQuadHex0)+nhp0+neh0; assert(ihc0473>=nhp0+neh0&&ihc0473<nhp0+neh0+nfh0);
    int ihc01234567 = ih + nhp0 + neh0 + nfh0;
    aVox1.push_back(ihc0);
    aVox1.push_back(ihc01);
    aVox1.push_back(ihc02);
    aVox1.push_back(ihc03);
    aVox1.push_back(ihc04);
    aVox1.push_back(ihc0154);
    aVox1.push_back(ihc01234567);
    aVox1.push_back(ihc0473);
    /////
    aVox1.push_back(ihc01);
    aVox1.push_back(ihc1);
    aVox1.push_back(ihc12);
    aVox1.push_back(ihc0321);
    aVox1.push_back(ihc0154);
    aVox1.push_back(ihc15);
    aVox1.push_back(ihc1265);
    aVox1.push_back(ihc01234567);
    /////
    aVox1.push_back(ihc0321);
    aVox1.push_back(ihc12);
    aVox1.push_back(ihc2);
    aVox1.push_back(ihc32);
    aVox1.push_back(ihc01234567);
    aVox1.push_back(ihc1265);
    aVox1.push_back(ihc26);
    aVox1.push_back(ihc3762);
    /////
    aVox1.push_back(ihc03);
    aVox1.push_back(ihc0321);
    aVox1.push_back(ihc32);
    aVox1.push_back(ihc3);
    aVox1.push_back(ihc0473);
    aVox1.push_back(ihc01234567);
    aVox1.push_back(ihc3762);
    aVox1.push_back(ihc37);
    /////
    aVox1.push_back(ihc04);
    aVox1.push_back(ihc0154);
    aVox1.push_back(ihc01234567);
    aVox1.push_back(ihc0473);
    aVox1.push_back(ihc4);
    aVox1.push_back(ihc45);
    aVox1.push_back(ihc4567);
    aVox1.push_back(ihc47);
    ////
    aVox1.push_back(ihc0154);
    aVox1.push_back(ihc15);
    aVox1.push_back(ihc1265);
    aVox1.push_back(ihc01234567);
    aVox1.push_back(ihc45);
    aVox1.push_back(ihc5);
    aVox1.push_back(ihc56);
    aVox1.push_back(ihc4567);
    /////
    aVox1.push_back(ihc01234567);
    aVox1.push_back(ihc1265);
    aVox1.push_back(ihc26);
    aVox1.push_back(ihc3762);
    aVox1.push_back(ihc4567);
    aVox1.push_back(ihc56);
    aVox1.push_back(ihc6);
    aVox1.push_back(ihc76);
    /////
    aVox1.push_back(ihc0473);
    aVox1.push_back(ihc01234567);
    aVox1.push_back(ihc3762);
    aVox1.push_back(ihc37);
    aVox1.push_back(ihc47);
    aVox1.push_back(ihc4567);
    aVox1.push_back(ihc76);
    aVox1.push_back(ihc7);
  }
}
*/


// -------------------------------------

void delfem2::SubdivTopo_MeshHex(
    std::vector<unsigned int> &aHex1,
    std::vector<unsigned int> &psupIndHex0,
    std::vector<unsigned int> &psupHex0,
    std::vector<unsigned int> &aQuadHex0,
    //
    const unsigned int *aHex0,
    size_t nHex0,
    const size_t nHexPoint0) {
  const auto nhp0 = static_cast<unsigned int>(nHexPoint0);
  std::vector<unsigned int> elsupIndHex0, elsupHex0;
  JArray_ElSuP_MeshElem(
      elsupIndHex0, elsupHex0,
      aHex0, nHex0, 8, nhp0);

  //edge
  JArrayEdge_MeshElem(
      psupIndHex0, psupHex0,
      aHex0, MESHELEM_HEX, elsupIndHex0, elsupHex0,
      false); // is_directional = false

  //face
  aQuadHex0.clear();
  {
    std::vector<unsigned int> aHexSuHex0;
    ElSuEl_MeshElem(
        aHexSuHex0,
        aHex0, nHex0, 8,
        elsupIndHex0, elsupHex0,
        nFaceElem(MESHELEM_HEX),
        nNodeElemFace(MESHELEM_HEX, 0),
        noelElemFace(MESHELEM_HEX));
    for (unsigned int ih = 0; ih < (unsigned int) nHex0; ++ih) {
      for (int ifh = 0; ifh < 6; ++ifh) {
        unsigned int jh0 = aHexSuHex0[ih * 6 + ifh];
        if (jh0 != UINT_MAX && ih > jh0) continue;
        for (int inofa = 0; inofa < 4; ++inofa) {
          int inoel0 = noelElemFace_Hex[ifh][inofa];
          unsigned int igp0 = aHex0[ih * 8 + inoel0];
          aQuadHex0.push_back(igp0);
        }
      }
    }
  }
  std::vector<unsigned int> elsupIndQuadHex0, elsupQuadHex0;
  JArray_ElSuP_MeshElem(
      elsupIndQuadHex0, elsupQuadHex0,
      aQuadHex0.data(), aQuadHex0.size() / 4,
      4, nhp0);

  const auto neh0 = static_cast<unsigned int>(psupHex0.size());
  const auto nfh0 = static_cast<unsigned int>(aQuadHex0.size() / 4);
//  std::cout << nfh0 << " " << aQuadHex0.size() << std::endl;

  /*
  const int aNoelEdge[12][2] = {
    {0,1},{1,2},{2,3},{3,0},
    {4,5},{5,6},{6,7},{7,4},
    {0,4},{1,5},{2,6},{3,7} };
   */
  /*
  const int noelElemFace_Hex[8][4] = { // this numbering is corresponds to VTK_VOXEL
    { 0, 4, 7, 3 }, // -x
    { 1, 2, 6, 5 }, // +x
    { 0, 1, 5, 4 }, // -y
    { 3, 7, 6, 2 }, // +y
    { 0, 3, 2, 1 }, // -z
    { 4, 5, 6, 7 }  // +z
  };
   */

  // making hex
  aHex1.clear();
  for (unsigned int ih = 0; ih < nHex0; ++ih) {
    unsigned int ihc0 = aHex0[ih * 8 + 0];
    unsigned int ihc1 = aHex0[ih * 8 + 1];
    unsigned int ihc2 = aHex0[ih * 8 + 2];
    unsigned int ihc3 = aHex0[ih * 8 + 3];
    unsigned int ihc4 = aHex0[ih * 8 + 4];
    unsigned int ihc5 = aHex0[ih * 8 + 5];
    unsigned int ihc6 = aHex0[ih * 8 + 6];
    unsigned int ihc7 = aHex0[ih * 8 + 7];
    const unsigned int ihc01 = findEdge(ihc0, ihc1, psupIndHex0, psupHex0) + nhp0;
    assert(ihc01 >= nhp0 && ihc01 < nhp0 + neh0);
    const unsigned int ihc12 = findEdge(ihc1, ihc2, psupIndHex0, psupHex0) + nhp0;
    assert(ihc12 >= nhp0 && ihc12 < nhp0 + neh0);
    const unsigned int ihc23 = findEdge(ihc2, ihc3, psupIndHex0, psupHex0) + nhp0;
    assert(ihc23 >= nhp0 && ihc23 < nhp0 + neh0);
    const unsigned int ihc30 = findEdge(ihc3, ihc0, psupIndHex0, psupHex0) + nhp0;
    assert(ihc30 >= nhp0 && ihc30 < nhp0 + neh0);
    const unsigned int ihc45 = findEdge(ihc4, ihc5, psupIndHex0, psupHex0) + nhp0;
    assert(ihc45 >= nhp0 && ihc45 < nhp0 + neh0);
    const unsigned int ihc56 = findEdge(ihc5, ihc6, psupIndHex0, psupHex0) + nhp0;
    assert(ihc56 >= nhp0 && ihc56 < nhp0 + neh0);
    const unsigned int ihc67 = findEdge(ihc6, ihc7, psupIndHex0, psupHex0) + nhp0;
    assert(ihc67 >= nhp0 && ihc67 < nhp0 + neh0);
    const unsigned int ihc74 = findEdge(ihc7, ihc4, psupIndHex0, psupHex0) + nhp0;
    assert(ihc74 >= nhp0 && ihc74 < nhp0 + neh0);
    const unsigned int ihc04 = findEdge(ihc0, ihc4, psupIndHex0, psupHex0) + nhp0;
    assert(ihc04 >= nhp0 && ihc04 < nhp0 + neh0);
    const unsigned int ihc15 = findEdge(ihc1, ihc5, psupIndHex0, psupHex0) + nhp0;
    assert(ihc15 >= nhp0 && ihc15 < nhp0 + neh0);
    const unsigned int ihc26 = findEdge(ihc2, ihc6, psupIndHex0, psupHex0) + nhp0;
    assert(ihc26 >= nhp0 && ihc26 < nhp0 + neh0);
    const unsigned int ihc37 = findEdge(ihc3, ihc7, psupIndHex0, psupHex0) + nhp0;
    assert(ihc37 >= nhp0 && ihc37 < nhp0 + neh0);
    unsigned int ihc0473 = findFace(ihc0, ihc4, ihc7, ihc3, aQuadHex0, elsupIndQuadHex0, elsupQuadHex0) + nhp0 + neh0;
    assert(ihc0473 >= nhp0 + neh0 && ihc0473 < nhp0 + neh0 + nfh0);
    unsigned int ihc1265 = findFace(ihc1, ihc2, ihc6, ihc5, aQuadHex0, elsupIndQuadHex0, elsupQuadHex0) + nhp0 + neh0;
    assert(ihc1265 >= nhp0 + neh0 && ihc1265 < nhp0 + neh0 + nfh0);
    unsigned int ihc0154 = findFace(ihc0, ihc1, ihc5, ihc4, aQuadHex0, elsupIndQuadHex0, elsupQuadHex0) + nhp0 + neh0;
    assert(ihc0154 >= nhp0 + neh0 && ihc0154 < nhp0 + neh0 + nfh0);
    unsigned int ihc3762 = findFace(ihc3, ihc7, ihc6, ihc2, aQuadHex0, elsupIndQuadHex0, elsupQuadHex0) + nhp0 + neh0;
    assert(ihc3762 >= nhp0 + neh0 && ihc3762 < nhp0 + neh0 + nfh0);
    unsigned int ihc0321 = findFace(ihc0, ihc3, ihc2, ihc1, aQuadHex0, elsupIndQuadHex0, elsupQuadHex0) + nhp0 + neh0;
    assert(ihc0321 >= nhp0 + neh0 && ihc0321 < nhp0 + neh0 + nfh0);
    unsigned int ihc4567 = findFace(ihc4, ihc5, ihc6, ihc7, aQuadHex0, elsupIndQuadHex0, elsupQuadHex0) + nhp0 + neh0;
    assert(ihc4567 >= nhp0 + neh0 && ihc4567 < nhp0 + neh0 + nfh0);
    unsigned int ihc01234567 = ih + nhp0 + neh0 + nfh0;
    //0
    aHex1.push_back(ihc0);
    aHex1.push_back(ihc01);
    aHex1.push_back(ihc0321);
    aHex1.push_back(ihc30);
    aHex1.push_back(ihc04);
    aHex1.push_back(ihc0154);
    aHex1.push_back(ihc01234567);
    aHex1.push_back(ihc0473);
    //1
    aHex1.push_back(ihc01);
    aHex1.push_back(ihc1);
    aHex1.push_back(ihc12);
    aHex1.push_back(ihc0321);
    aHex1.push_back(ihc0154);
    aHex1.push_back(ihc15);
    aHex1.push_back(ihc1265);
    aHex1.push_back(ihc01234567);
    //2
    aHex1.push_back(ihc0321);
    aHex1.push_back(ihc12);
    aHex1.push_back(ihc2);
    aHex1.push_back(ihc23);
    aHex1.push_back(ihc01234567);
    aHex1.push_back(ihc1265);
    aHex1.push_back(ihc26);
    aHex1.push_back(ihc3762);
    //3
    aHex1.push_back(ihc30);
    aHex1.push_back(ihc0321);
    aHex1.push_back(ihc23);
    aHex1.push_back(ihc3);
    aHex1.push_back(ihc0473);
    aHex1.push_back(ihc01234567);
    aHex1.push_back(ihc3762);
    aHex1.push_back(ihc37);
    //4
    aHex1.push_back(ihc04);
    aHex1.push_back(ihc0154);
    aHex1.push_back(ihc01234567);
    aHex1.push_back(ihc0473);
    aHex1.push_back(ihc4);
    aHex1.push_back(ihc45);
    aHex1.push_back(ihc4567);
    aHex1.push_back(ihc74);
    //5
    aHex1.push_back(ihc0154);
    aHex1.push_back(ihc15);
    aHex1.push_back(ihc1265);
    aHex1.push_back(ihc01234567);
    aHex1.push_back(ihc45);
    aHex1.push_back(ihc5);
    aHex1.push_back(ihc56);
    aHex1.push_back(ihc4567);
    //6
    aHex1.push_back(ihc01234567);
    aHex1.push_back(ihc1265);
    aHex1.push_back(ihc26);
    aHex1.push_back(ihc3762);
    aHex1.push_back(ihc4567);
    aHex1.push_back(ihc56);
    aHex1.push_back(ihc6);
    aHex1.push_back(ihc67);
    //7
    aHex1.push_back(ihc0473);
    aHex1.push_back(ihc01234567);
    aHex1.push_back(ihc3762);
    aHex1.push_back(ihc37);
    aHex1.push_back(ihc74);
    aHex1.push_back(ihc4567);
    aHex1.push_back(ihc67);
    aHex1.push_back(ihc7);
  }
}

// TODO: make this handle open surface (average face & edge independently)
void delfem2::SubdivisionPoints_QuadCatmullClark(
    std::vector<double> &vtx_xyz1,
    // ------------------------
    const std::vector<unsigned int> &quad_vtxidx1,
    const std::vector<unsigned int> &aEdgeFace0,
    const std::vector<unsigned int> &psupIndQuad0,
    const std::vector<unsigned int> &psupQuad0,
    const unsigned int *quad_vtxidx0,
    size_t num_quad0,
    const double *vtx_xyz0,
    size_t num_vtx0) {
  const std::size_t nv0 = num_vtx0;
  const std::size_t ne0 = psupQuad0.size();
  const size_t nq0 = num_quad0;
  assert(aEdgeFace0.size() == ne0 * 4);
  const std::size_t nv1 = nv0 + ne0 + nq0;
  vtx_xyz1.resize(nv1 * 3);
  std::vector<unsigned int> aNFace(nv0, 0); // number of faces touching vertex
  for (unsigned int iv = 0; iv < nv0; ++iv) {
    vtx_xyz1[iv * 3 + 0] = 0;
    vtx_xyz1[iv * 3 + 1] = 0;
    vtx_xyz1[iv * 3 + 2] = 0;
  }
  for (unsigned int iq = 0; iq < nq0; ++iq) { // face
    const unsigned int iv0 = quad_vtxidx0[iq * 4 + 0];
    const unsigned int iv1 = quad_vtxidx0[iq * 4 + 1];
    const unsigned int iv2 = quad_vtxidx0[iq * 4 + 2];
    const unsigned int iv3 = quad_vtxidx0[iq * 4 + 3];
    const double p0x = (vtx_xyz0[iv0 * 3 + 0] + vtx_xyz0[iv1 * 3 + 0] + vtx_xyz0[iv2 * 3 + 0] + vtx_xyz0[iv3 * 3 + 0]) * 0.25;
    const double p0y = (vtx_xyz0[iv0 * 3 + 1] + vtx_xyz0[iv1 * 3 + 1] + vtx_xyz0[iv2 * 3 + 1] + vtx_xyz0[iv3 * 3 + 1]) * 0.25;
    const double p0z = (vtx_xyz0[iv0 * 3 + 2] + vtx_xyz0[iv1 * 3 + 2] + vtx_xyz0[iv2 * 3 + 2] + vtx_xyz0[iv3 * 3 + 2]) * 0.25;
    vtx_xyz1[(nv0 + ne0 + iq) * 3 + 0] = p0x;
    vtx_xyz1[(nv0 + ne0 + iq) * 3 + 1] = p0y;
    vtx_xyz1[(nv0 + ne0 + iq) * 3 + 2] = p0z;
    const unsigned int aIV[4] = {iv0, iv1, iv2, iv3};
    for (unsigned int jv0 : aIV) {
      vtx_xyz1[jv0 * 3 + 0] += p0x;
      vtx_xyz1[jv0 * 3 + 1] += p0y;
      vtx_xyz1[jv0 * 3 + 2] += p0z;
      aNFace[jv0] += 1;
    }
  }
  for (unsigned int ie = 0; ie < ne0; ++ie) { // edge
    const unsigned int iv0 = aEdgeFace0[ie * 4 + 0];
    const unsigned int iv1 = aEdgeFace0[ie * 4 + 1];
    const unsigned int iq0 = aEdgeFace0[ie * 4 + 2];
    const unsigned int iq1 = aEdgeFace0[ie * 4 + 3];
    if (iq1 != UINT_MAX) {
      const size_t iv1e = nv0 + ie;
      assert(iv1e < nv1);
      const size_t iv1q0 = nv0 + ne0 + iq0;
      assert(iv1q0 < nv1);
      const size_t iv1q1 = nv0 + ne0 + iq1;
      assert(iv1q1 < nv1);
      vtx_xyz1[iv1e * 3 + 0] =
          (vtx_xyz0[iv0 * 3 + 0] + vtx_xyz0[iv1 * 3 + 0] + vtx_xyz1[iv1q0 * 3 + 0] + vtx_xyz1[iv1q1 * 3 + 0]) * 0.25;
      vtx_xyz1[iv1e * 3 + 1] =
          (vtx_xyz0[iv0 * 3 + 1] + vtx_xyz0[iv1 * 3 + 1] + vtx_xyz1[iv1q0 * 3 + 1] + vtx_xyz1[iv1q1 * 3 + 1]) * 0.25;
      vtx_xyz1[iv1e * 3 + 2] =
          (vtx_xyz0[iv0 * 3 + 2] + vtx_xyz0[iv1 * 3 + 2] + vtx_xyz1[iv1q0 * 3 + 2] + vtx_xyz1[iv1q1 * 3 + 2]) * 0.25;
    } else {
      const size_t iv1e = nv0 + ie;
      assert(iv1e < nv1);
      const size_t iv1q0 = nv0 + ne0 + iq0;
      assert(iv1q0 < nv1);
      vtx_xyz1[iv1e * 3 + 0] = (vtx_xyz0[iv0 * 3 + 0] + vtx_xyz0[iv1 * 3 + 0]) * 0.5;
      vtx_xyz1[iv1e * 3 + 1] = (vtx_xyz0[iv0 * 3 + 1] + vtx_xyz0[iv1 * 3 + 1]) * 0.5;
      vtx_xyz1[iv1e * 3 + 2] = (vtx_xyz0[iv0 * 3 + 2] + vtx_xyz0[iv1 * 3 + 2]) * 0.5;
    }
    vtx_xyz1[iv0 * 3 + 0] += vtx_xyz0[iv0 * 3 + 0] + vtx_xyz0[iv1 * 3 + 0];
    vtx_xyz1[iv0 * 3 + 1] += vtx_xyz0[iv0 * 3 + 1] + vtx_xyz0[iv1 * 3 + 1];
    vtx_xyz1[iv0 * 3 + 2] += vtx_xyz0[iv0 * 3 + 2] + vtx_xyz0[iv1 * 3 + 2];
    vtx_xyz1[iv1 * 3 + 0] += vtx_xyz0[iv0 * 3 + 0] + vtx_xyz0[iv1 * 3 + 0];
    vtx_xyz1[iv1 * 3 + 1] += vtx_xyz0[iv0 * 3 + 1] + vtx_xyz0[iv1 * 3 + 1];
    vtx_xyz1[iv1 * 3 + 2] += vtx_xyz0[iv0 * 3 + 2] + vtx_xyz0[iv1 * 3 + 2];
  }
  for (unsigned int iv = 0; iv < nv0; ++iv) {
    const unsigned int nf = aNFace[iv]; // number of faces touching this vertex
    if (nf == 0) { continue; }
    // add face
    const double tmp0 = 1.0 / (nf * nf);
    vtx_xyz1[iv * 3 + 0] *= tmp0;
    vtx_xyz1[iv * 3 + 1] *= tmp0;
    vtx_xyz1[iv * 3 + 2] *= tmp0;
    // add vertex
    const double tmp1 = (nf - 3.0) / (nf);
    vtx_xyz1[iv * 3 + 0] += tmp1 * vtx_xyz0[iv * 3 + 0];
    vtx_xyz1[iv * 3 + 1] += tmp1 * vtx_xyz0[iv * 3 + 1];
    vtx_xyz1[iv * 3 + 2] += tmp1 * vtx_xyz0[iv * 3 + 2];
  }
  for (unsigned int ie = 0; ie < ne0; ++ie) { // edge
    const unsigned int iv0 = aEdgeFace0[ie * 4 + 0];
    const unsigned int iv1 = aEdgeFace0[ie * 4 + 1];
    const unsigned int ie1 = aEdgeFace0[ie * 4 + 3];
    if (ie1 != UINT_MAX) { continue; }
    vtx_xyz1[iv0 * 3 + 0] = vtx_xyz0[iv0 * 3 + 0];
    vtx_xyz1[iv0 * 3 + 1] = vtx_xyz0[iv0 * 3 + 1];
    vtx_xyz1[iv0 * 3 + 2] = vtx_xyz0[iv0 * 3 + 2];
    vtx_xyz1[iv1 * 3 + 0] = vtx_xyz0[iv1 * 3 + 0];
    vtx_xyz1[iv1 * 3 + 1] = vtx_xyz0[iv1 * 3 + 1];
    vtx_xyz1[iv1 * 3 + 2] = vtx_xyz0[iv1 * 3 + 2];
  }
}

void delfem2::SubdivPoints3_MeshQuad(
    std::vector<double> &aXYZ1,
    //
    const std::vector<int> &aEdgeFace0,
    const std::vector<unsigned int> &aQuad0,
    const std::vector<double> &aXYZ0) {
  const std::size_t nv0 = aXYZ0.size() / 3;
  const std::size_t ne0 = aEdgeFace0.size() / 4;
  const std::size_t nq0 = aQuad0.size() / 4;
  assert(aEdgeFace0.size() == ne0 * 4);
  aXYZ1.resize((nv0 + ne0 + nq0) * 3);
  for (unsigned int iv = 0; iv < nv0; ++iv) {
    aXYZ1[iv * 3 + 0] = aXYZ0[iv * 3 + 0];
    aXYZ1[iv * 3 + 1] = aXYZ0[iv * 3 + 1];
    aXYZ1[iv * 3 + 2] = aXYZ0[iv * 3 + 2];
  }
  for (unsigned int ie = 0; ie < ne0; ++ie) {
    const int iv0 = aEdgeFace0[ie * 4 + 0];
    const int iv1 = aEdgeFace0[ie * 4 + 1];
    aXYZ1[(nv0 + ie) * 3 + 0] = (aXYZ0[iv0 * 3 + 0] + aXYZ0[iv1 * 3 + 0]) * 0.5;
    aXYZ1[(nv0 + ie) * 3 + 1] = (aXYZ0[iv0 * 3 + 1] + aXYZ0[iv1 * 3 + 1]) * 0.5;
    aXYZ1[(nv0 + ie) * 3 + 2] = (aXYZ0[iv0 * 3 + 2] + aXYZ0[iv1 * 3 + 2]) * 0.5;
  }
  for (unsigned int iq = 0; iq < nq0; ++iq) {
    const unsigned int iv0 = aQuad0[iq * 4 + 0];
    const unsigned int iv1 = aQuad0[iq * 4 + 1];
    const unsigned int iv2 = aQuad0[iq * 4 + 2];
    const unsigned int iv3 = aQuad0[iq * 4 + 3];
    aXYZ1[(nv0 + ne0 + iq) * 3 + 0] =
        (aXYZ0[iv0 * 3 + 0] + aXYZ0[iv1 * 3 + 0] + aXYZ0[iv2 * 3 + 0] + aXYZ0[iv3 * 3 + 0]) * 0.25;
    aXYZ1[(nv0 + ne0 + iq) * 3 + 1] =
        (aXYZ0[iv0 * 3 + 1] + aXYZ0[iv1 * 3 + 1] + aXYZ0[iv2 * 3 + 1] + aXYZ0[iv3 * 3 + 1]) * 0.25;
    aXYZ1[(nv0 + ne0 + iq) * 3 + 2] =
        (aXYZ0[iv0 * 3 + 2] + aXYZ0[iv1 * 3 + 2] + aXYZ0[iv2 * 3 + 2] + aXYZ0[iv3 * 3 + 2]) * 0.25;
  }
}

void delfem2::SubdivisionPoints_Hex(
    std::vector<double> &aXYZ1,
    //
    const std::vector<unsigned int> &psupIndHex0,
    const std::vector<unsigned int> &psupHex0,
    const std::vector<unsigned int> &aQuadHex0,
    const unsigned int *aHex0,
    unsigned int nHex0,
    const double *aXYZ0,
    unsigned int nXYZ0) {
  const unsigned int nv0 = nXYZ0;
  const std::size_t ne0 = psupHex0.size();
  const std::size_t nq0 = aQuadHex0.size() / 4;
  const unsigned int nh0 = nHex0;
  aXYZ1.resize((nv0 + ne0 + nq0 + nh0) * 3);
  for (unsigned int iv = 0; iv < nv0; ++iv) {
    aXYZ1[iv * 3 + 0] = aXYZ0[iv * 3 + 0];
    aXYZ1[iv * 3 + 1] = aXYZ0[iv * 3 + 1];
    aXYZ1[iv * 3 + 2] = aXYZ0[iv * 3 + 2];
  }
  for (unsigned int iv = 0; iv < nv0; ++iv) {
    for (unsigned int ipsup = psupIndHex0[iv]; ipsup < psupIndHex0[iv + 1]; ++ipsup) {
      unsigned int jv = psupHex0[ipsup];
      aXYZ1[(nv0 + ipsup) * 3 + 0] = (aXYZ0[iv * 3 + 0] + aXYZ0[jv * 3 + 0]) * 0.5;
      aXYZ1[(nv0 + ipsup) * 3 + 1] = (aXYZ0[iv * 3 + 1] + aXYZ0[jv * 3 + 1]) * 0.5;
      aXYZ1[(nv0 + ipsup) * 3 + 2] = (aXYZ0[iv * 3 + 2] + aXYZ0[jv * 3 + 2]) * 0.5;
    }
  }
  for (unsigned int iq = 0; iq < nq0; ++iq) {
    const unsigned int iv0 = aQuadHex0[iq * 4 + 0];
    const unsigned int iv1 = aQuadHex0[iq * 4 + 1];
    const unsigned int iv2 = aQuadHex0[iq * 4 + 2];
    const unsigned int iv3 = aQuadHex0[iq * 4 + 3];
    aXYZ1[(nv0 + ne0 + iq) * 3 + 0] =
        (aXYZ0[iv0 * 3 + 0] + aXYZ0[iv1 * 3 + 0] + aXYZ0[iv2 * 3 + 0] + aXYZ0[iv3 * 3 + 0]) * 0.25;
    aXYZ1[(nv0 + ne0 + iq) * 3 + 1] =
        (aXYZ0[iv0 * 3 + 1] + aXYZ0[iv1 * 3 + 1] + aXYZ0[iv2 * 3 + 1] + aXYZ0[iv3 * 3 + 1]) * 0.25;
    aXYZ1[(nv0 + ne0 + iq) * 3 + 2] =
        (aXYZ0[iv0 * 3 + 2] + aXYZ0[iv1 * 3 + 2] + aXYZ0[iv2 * 3 + 2] + aXYZ0[iv3 * 3 + 2]) * 0.25;
  }
  for (unsigned int ih = 0; ih < nh0; ++ih) {
    const unsigned int iv0 = aHex0[ih * 8 + 0];
    const unsigned int iv1 = aHex0[ih * 8 + 1];
    const unsigned int iv2 = aHex0[ih * 8 + 2];
    const unsigned int iv3 = aHex0[ih * 8 + 3];
    const unsigned int iv4 = aHex0[ih * 8 + 4];
    const unsigned int iv5 = aHex0[ih * 8 + 5];
    const unsigned int iv6 = aHex0[ih * 8 + 6];
    const unsigned int iv7 = aHex0[ih * 8 + 7];
    aXYZ1[(nv0 + ne0 + nq0 + ih) * 3 + 0] =
        (aXYZ0[iv0 * 3 + 0] + aXYZ0[iv1 * 3 + 0] + aXYZ0[iv2 * 3 + 0] + aXYZ0[iv3 * 3 + 0] + aXYZ0[iv4 * 3 + 0]
            + aXYZ0[iv5 * 3 + 0] + aXYZ0[iv6 * 3 + 0] + aXYZ0[iv7 * 3 + 0]) * 0.125;
    aXYZ1[(nv0 + ne0 + nq0 + ih) * 3 + 1] =
        (aXYZ0[iv0 * 3 + 1] + aXYZ0[iv1 * 3 + 1] + aXYZ0[iv2 * 3 + 1] + aXYZ0[iv3 * 3 + 1] + aXYZ0[iv4 * 3 + 1]
            + aXYZ0[iv5 * 3 + 1] + aXYZ0[iv6 * 3 + 1] + aXYZ0[iv7 * 3 + 1]) * 0.125;
    aXYZ1[(nv0 + ne0 + nq0 + ih) * 3 + 2] =
        (aXYZ0[iv0 * 3 + 2] + aXYZ0[iv1 * 3 + 2] + aXYZ0[iv2 * 3 + 2] + aXYZ0[iv3 * 3 + 2] + aXYZ0[iv4 * 3 + 2]
            + aXYZ0[iv5 * 3 + 2] + aXYZ0[iv6 * 3 + 2] + aXYZ0[iv7 * 3 + 2]) * 0.125;
  }
}
