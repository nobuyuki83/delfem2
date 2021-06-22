/*
 * Copyright (c) 2019 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */

#include <vector>
#include <cassert>
#include <stack>
#include <set>
#include <climits>

#include "delfem2/mshuni.h"

// ---------------------------------------------

DFM2_INLINE void delfem2::JArray_ElSuP_MeshElem(
    std::vector<unsigned int> &elsup_ind,
    std::vector<unsigned int> &elsup,
    // ----------
    const unsigned int* pElem,
    size_t nElem,
    unsigned int nPoEl,
    size_t nPo)
{
  //  const int nElem = (int)aElem.size()/nPoEl;
  elsup_ind.assign(nPo+1,0);
  for(unsigned int ielem=0;ielem<nElem;ielem++){
    for(unsigned int inoel=0;inoel<nPoEl;inoel++){
      const unsigned int ino1 = pElem[ielem*nPoEl+inoel];
      elsup_ind[ino1+1] += 1;
    }
  }
  for(unsigned int ino=0;ino<nPo;++ino){
    elsup_ind[ino+1] += elsup_ind[ino];
  }
  unsigned int nelsup = elsup_ind[nPo];
  elsup.resize(nelsup);
  for(unsigned int ielem=0;ielem<nElem;ielem++){
    for(unsigned int inoel=0;inoel<nPoEl;inoel++){
      int unsigned ino1 = pElem[ielem*nPoEl+inoel];
      int ind1 = elsup_ind[ino1];
      elsup[ind1] = ielem;
      elsup_ind[ino1] += 1;
    }
  }
  for(int ino=(int)nPo;ino>=1;ino--){
    elsup_ind[ino] = elsup_ind[ino-1];
  }
  elsup_ind[0] = 0;
}

DFM2_INLINE unsigned int delfem2::FindAdjEdgeIndex(
    unsigned int itri0,
    unsigned int ied0,
    unsigned int jtri0,
    const unsigned int* aTri)
{
  const unsigned int iv0 = aTri[itri0*3+(ied0+1)%3];
  const unsigned int iv1 = aTri[itri0*3+(ied0+2)%3];
  assert( iv0 != iv1 );
  assert( jtri0 != UINT_MAX );
  if( aTri[jtri0*3+1] == iv1 && aTri[jtri0*3+2] == iv0 ){ return 0; }
  if( aTri[jtri0*3+2] == iv1 && aTri[jtri0*3+0] == iv0 ){ return 1; }
  if( aTri[jtri0*3+0] == iv1 && aTri[jtri0*3+1] == iv0 ){ return 2; }
  return UINT_MAX;
}

DFM2_INLINE void delfem2::ElemQuad_DihedralTri(
    std::vector<unsigned int>& aQuad,
    const unsigned int* aTri,
    const unsigned int nTri,
    const unsigned int np)
{
  std::vector<unsigned int> aElSuEl;
  ElSuEl_MeshElem(aElSuEl,
      aTri, nTri, MESHELEM_TRI,
      np);
  assert( aElSuEl.size() == nTri*3 );
  for(unsigned int itri=0; itri<nTri; ++itri){
    for(int iedtri=0;iedtri<3;++iedtri){
      const unsigned int jtri = aElSuEl[itri*3+iedtri];
      if( jtri == UINT_MAX ) continue; // on the boundary
      if( jtri < itri ) continue;
      const unsigned int jedtri = FindAdjEdgeIndex(itri, iedtri, jtri, aTri);
      assert( jedtri != UINT_MAX );
      const unsigned int ipo0 = aTri[itri*3+iedtri];
      const unsigned int ipo1 = aTri[jtri*3+jedtri];
      const unsigned int ipo2 = aTri[itri*3+(iedtri+1)%3];
      const unsigned int ipo3 = aTri[itri*3+(iedtri+2)%3];
      assert( aTri[jtri*3+(jedtri+2)%3] == ipo2 );
      assert( aTri[jtri*3+(jedtri+1)%3] == ipo3 );
      aQuad.push_back(ipo0);
      aQuad.push_back(ipo1);
      aQuad.push_back(ipo2);
      aQuad.push_back(ipo3);
    }
  }
}

// ---------------------------------

DFM2_INLINE void delfem2::convert2Tri_Quad(
    std::vector<unsigned int>& aTri,
    const std::vector<unsigned int>& aQuad)
{
  const size_t nq = aQuad.size()/4;
  aTri.resize(nq*6);
  for(unsigned int iq=0;iq<nq;++iq){
    const unsigned int i0 = aQuad[iq*4+0];
    const unsigned int i1 = aQuad[iq*4+1];
    const unsigned int i2 = aQuad[iq*4+2];
    const unsigned int i3 = aQuad[iq*4+3];
    aTri[iq*6+0] = i0;  aTri[iq*6+1] = i1;  aTri[iq*6+2] = i2;
    aTri[iq*6+3] = i2;  aTri[iq*6+4] = i3;  aTri[iq*6+5] = i0;
  }
}


DFM2_INLINE void delfem2::FlipElement_Tri(std::vector<unsigned int>& aTri)
{
  for (std::size_t itri = 0; itri<aTri.size()/3; itri++){
    //    int i0 = aTri[itri*3+0];
    int i1 = aTri[itri*3+1];
    int i2 = aTri[itri*3+2];
    aTri[itri*3+1] = i2;
    aTri[itri*3+2] = i1;
  }
}


// -------------------------------------

DFM2_INLINE void delfem2::JArray_ElSuP_MeshTri(
    std::vector<unsigned int> &elsup_ind,
    std::vector<unsigned int> &elsup,
    // --
    const std::vector<unsigned int>& aTri,
    int nXYZ)
{
  JArray_ElSuP_MeshElem(
	  elsup_ind, elsup,
	  aTri.data(), static_cast<unsigned int>(aTri.size()/3), 3, 
	  nXYZ);
}

// ----------------------------------------------------------------------------------------------------------

DFM2_INLINE void delfem2::ElSuEl_MeshElem(
    std::vector<unsigned int>& aElSuEl,
    const unsigned int* aEl,
    size_t nEl,
    int nNoEl,
    const std::vector<unsigned int> &elsup_ind,
    const std::vector<unsigned int> &elsup,
    const int nfael,
    const int nnofa,
    const int (*noelElemFace)[4])
{
  assert( elsup_ind.size()>=1 );
  const std::size_t np = elsup_ind.size()-1;
  
  aElSuEl.assign(nEl*nfael,UINT_MAX);
  
  std::vector<int> flg_point(np,0);
  std::vector<unsigned int> inpofa(nnofa);
  for (unsigned int iel = 0; iel<nEl; iel++){
    for (int ifael=0; ifael<nfael; ifael++){
      for (int ipofa=0; ipofa<nnofa; ipofa++){
        int int0 = noelElemFace[ifael][ipofa];
        const unsigned int ip = aEl[iel*nNoEl+int0];
        assert( ip<np );
        inpofa[ipofa] = ip;
        flg_point[ip] = 1;
      }
      const int ipoin0 = inpofa[0];
      bool iflg = false;
      for (unsigned int ielsup = elsup_ind[ipoin0]; ielsup<elsup_ind[ipoin0+1]; ielsup++){
        const unsigned int jelem0 = elsup[ielsup];
        if (jelem0==iel) continue;
        for (int jfael = 0; jfael<nfael; jfael++){
          iflg = true;
          for (int jpofa = 0; jpofa<nnofa; jpofa++){
            int jnt0 = noelElemFace[jfael][jpofa];
            const unsigned int jpoin0 = aEl[jelem0*nNoEl+jnt0];
            if (flg_point[jpoin0]==0){ iflg = false; break; }
          }
          if (iflg){
            aElSuEl[iel*nfael+ifael] = jelem0;
            break;
          }
        }
        if (iflg) break;
      }
      if (!iflg){
        aElSuEl[iel*nfael+ifael] = UINT_MAX;
      }
      for (int ipofa = 0; ipofa<nnofa; ipofa++){
        flg_point[inpofa[ipofa]] = 0;
      }
    }
  }
}

/*
void makeSurroundingRelationship
(std::vector<int>& aElSurRel,
 const int* aEl, int nEl, int nNoEl,
 MESHELEM_TYPE type,
 const std::vector<int>& elsup_ind,
 const std::vector<int>& elsup)
{
  const int nfael = nFaceElem(type);
  const int nnofa = nNodeElemFace(type, 0);
  assert( nNoEl == nNodeElem(type) );
  makeSurroundingRelationship(aElSurRel,
                              aEl, nEl, nNoEl,
                              elsup_ind, elsup,
                              nfael, nnofa, noelElemFace(type));
}
 */

DFM2_INLINE void delfem2::ElSuEl_MeshElem(
    std::vector<unsigned int>& aElSuEl,
    const unsigned int* aElem,
    size_t nElem,
    MESHELEM_TYPE type,
    const size_t nXYZ)
{
  const int nNoEl = nNodeElem(type);
  std::vector<unsigned int> elsup_ind, elsup;
  JArray_ElSuP_MeshElem(elsup_ind, elsup,
      aElem, nElem, nNoEl, nXYZ);
  const int nfael = nFaceElem(type);
  const int nnofa = nNodeElemFace(type, 0);
  ElSuEl_MeshElem(aElSuEl,
      aElem, nElem, nNoEl,
      elsup_ind,elsup,
      nfael, nnofa, noelElemFace(type));
  assert( aElSuEl.size() == nElem*nfael );
}


// -------------------------------------------------------------------------

DFM2_INLINE void delfem2::JArrayPointSurPoint_MeshOneRingNeighborhood(
    std::vector<unsigned int>& psup_ind,
    std::vector<unsigned int>& psup,
    //
    const unsigned int* pElem,
    const std::vector<unsigned int> &elsup_ind,
    const std::vector<unsigned int> &elsup,
    unsigned int nnoel,
    size_t nPoint)
{
  std::vector<int> aflg(nPoint,-1);
  psup_ind.assign(nPoint+1,0);
  for(unsigned int ipoint=0;ipoint<nPoint;ipoint++){
    aflg[ipoint] = ipoint;
    for(unsigned int ielsup=elsup_ind[ipoint];ielsup<elsup_ind[ipoint+1];ielsup++){
      unsigned int jelem = elsup[ielsup];
      for(unsigned int jnoel=0;jnoel<nnoel;jnoel++){
        unsigned int jnode = pElem[jelem*nnoel+jnoel];
        if( aflg[jnode] != (int)ipoint ){
          aflg[jnode] = ipoint;
          psup_ind[ipoint+1]++;
        }
      }
    }
  }
  for(unsigned int ipoint=0;ipoint<nPoint;ipoint++){
    psup_ind[ipoint+1] += psup_ind[ipoint];
  }
  const int npsup = psup_ind[nPoint];
  psup.resize(npsup);
  for(unsigned int ipoint=0;ipoint<nPoint;ipoint++){ aflg[ipoint] = -1; }
  for(unsigned int ipoint=0;ipoint<nPoint;ipoint++){
    aflg[ipoint] = ipoint;
    for(unsigned int ielsup=elsup_ind[ipoint];ielsup<elsup_ind[ipoint+1];ielsup++){
      unsigned int jelem = elsup[ielsup];
      for(unsigned int jnoel=0;jnoel<nnoel;jnoel++){
        unsigned int jnode = pElem[jelem*nnoel+jnoel];
        if( aflg[jnode] != (int)ipoint ){
          aflg[jnode] = ipoint;
          const int ind = psup_ind[ipoint];
          psup[ind] = jnode;
          psup_ind[ipoint]++;
        }
      }
    }
  }
  for(int ipoint=(int)nPoint;ipoint>0;ipoint--){
    psup_ind[ipoint] = psup_ind[ipoint-1];
  }
  psup_ind[0] = 0;
}

DFM2_INLINE void delfem2::JArray_PSuP_MeshElem(
    std::vector<unsigned int>& psup_ind,
    std::vector<unsigned int>& psup,
    //
    const unsigned int* pElem,
    size_t nEl,
    unsigned int nPoEl,
    size_t nPo)
{
  std::vector<unsigned int> elsup_ind, elsup;
  JArray_ElSuP_MeshElem(
	  elsup_ind, elsup,
      pElem, nEl, nPoEl, nPo);
  JArrayPointSurPoint_MeshOneRingNeighborhood(
	  psup_ind, psup,
	  pElem, elsup_ind,elsup, nPoEl, nPo);
}

DFM2_INLINE void delfem2::makeOneRingNeighborhood_TriFan(
    std::vector<int>& psup_ind,
    std::vector<int>& psup,
    // ----------------------
    const std::vector<int>& aTri,
    const std::vector<int>& aTriSurRel,
    const std::vector<int>& elsup_ind,
    const std::vector<int>& elsup,
    int npoint)
{
  psup_ind.resize(npoint+1);
  psup_ind[0] = 0;
  psup.clear();
  for(int ipoint=0;ipoint<npoint;++ipoint){
    int iel0 = -1;
    int inoel0 = -1;
    {
      int ielsup0 = elsup_ind[ipoint];
      iel0 = elsup[ielsup0];
      if( aTri[iel0*3+0] == ipoint ){ inoel0 = 0; }
      if( aTri[iel0*3+1] == ipoint ){ inoel0 = 1; }
      if( aTri[iel0*3+2] == ipoint ){ inoel0 = 2; }
      assert( inoel0 != -1 );
    }
    int iel_cur = iel0;
    int inoel_cur = inoel0;
    for(;;){
      int jnoel_cur = (inoel_cur+1)%3;
      int jp0 = aTri[iel_cur*3+jnoel_cur];
      psup.push_back(jp0);
      int iel_next = aTriSurRel[iel_cur*6+2*jnoel_cur+0];
      int inoel_next = -1;
      if( aTri[iel_next*3+0] == ipoint ){ inoel_next = 0; }
      if( aTri[iel_next*3+1] == ipoint ){ inoel_next = 1; }
      if( aTri[iel_next*3+2] == ipoint ){ inoel_next = 2; }
      assert( inoel_next != -1 );
      if( iel_next == iel0 ) break;
      iel_cur = iel_next;
      inoel_cur = inoel_next;
    }
    psup_ind[ipoint+1] = (int)psup.size();
  }
}

DFM2_INLINE void delfem2::JArrayEdge_MeshElem(
    std::vector<unsigned int> &edge_ind,
    std::vector<unsigned int> &edge,
    //
    const unsigned int* aElm0,
    MESHELEM_TYPE elem_type,
    const std::vector<unsigned int> &elsup_ind,
    const std::vector<unsigned int> &elsup,
    bool is_bidirectional)
{
  const int neElm = mapMeshElemType2NEdgeElem[elem_type];
  const int nnoelElm = mapMeshElemType2NNodeElem[elem_type];
  const int (*aNoelEdge)[2] = noelElemEdge(elem_type);
  const std::size_t nPoint0 = elsup_ind.size()-1;
  edge_ind.resize(nPoint0+1);
  edge_ind[0] = 0;
  for(unsigned int ip=0;ip<nPoint0;++ip){
    std::set<int> setIP;
    for(unsigned int ielsup=elsup_ind[ip];ielsup<elsup_ind[ip+1];++ielsup){
      int iq0 = elsup[ielsup];
      for(int ie=0;ie<neElm;++ie){
        int inoel0 = aNoelEdge[ie][0];
        int inoel1 = aNoelEdge[ie][1];
        unsigned int ip0 = aElm0[iq0*nnoelElm+inoel0];
        unsigned int ip1 = aElm0[iq0*nnoelElm+inoel1];
        if( ip0 != ip && ip1 != ip ) continue;
        if( ip0 == ip ){
          if( is_bidirectional || ip1 > ip ){ setIP.insert(ip1); }
        }
        else{
          if( is_bidirectional || ip0 > ip ){ setIP.insert(ip0); }
        }
      }
    }
    for(int itr : setIP){
      edge.push_back(itr);
    }
    edge_ind[ip+1] = edge_ind[ip] + (int)setIP.size();
  }
}


DFM2_INLINE void delfem2::MeshLine_JArrayEdge(
    std::vector<unsigned int>& aLine,
    //
    const std::vector<unsigned int> &psup_ind,
    const std::vector<unsigned int> &psup)
{
  aLine.reserve(psup.size()*2);
  const std::size_t np = psup_ind.size()-1;
  for(unsigned int ip=0;ip<np;++ip){
    for(unsigned int ipsup=psup_ind[ip];ipsup<psup_ind[ip+1];++ipsup){
      unsigned int jp = psup[ipsup];
      aLine.push_back(ip);
      aLine.push_back(jp);
    }
  }
}

DFM2_INLINE void delfem2::MeshLine_MeshElem(
    std::vector<unsigned int>& aLine,
    const unsigned int* aElm0,
    unsigned int nElem,
    MESHELEM_TYPE elem_type,
    unsigned int nPo)
{
  std::vector<unsigned int> elsup_ind,elsup;
  const unsigned int nPoEl = mapMeshElemType2NNodeElem[elem_type];
  JArray_ElSuP_MeshElem(elsup_ind, elsup,
      aElm0, nElem, nPoEl, nPo);
  std::vector<unsigned int> edge_ind, edge;
  JArrayEdge_MeshElem(edge_ind, edge,
      aElm0,
      elem_type,
      elsup_ind,elsup,false);
  MeshLine_JArrayEdge(aLine,
      edge_ind,edge);
}


// -----------------------------------------

DFM2_INLINE void
delfem2::JArray_AddMasterSlavePattern(
    std::vector<unsigned int> &index,
    std::vector<unsigned int> &array,
    const unsigned int* aMSFlag,
    int ndim,
    const unsigned int *psup_ind0,
    int npsup_ind0,
    const unsigned int *psup0)
{
  assert(npsup_ind0>0);
  const int nno = npsup_ind0-1;
  //assert( aMSFlag.size() == nno*ndim );
  std::vector< std::vector<int> > mapM2S(nno);
  for(int ino1=0;ino1<nno;++ino1){
    for(int idim1=0;idim1<ndim;++idim1){
      int idof0 = aMSFlag[ino1*ndim+idim1];
      if( idof0 == -1 ){ continue; }
      int ino0 = idof0/ndim;
//      int idim0 = idof0 - ino0*ndim;
      assert( ino0 < nno && idof0 - ino0*ndim < ndim );
//      std::cout << idim1 << " " << idim0 << " " << ndim << std::endl;
      assert( idim1 == idof0 - ino0*ndim );
      mapM2S[ino0].push_back(ino1);
    }
  }
  //
  index.assign(nno+1,0);
  array.clear();
  std::vector<int> aflg(nno,-1);
  ///
  for(int ino0=0;ino0<nno;++ino0){
    aflg[ino0] = ino0;
    for(unsigned int icrs=psup_ind0[ino0];icrs<psup_ind0[ino0+1];++icrs){
      const unsigned int jno = psup0[icrs];
      if( aflg[jno] == ino0 ){ continue; }
      aflg[jno] = ino0;
      index[ino0+1]++;
    }
    for(int iino1=0;iino1<(int)mapM2S[ino0].size();++iino1){
      const int ino1 = mapM2S[ino0][iino1];
      if( aflg[ino1] != ino0 ){
        aflg[ino1] = ino0;
        index[ino0+1]++;
      }
      for(unsigned int jcrs=psup_ind0[ino1];jcrs<psup_ind0[ino1+1];++jcrs){
        const unsigned int jno1 = psup0[jcrs];
        if( aflg[jno1] == ino0 ){ continue; }
        aflg[jno1] = ino0;
        index[ino0+1]++;
      }
    }
    for(unsigned int icrs=psup_ind0[ino0];icrs<psup_ind0[ino0+1];++icrs){
      const unsigned int jno = psup0[icrs];
      for(int jdim=0;jdim<ndim;++jdim){
        int kdof = aMSFlag[jno*ndim+jdim];
        if( kdof == -1 ) continue;
        int kno = kdof/ndim;
        if( aflg[kno] == ino0 ){ continue; }
        aflg[kno] = ino0;
        index[ino0+1]++;
      }
    }
  }
  ////
  for(int ino=0;ino<nno;ino++){ index[ino+1] += index[ino]; }
  const int narray = index[nno];
  array.resize(narray);
  for(int ino=0;ino<nno;ino++){ aflg[ino] = -1; }
  ////
  for(int ino0=0;ino0<nno;++ino0){
    aflg[ino0] = ino0;
    for(unsigned int icrs=psup_ind0[ino0];icrs<psup_ind0[ino0+1];++icrs){
      const unsigned int jno = psup0[icrs];
      if( aflg[jno] == ino0 ){ continue; }
      aflg[jno] = ino0;
      const int ind = index[ino0];
      array[ind] = jno;
      index[ino0]++;
    }
    for(std::size_t jjno=0;jjno<mapM2S[ino0].size();++jjno){
      const int jno = mapM2S[ino0][jjno];
      if( aflg[jno] != ino0 ){
        aflg[jno] = ino0;
        const int ind = index[ino0];
        array[ind] = jno;
        index[ino0]++;
      }
      for(unsigned int jcrs=psup_ind0[jno];jcrs<psup_ind0[jno+1];++jcrs){
        const unsigned int kno = psup0[jcrs];
        if( aflg[kno] == ino0 ){ continue; }
        aflg[kno] = ino0;
        const int ind = index[ino0];
        array[ind] = kno;
        index[ino0]++;
      }
    }
    for(unsigned int icrs=psup_ind0[ino0];icrs<psup_ind0[ino0+1];++icrs){
      const unsigned int jno = psup0[icrs];
      for(int jdim=0;jdim<ndim;++jdim){
        int kdof = aMSFlag[jno*ndim+jdim];
        if( kdof == -1 ) continue;
        int kno = kdof/ndim;
        if( aflg[kno] == ino0 ){ continue; }
        aflg[kno] = ino0;
        const int ind = index[ino0];
        array[ind] = kno;
        index[ino0]++;
      }
    }
  }
  // ---------
  for(int ino=nno;ino>0;ino--){ index[ino] = index[ino-1]; }
  index[0] = 0;
}


// ---------------------------------------

DFM2_INLINE void delfem2::MarkConnectedElements(
    std::vector<unsigned int>& aFlagElem,
    unsigned int itri_ker,
    int igroup,
    const std::vector<unsigned int>& aElSuEl)
{
  const unsigned int nel = aFlagElem.size();
  const unsigned int nfael = aElSuEl.size()/nel;
  aFlagElem[itri_ker] = igroup;
  std::stack<int> next;
  next.push(itri_ker);
  while(!next.empty()){
    int itri0 = next.top();
    next.pop();
    for(unsigned int ie=0;ie<nfael;++ie){
      const unsigned int ita = aElSuEl[itri0*nfael+ie];
      if( ita == UINT_MAX ) continue;
      if( aFlagElem[ita] != igroup ){
        aFlagElem[ita] = igroup;
        next.push(ita);
      }
    }
  }
}

DFM2_INLINE void delfem2::MakeGroupElem(
    int& ngroup,
    std::vector<unsigned int>& aIndGroup,
    const std::vector<unsigned int>& aTri,
    const std::vector<unsigned int>& aTriSurRel,
    const int nfael,
    const int nnoel)
{
  const std::size_t nelem = aTri.size()/nnoel;
  aIndGroup.assign(nelem,-1);
  int igroup = -1;
  for(;;){
    unsigned int itri_ker = 0;
    for(;itri_ker<nelem;++itri_ker){
      if( aIndGroup[itri_ker]==-1) break;
    }
    if( itri_ker == nelem ) break;
    igroup++;
    MarkConnectedElements(aIndGroup, itri_ker, igroup, aTriSurRel);
  }
  ngroup = igroup+1;
}



