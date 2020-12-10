

#include "delfem2/mshmix.h"
#include <stack>
#include <cassert>

DFM2_INLINE void delfem2::AddElement
 (const delfem2::MESHELEM_TYPE& femelem_type,
  const std::vector<int>& aElemIn,
  //
  std::vector<unsigned int>& aElemInd,
  std::vector<unsigned int>& aElem,
  std::vector<MESHELEM_TYPE>& aElemType)
{
  const int nnoel = nNodeElem(femelem_type);
  const std::size_t nElemIn = aElemIn.size()/nnoel;
  aElemType.resize(aElemType.size()+nElemIn,femelem_type);
  //
  std::copy(aElemIn.begin(), aElemIn.end(), std::back_inserter(aElem));
  //
  assert( aElemInd.size() >= 2 );
  const std::size_t nelem0 = aElemInd.size()-1;
  const int nei0 = aElemInd[nelem0];
  aElemInd.reserve(aElemInd.size()+nElemIn);
  for(unsigned int ie=0;ie<nElemIn;++ie){
    aElemInd.push_back((ie+1)*nnoel+nei0);
  }
}


DFM2_INLINE void delfem2::JArray_ElSuP_MeshMix
 (std::vector<unsigned int> &elsup_ind,
  std::vector<unsigned int> &elsup,
  // ---
  const std::vector<unsigned int>& aElemInd,
  const std::vector<unsigned int>& aElem,
  const int nPo)
{
  const std::size_t nElem = aElemInd.size()-1;
  elsup_ind.assign(nPo+1,0);
  for(unsigned int ielem=0;ielem<nElem;ielem++){
    for(unsigned int iino=aElemInd[ielem];iino<aElemInd[ielem+1];iino++){
      int ino1 = aElem[iino];
      if( ino1 == -1 ){ break; }
      elsup_ind[ino1+1] += 1;
    }
  }
  for(int ipoint=0;ipoint<nPo;++ipoint){
    elsup_ind[ipoint+1] += elsup_ind[ipoint];
  }
  unsigned int nelsup = elsup_ind[nPo];
  elsup.resize(nelsup);
  for(unsigned int ielem=0;ielem<nElem;ielem++){
    for(unsigned int iino=aElemInd[ielem];iino<aElemInd[ielem+1];iino++){
      int ino1 = aElem[iino];
      if( ino1 == -1 ){ break; }
      int ind1 = elsup_ind[ino1];
      elsup[ind1] = ielem;
      elsup_ind[ino1] += 1;
    }
  }
  for(auto ipoint=(int)nPo;ipoint>=1;ipoint--){
    elsup_ind[ipoint] = elsup_ind[ipoint-1];
  }
  elsup_ind[0] = 0;
}


DFM2_INLINE void delfem2::ElSuEl_MeshMix
 (std::vector<int>& aElemFaceInd,
  std::vector<int>& aElemFaceRel,
  const std::vector<unsigned int>& aElemInd,
  const std::vector<unsigned int>& aElem,
  const std::vector<MESHELEM_TYPE>& aElemType,
  const int nXYZ)
{
  std::vector<unsigned int> elsup_ind, elsup;
  JArray_ElSuP_MeshMix(elsup_ind, elsup,
                       aElemInd,aElem,
                       nXYZ);
  ElSuEl_MeshMix(aElemFaceInd,aElemFaceRel,
                 aElemInd, aElem, aElemType,
                 elsup_ind,elsup);
}

DFM2_INLINE void delfem2::ElSuEl_MeshMix
 (std::vector<int>& aElemFaceInd,
  std::vector<int>& aElemFaceRel,
  //
  const std::vector<unsigned int>& aElemInd,
  const std::vector<unsigned int>& aElem,
  const std::vector<MESHELEM_TYPE> &aElemType,
  const std::vector<unsigned int> &elsup_ind,
  const std::vector<unsigned int> &elsup)
{
  assert(!aElemInd.empty());
  const std::size_t nelem = aElemInd.size()-1;
  const std::size_t np = elsup_ind.size()-1;
  assert( aElemType.size() == nelem );
  aElemFaceInd.assign(nelem+1,0);
  for(unsigned int ielem=0;ielem<nelem;++ielem){
    aElemFaceInd[ielem+1] = nFaceElem(aElemType[ielem]);
  }
  for(unsigned int ielem=0;ielem<nelem;++ielem){
    aElemFaceInd[ielem+1] += aElemFaceInd[ielem];
  }
  const int nface = aElemFaceInd[nelem];
  aElemFaceRel.assign(nface*2,-1);
  std::vector<int> aFlg(np,-1);
  for(unsigned int ielem=0;ielem<nelem;++ielem){
    const MESHELEM_TYPE type_i = aElemType[ielem];
    assert( aElemFaceInd[ielem+1]-aElemFaceInd[ielem] == nFaceElem(type_i) );
    for(int iiface=aElemFaceInd[ielem];iiface<aElemFaceInd[ielem+1];++iiface){
      const int iface = iiface-aElemFaceInd[ielem];
      const int nnofa_i = nNodeElemFace(type_i,iface);
      int ip0=-1;
      for(int inofa=0;inofa<nnofa_i;++inofa){
        const int ino0 = noelElemFace(type_i)[iface][inofa];
        assert(ino0!=-1);
        ip0 = aElem[ aElemInd[ielem]+ino0 ];
        assert(ip0>=0&&ip0<(int)np);
        aFlg[ip0] =  1;
      }
      for(unsigned int jelsup=elsup_ind[ip0];jelsup<elsup_ind[ip0+1];++jelsup){
        const int je0 = elsup[jelsup];
        if( (int)ielem == je0 ) continue;
        const MESHELEM_TYPE type_j = aElemType[je0];
        for(int ijface=aElemFaceInd[je0];ijface<aElemFaceInd[je0+1];++ijface){
          const int jface = ijface-aElemFaceInd[je0];
          const int nnofa_j = nNodeElemFace(type_j,jface);
          if( nnofa_i != nnofa_j ) continue;
          bool is_ok = true;
          for(int jnofa=0;jnofa<nnofa_j;++jnofa){
            const int jno0 = noelElemFace(type_j)[jface][jnofa];
            int jp0 = aElem[ aElemInd[je0]+jno0 ];
            if( aFlg[jp0] != 1 ){ is_ok=false; break; }
          }
          if( !is_ok ){ continue; }
          aElemFaceRel[iiface*2+0] = je0;
          aElemFaceRel[iiface*2+1] = jface;
          break;
        }
        if( aElemFaceRel[iiface*2+0] != -1 ) break;
      }
      for(int inofa=0;inofa<nnofa_i;++inofa){
        const int ino0 = noelElemFace(type_i)[iface][inofa];
        ip0 = aElem[ aElemInd[ielem]+ino0 ];
        aFlg[ip0] = -1;
      }
    }
  }
}

DFM2_INLINE void delfem2::Boundary_MeshMix
 (std::vector<unsigned int>& aElemInd_Bound,
  std::vector<unsigned int>& aElem_Bound,
  std::vector<MESHELEM_TYPE>& aElemType_Bound,
  //
  const std::vector<unsigned int>& aElemInd,
  const std::vector<unsigned int>& aElem,
  const std::vector<MESHELEM_TYPE>& aElemType,
  const std::vector<int>& aElemFaceInd,
  const std::vector<int>& aElemFaceRel)
{
  aElemType_Bound.clear();
  aElem_Bound.clear();
  aElemInd_Bound.clear();
  aElemInd_Bound.push_back(0);
  const std::size_t nelem = aElemInd.size()-1;
  for(unsigned int ielem=0;ielem<nelem;++ielem){
    const MESHELEM_TYPE type_i = aElemType[ielem];
    assert( aElemFaceInd[ielem+1]-aElemFaceInd[ielem]==nFaceElem(type_i) );
    for(int iiface=aElemFaceInd[ielem];iiface<aElemFaceInd[ielem+1];++iiface){
      if( aElemFaceRel[iiface*2+0] != -1 ) continue;
      const int iface = iiface-aElemFaceInd[ielem];
      const int nnofa_i = nNodeElemFace(type_i,iface);
      if(      nnofa_i == 3 ){ aElemType_Bound.push_back(MESHELEM_TRI ); }
      else if( nnofa_i == 4 ){ aElemType_Bound.push_back(MESHELEM_QUAD); }
      aElemInd_Bound.push_back(nnofa_i);
      for(int inofa=0;inofa<nnofa_i;++inofa){
        const int ino0 = noelElemFace(type_i)[iface][inofa];
        //          std::cout << "   " << iface << " " << inofa << " " << ino0 << std::endl;
        const int ip0 = aElem[ aElemInd[ielem]+ino0 ];
        aElem_Bound.push_back(ip0);
      }
    }
  }
  const std::size_t neb = aElemInd_Bound.size()-1;
  for(unsigned int ieb=0;ieb<neb;++ieb){
    aElemInd_Bound[ieb+1] += aElemInd_Bound[ieb];
  }
}


DFM2_INLINE void delfem2::MakeGroupElem_MeshMix
 (int& ngroup,
  std::vector<int>& aIndGroup,
  //
  const std::vector<unsigned int>& aElemInd,
  const std::vector<unsigned int>& aElem,
  const std::vector<MESHELEM_TYPE>& aElemType,
  int nPo)
{
  std::vector<unsigned int> elsup_ind, elsup;
  JArray_ElSuP_MeshMix(elsup_ind, elsup,
                       aElemInd,aElem,nPo);
  std::vector<int> aElemFaceInd, aElemFaceRel;
  ElSuEl_MeshMix(aElemFaceInd, aElemFaceRel,
                 aElemInd,aElem,aElemType,
                 elsup_ind, elsup);
  MakeGroupElem(ngroup, aIndGroup,
                aElemInd,aElem,aElemFaceInd,aElemFaceRel);
}

DFM2_INLINE void delfem2::ClipGroup_MeshMix
(std::vector<unsigned int>& aElemInd1,
 std::vector<unsigned int>& aElem1,
 std::vector<MESHELEM_TYPE>& aElemType1,
 //
 const std::vector<unsigned int>& aElemInd,
 const std::vector<unsigned int>& aElem,
 const std::vector<MESHELEM_TYPE>& aElemType,
 int igroup,
 const std::vector<int>& aIndGroup)
{
  aElem1.clear();
  aElemType1.clear();
  aElemInd1.clear();
  aElemInd1.push_back(0);
  assert( aElemInd.size() >= 2 );
  std::size_t nelem = aElemInd.size()-1;
  for(std::size_t ie=0;ie<nelem;++ie){
    if( aIndGroup[ie] != igroup ) continue;
    MESHELEM_TYPE type = aElemType[ie];
    aElemType1.push_back(type);
    aElemInd1.push_back( nNodeElem(type) );
    for(unsigned int iip=aElemInd[ie];iip<aElemInd[ie+1];++iip){
      int ip0 = aElem[iip];
      aElem1.push_back(ip0);
    }
  }
  assert( aElemInd1.size() >= 2 );
  const std::size_t ne = aElemInd1.size()-1;
  for(unsigned int ie=0;ie<ne;++ie){
    aElemInd1[ie+1] += aElemInd1[ie];
  }
}



DFM2_INLINE void delfem2::Convert2Tri_MeshMix
 (std::vector<unsigned int>& aTri,
  //
  const std::vector<unsigned int>& aElemInd,
  const std::vector<unsigned int>& aElem,
  const std::vector<MESHELEM_TYPE>& aElemType)
{
  const std::size_t nElem0 = aElemInd.size()-1;
  aTri.clear();
  aTri.reserve(nElem0*6);
  for(unsigned int ie=0;ie<nElem0;++ie){
    const int nnoel = aElemInd[ie+1]-aElemInd[ie];
    const int iip0 = aElemInd[ie];
    assert( nnoel == 3 || nnoel == 4 );
    aTri.push_back(aElem[iip0+0]);
    aTri.push_back(aElem[iip0+1]);
    aTri.push_back(aElem[iip0+2]);
    if( nnoel == 4 ){
      aTri.push_back(aElem[iip0+2]);
      aTri.push_back(aElem[iip0+3]);
      aTri.push_back(aElem[iip0+0]);
    }
  }
}

DFM2_INLINE void delfem2::FlipElement_MeshMix
 (std::vector<int>& aElem_Flip,
  // ----------
  const std::vector<unsigned int>& aElemInd,
  const std::vector<unsigned int>& aElem,
  const std::vector<MESHELEM_TYPE>& aElemType)
{
  aElem_Flip.resize(aElem.size());
  assert(!aElemInd.empty());
  const unsigned long nelem = aElemInd.size()-1;
  for(unsigned int ie=0;ie<nelem;++ie){
    const int nnoel = aElemInd[ie+1]-aElemInd[ie];
    assert( nnoel == 3 || nnoel == 4 );
    for(int inoel=0;inoel<nnoel;++inoel){
      const int ip0 = aElem[ aElemInd[ie]+inoel ];
      const int jnoel = nnoel-inoel-1;
      aElem_Flip[ aElemInd[ie]+jnoel ] = ip0;
    }
  }
}


DFM2_INLINE void delfem2::MarkConnectedElements
 (std::vector<int>& aIndGroup,
  unsigned int itri_ker,
  int igroup,
  const std::vector<int>& aElemFaceInd,
  const std::vector<int>& aElemFaceRel)
{
  assert( itri_ker < aIndGroup.size() );
  aIndGroup[itri_ker] = igroup;
  std::stack<int> next;
  next.push(itri_ker);
  while(!next.empty()){
    int ie0 = next.top();
    next.pop();
    for(int iiface=aElemFaceInd[ie0];iiface<aElemFaceInd[ie0+1];++iiface){
      assert( iiface*2 < (int)aElemFaceRel.size() );
      int je0 = aElemFaceRel[iiface*2+0];
      if( je0 == -1 ) continue;
      if( aIndGroup[je0] != igroup ){
        aIndGroup[je0] = igroup;
        next.push(je0);
      }
    }
  }
}

DFM2_INLINE void delfem2::MakeGroupElem
 (int& ngroup,
  std::vector<int>& aIndGroup,
  // -----------
  const std::vector<unsigned int>& aElemInd,
  const std::vector<unsigned int>& aElem,
  const std::vector<int>& aElemFaceInd,
  const std::vector<int>& aElemFaceRel)
{
  assert(!aElemInd.empty());
  const int unsigned nelem = aElemInd.size()-1;
  assert( aElemFaceInd.size() == (nelem+1) );
  aIndGroup.assign(nelem,-1);
  int igroup = -1;
  for(;;){
    unsigned int ielem_ker = 0;
    for(;ielem_ker<nelem;++ielem_ker){
      if( aIndGroup[ielem_ker]==-1) break;
    }
    if( ielem_ker == nelem ) break;
    igroup++;
    MarkConnectedElements(aIndGroup,
                          ielem_ker, igroup,
                          aElemFaceInd,aElemFaceRel);
  }
  ngroup = igroup+1;
}


