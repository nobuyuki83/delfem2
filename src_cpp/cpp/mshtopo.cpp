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
#include <iostream>

#include "delfem2/mshtopo.h"

void JArray_Print
(const std::vector<int>& index,
 const std::vector<int>& array)
{
  const int np = index.size()-1;
  for(int ip=0;ip<np;++ip){
    std::cout << ip << " --> ";
    for(int ipsup=index[ip];ipsup<index[ip+1];++ipsup){
      std::cout << array[ipsup] << " ";
    }
    std::cout << std::endl;
  }
}

void JArray_Sort
(const std::vector<int>& index,
 std::vector<int>& array)
{
  if( index.size() == 0 ) return;
  const int size = (int)index.size()-1;
  for(int ipoin=0;ipoin<size;ipoin++){
    const int is = index[ipoin  ];
    const int ie = index[ipoin+1];
    if( is == ie ) continue;
    assert( is < ie );
    int itmp;
    for(int i=is;i<ie-1;i++){
      for(int j=ie-1;j>i;j--){
        if( array[j] < array[j-1] ){
          itmp = array[j];
          array[j] = array[j-1];
          array[j-1] = itmp;
        }
      }
    }
  }
}

void JArray_Sort
(const int* index, const int size,
 int* array)
{
  if( size == 0 ) return;
//  if( index.size() == 0 ) return;
//  const int size = (int)index.size()-1;
  for(int ipoin=0;ipoin<size;ipoin++){
    const int is = index[ipoin  ];
    const int ie = index[ipoin+1];
    if( is == ie ) continue;
    assert( is < ie );
    int itmp;
    for(int i=is;i<ie-1;i++){
      for(int j=ie-1;j>i;j--){
        if( array[j] < array[j-1] ){
          itmp = array[j];
          array[j] = array[j-1];
          array[j-1] = itmp;
        }
      }
    }
  }
}

void JArray_AddDiagonal
(std::vector<int >& psup_ind1,
 std::vector<int >& psup1,
 const int* psup_ind0, int npsup_ind0,
 const int* psup0, int npsup0)
{
  const int np = npsup_ind0-1;
  std::vector<int> tmp(np,-1);
  psup_ind1.assign(np+1,0);
  for(int ip=0;ip<np;++ip){
    for(int ipsup=psup_ind0[ip];ipsup<psup_ind0[ip+1];++ipsup){
      const int jp = psup0[ipsup];
      assert( tmp[jp] != ip );
      tmp[jp] = ip;
      psup_ind1[ip+1] += 1;
    }
    if( tmp[ip] != ip ){
      tmp[ip] = ip;
      psup_ind1[ip+1] += 1;
    }
  }
  /////
  for(int ip=0;ip<np;++ip){
    psup_ind1[ip+1] += psup_ind1[ip];
  }
  const int npsup = psup_ind1[np];
  psup1.resize(npsup);
  tmp.assign(np,-1);
  /////
  for(int ip=0;ip<np;++ip){
    for(int ipsup=psup_ind0[ip];ipsup<psup_ind0[ip+1];++ipsup){
      const int jp = psup0[ipsup];
      assert( tmp[jp] != ip );
      tmp[jp] = ip;
      int iclstr  = psup_ind1[ip];
      psup1[ iclstr ] = jp;
      psup_ind1[ip] += 1;
    }
    if( tmp[ip] != ip ){
      int iclstr  = psup_ind1[ip];
      psup1[ iclstr ] = ip;
      psup_ind1[ip] += 1;
    }
  }
  //////
  for(int ip=np-1;ip>=0;--ip){
    psup_ind1[ip+1] = psup_ind1[ip];
  }
  psup_ind1[0] = 0;
}

// in the edge ip -> jp, it holds (ip < jp)
void JArrayEdgeUnidir_PointSurPoint
(std::vector<int>& edge_ind,
 std::vector<int>& edge,
 /////
 const std::vector<int>& psup_ind,
 const std::vector<int>& psup)
{
  const int np = psup_ind.size()-1;
  edge_ind.resize(np+1);
  edge_ind[0] = 0;
  edge.clear();
  ////
  for(int ip=0;ip<np;++ip){
    for(int ipsup=psup_ind[ip];ipsup<psup_ind[ip+1];++ipsup){
      int ip0 = psup[ipsup];
      if( ip0 <= ip ) continue;
      edge_ind[ip+1]++;
    }
  }
  for(int ip=0;ip<np;ip++){
    edge_ind[ip+1] += edge_ind[ip];
  }
  const int nedge = edge_ind[np];
  edge.resize(nedge);
  for(int ip=0;ip<np;++ip){
    for(int ipsup=psup_ind[ip];ipsup<psup_ind[ip+1];++ipsup){
      const int ip0 = psup[ipsup];
      if( ip0 <= ip ) continue;
      const int iedge = edge_ind[ip];
      edge[iedge] = ip0;
      edge_ind[ip]++;
    }
  }
  for(int ip=np;ip>0;ip--){
    edge_ind[ip] = edge_ind[ip-1];
  }
  edge_ind[0] = 0;
}

void JArrayElemSurPoint_MeshElem
(std::vector<int>& elsup_ind,
 std::vector<int>& elsup,
 ////
 const unsigned int* pElem,
 int nElem,
 int nPoEl,
 int nPo)
{
  //  const int nElem = (int)aElem.size()/nPoEl;
  elsup_ind.assign(nPo+1,0);
  for(int ielem=0;ielem<nElem;ielem++){
    for(int inoel=0;inoel<nPoEl;inoel++){
      const int ino1 = pElem[ielem*nPoEl+inoel];
      if( ino1 == -1 ){ break; }
      elsup_ind[ino1+1] += 1;
    }
  }
  for(int ino=0;ino<nPo;++ino){
    elsup_ind[ino+1] += elsup_ind[ino];
  }
  int nelsup = elsup_ind[nPo];
  elsup.resize(nelsup);
  for(int ielem=0;ielem<nElem;ielem++){
    for(int inoel=0;inoel<nPoEl;inoel++){
      int ino1 = pElem[ielem*nPoEl+inoel];
      if( ino1 == -1 ){ break; }
      int ind1 = elsup_ind[ino1];
      elsup[ind1] = ielem;
      elsup_ind[ino1] += 1;
    }
  }
  for(int ino=nPo;ino>=1;ino--){
    elsup_ind[ino] = elsup_ind[ino-1];
  }
  elsup_ind[0] = 0;
}

///////////////////////////////////////////////////////////////////////////////

void ElemQuad_DihedralTri
(std::vector<unsigned int>& aQuad,
 const unsigned int* aTri, int nTri,
 int np)
{
  std::vector<int> aElemSurRel;
  makeSurroundingRelationship(aElemSurRel,
                              aTri, nTri,
                              MESHELEM_TRI, np);
  ////
  for(int itri=0; itri<nTri; ++itri){
    for(int iedtri=0;iedtri<3;++iedtri){
      int jtri = aElemSurRel[itri*6+iedtri*2+0];
      if( jtri == -1 ) continue;
      if( jtri < itri ) continue;
      int jedtri = aElemSurRel[itri*6+iedtri*2+1];
      assert( itri == aElemSurRel[jtri*6+jedtri*2+0] );
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

///////////////////////////////////////////////////////////////////////////////

void convert2Tri_Quad
(std::vector<unsigned int>& aTri,
 const std::vector<unsigned int>& aQuad)
{
  const unsigned long nq = aQuad.size()/4;
  aTri.resize(nq*6);
  for(unsigned int iq=0;iq<nq;++iq){
    const int i0 = aQuad[iq*4+0];
    const int i1 = aQuad[iq*4+1];
    const int i2 = aQuad[iq*4+2];
    const int i3 = aQuad[iq*4+3];
    aTri[iq*6+0] = i0;  aTri[iq*6+1] = i1;  aTri[iq*6+2] = i2;
    aTri[iq*6+3] = i2;  aTri[iq*6+4] = i3;  aTri[iq*6+5] = i0;
  }
}

void convert2Tri
(std::vector<int>& aTri,
 ////
 const std::vector<int>& aElemInd,
 const std::vector<int>& aElem,
 const std::vector<MESHELEM_TYPE>& aElemType)
{
  const long nElem0 = aElemInd.size()-1;
  aTri.clear();
  aTri.reserve(nElem0*6);
  for(int ie=0;ie<nElem0;++ie){
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

void FlipElement_Tri(std::vector<int>& aTri)
{
  for (unsigned int itri = 0; itri<aTri.size()/3; itri++){
    //    int i0 = aTri[itri*3+0];
    int i1 = aTri[itri*3+1];
    int i2 = aTri[itri*3+2];
    aTri[itri*3+1] = i2;
    aTri[itri*3+2] = i1;
  }
}

void FlipElement
(std::vector<int>& aElem_Flip,
 ////
 const std::vector<int>& aElemInd,
 const std::vector<int>& aElem,
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


/////////////////////////////////////////////////////////////////////////

void AddElement
(const MESHELEM_TYPE& femelem_type,
 const std::vector<int>& aElemIn,
 ////
 std::vector<int>& aElemInd,
 std::vector<int>& aElem,
 std::vector<MESHELEM_TYPE>& aElemType)
{
  const int nnoel = nNodeElem(femelem_type);
  const int nElemIn = aElemIn.size()/nnoel;
  aElemType.resize(aElemType.size()+nElemIn,femelem_type);
  ////
  std::copy(aElemIn.begin(), aElemIn.end(), std::back_inserter(aElem));
  ////
  const int nelem0 = aElemInd.size()-1;
  const int nei0 = aElemInd[nelem0];
  aElemInd.reserve(aElemInd.size()+nElemIn);
  for(int ie=0;ie<nElemIn;++ie){ aElemInd.push_back((ie+1)*nnoel+nei0); }
}

/////////////////////////////////////////////////////////////////////////////////////



void JArrayElemSurPoint_MeshTri
(std::vector<int>& elsup_ind,
 std::vector<int>& elsup,
 ////
 const std::vector<unsigned int>& aTri,
 int nXYZ)
{
  JArrayElemSurPoint_MeshElem(elsup_ind, elsup,
                           aTri.data(), aTri.size()/3, 3, nXYZ);
}

void JArrayElemSurPoint_MeshMix
(std::vector<int>& elsup_ind,
 std::vector<int>& elsup,
 ////
 const std::vector<int>& aElemInd,
 const std::vector<int>& aElem,
 const int nPo)
{
  const int nElem = aElemInd.size()-1;
  elsup_ind.assign(nPo+1,0);
  for(int ielem=0;ielem<nElem;ielem++){
    for(int iino=aElemInd[ielem];iino<aElemInd[ielem+1];iino++){
      int ino1 = aElem[iino];
      if( ino1 == -1 ){ break; }
      elsup_ind[ino1+1] += 1;
    }
  }
  for(int ipoint=0;ipoint<nPo;++ipoint){
    elsup_ind[ipoint+1] += elsup_ind[ipoint];
  }
  int nelsup = elsup_ind[nPo];
  elsup.resize(nelsup);
  for(int ielem=0;ielem<nElem;ielem++){
    for(int iino=aElemInd[ielem];iino<aElemInd[ielem+1];iino++){
      int ino1 = aElem[iino];
      if( ino1 == -1 ){ break; }
      int ind1 = elsup_ind[ino1];
      elsup[ind1] = ielem;
      elsup_ind[ino1] += 1;
    }
  }
  for(int ipoint=nPo;ipoint>=1;ipoint--){
    elsup_ind[ipoint] = elsup_ind[ipoint-1];
  }
  elsup_ind[0] = 0;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void makeSurroundingRelationship
(std::vector<int>& aElSurRel,
 const unsigned int* aEl, int nEl, int nNoEl,
 const std::vector<int>& elsup_ind,
 const std::vector<int>& elsup,
 const int nfael,
 const int nnofa,
 const int noelElemFace[][4])
{
//  std::cout << nfael << " " << nnofa << " " << nnoel << std::endl;
  const int np = (int)elsup_ind.size()-1;
  
  aElSurRel.assign(nEl*nfael*2,-1);
  
  std::vector<int> tmp_poin(np,0);
  std::vector<int> inpofa(nnofa);
  for (int iel = 0; iel<nEl; iel++){
    for (int ifael=0; ifael<nfael; ifael++){
      for (int ipofa=0; ipofa<nnofa; ipofa++){
        int int0 = noelElemFace[ifael][ipofa];
        const int ip = aEl[iel*nNoEl+int0];
        assert( ip>=0 && ip<np );
        inpofa[ipofa] = ip;
        tmp_poin[ip] = 1;
      }
      const int ipoin0 = inpofa[0];
      bool iflg = false;
      for (int ielsup = elsup_ind[ipoin0]; ielsup<elsup_ind[ipoin0+1]; ielsup++){
        const int jelem0 = elsup[ielsup];
        if (jelem0==iel) continue;
        for (int jfael = 0; jfael<nfael; jfael++){
          iflg = true;
          for (int jpofa = 0; jpofa<nnofa; jpofa++){
            int jnt0 = noelElemFace[jfael][jpofa];
            const int jpoin0 = aEl[jelem0*nNoEl+jnt0];
            if (tmp_poin[jpoin0]==0){ iflg = false; break; }
          }
          if (iflg){
            aElSurRel[iel*nfael*2+ifael*2+0] = jelem0;
            aElSurRel[iel*nfael*2+ifael*2+1] = jfael;
            break;
          }
        }
        if (iflg) break;
      }
      if (!iflg){
        aElSurRel[iel*nfael*2+ifael*2+0] = -1;
        aElSurRel[iel*nfael*2+ifael*2+1] = -1;
      }
      for (int ipofa = 0; ipofa<nnofa; ipofa++){
        tmp_poin[inpofa[ipofa]] = 0;
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

void makeSurroundingRelationship
(std::vector<int>& aElemSurRel,
 const unsigned int* aElem, int nElem,
 MESHELEM_TYPE type,
 const int nXYZ)
{
  const int nNoEl = nNodeElem(type);
  std::vector<int> elsup_ind, elsup;
  JArrayElemSurPoint_MeshElem(elsup_ind, elsup,
                           aElem, nElem, nNoEl, nXYZ);
  const int nfael = nFaceElem(type);
  const int nnofa = nNodeElemFace(type, 0);
  makeSurroundingRelationship(aElemSurRel,
                              aElem, nElem, nNoEl,
                              elsup_ind,elsup,
                              nfael, nnofa, noelElemFace(type));
}


void makeSurroundingRelationship
(std::vector<int>& aElemFaceInd,
 std::vector<int>& aElemFaceRel,
 const std::vector<int>& aElemInd,
 const std::vector<int>& aElem,
 const std::vector<MESHELEM_TYPE>& aElemType,
 const int nXYZ)
{
  std::vector<int> elsup_ind, elsup;
  JArrayElemSurPoint_MeshMix(elsup_ind, elsup,
                           aElemInd,aElem,
                           nXYZ);
  makeSurroundingRelationship(aElemFaceInd,aElemFaceRel,
                              aElemInd, aElem, aElemType,
                              elsup_ind,elsup);
}

void makeSurroundingRelationship
(std::vector<int>& aElemFaceInd,
 std::vector<int>& aElemFaceRel,
 ///
 const std::vector<int>& aElemInd,
 const std::vector<int>& aElem,
 const std::vector<MESHELEM_TYPE>& aElemType,
 const std::vector<int>& elsup_ind,
 const std::vector<int>& elsup)
{
  assert(aElemInd.size()>0);
  const unsigned int nelem = aElemInd.size()-1;
  const int np = elsup_ind.size()-1;
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
        assert(ip0>=0&&ip0<np);
        aFlg[ip0] =  1;
      }
      for(int jelsup=elsup_ind[ip0];jelsup<elsup_ind[ip0+1];++jelsup){
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

void makeBoundary
(std::vector<int>& aElemInd_Bound,
 std::vector<int>& aElem_Bound,
 std::vector<MESHELEM_TYPE>& aElemType_Bound,
 ////
 const std::vector<int>& aElemInd,
 const std::vector<int>& aElem,
 const std::vector<MESHELEM_TYPE>& aElemType,
 const std::vector<int>& aElemFaceInd,
 const std::vector<int>& aElemFaceRel)
{
  aElemType_Bound.clear();
  aElem_Bound.clear();
  aElemInd_Bound.clear();
  aElemInd_Bound.push_back(0);
  const int nelem = aElemInd.size()-1;
  for(int ielem=0;ielem<nelem;++ielem){
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
  const int neb = aElemInd_Bound.size()-1;
  for(int ieb=0;ieb<neb;++ieb){
    aElemInd_Bound[ieb+1] += aElemInd_Bound[ieb];
  }
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void JArrayPointSurPoint_MeshOneRingNeighborhood
(std::vector<int>& psup_ind,
 std::vector<int>& psup,
 ////
 const unsigned int* pElem,
 const std::vector<int>& elsup_ind,
 const std::vector<int>& elsup,
 int nnoel,
 int nPoint)
{
  std::vector<int> aflg(nPoint,-1);
  psup_ind.assign(nPoint+1,0);
  for(int ipoint=0;ipoint<nPoint;ipoint++){
    aflg[ipoint] = ipoint;
    for(int ielsup=elsup_ind[ipoint];ielsup<elsup_ind[ipoint+1];ielsup++){
      int jelem = elsup[ielsup];
      for(int jnoel=0;jnoel<nnoel;jnoel++){
        int jnode = pElem[jelem*nnoel+jnoel];
        if( aflg[jnode] != ipoint ){
          aflg[jnode] = ipoint;
          psup_ind[ipoint+1]++;
        }
      }
    }
  }
  for(int ipoint=0;ipoint<nPoint;ipoint++){
    psup_ind[ipoint+1] += psup_ind[ipoint];
  }
  const int npsup = psup_ind[nPoint];
  psup.resize(npsup);
  for(int ipoint=0;ipoint<nPoint;ipoint++){ aflg[ipoint] = -1; }
  for(int ipoint=0;ipoint<nPoint;ipoint++){
    aflg[ipoint] = ipoint;
    for(int ielsup=elsup_ind[ipoint];ielsup<elsup_ind[ipoint+1];ielsup++){
      int jelem = elsup[ielsup];
      for(int jnoel=0;jnoel<nnoel;jnoel++){
        int jnode = pElem[jelem*nnoel+jnoel];
        if( aflg[jnode] != ipoint ){
          aflg[jnode] = ipoint;
          const int ind = psup_ind[ipoint];
          psup[ind] = jnode;
          psup_ind[ipoint]++;
        }
      }
    }
  }
  for(int ipoint=nPoint;ipoint>0;ipoint--){
    psup_ind[ipoint] = psup_ind[ipoint-1];
  }
  psup_ind[0] = 0;
}

void JArrayPointSurPoint_MeshOneRingNeighborhood
(std::vector<int>& psup_ind,
 std::vector<int>& psup,
 ////
 const unsigned int* pElem,
 int nEl,
 int nPoEl,
 int nPo)
{
  std::vector<int> elsup_ind, elsup;
  JArrayElemSurPoint_MeshElem(elsup_ind, elsup,
                           pElem, nEl, nPoEl, nPo);
  JArrayPointSurPoint_MeshOneRingNeighborhood(psup_ind, psup,
                          pElem, elsup_ind,elsup, nPoEl, nPo);
}

void makeOneRingNeighborhood_TriFan
(std::vector<int>& psup_ind,
 std::vector<int>& psup,
 ////
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

void JArrayEdge_MeshElem
(std::vector<int>& edge_ind,
 std::vector<int>& edge,
 ////
 const unsigned int* aElm0,
 MESHELEM_TYPE elem_type,
 const std::vector<int>& elsup_ind,
 const std::vector<int>& elsup,
 bool is_bidirectional)
{
  const int neElm = mapMeshElemType2NEdgeElem[elem_type];
  const int nnoelElm = mapMeshElemType2NNodeElem[elem_type];
  const int (*aNoelEdge)[2] = noelElemEdge(elem_type);
  const int nPoint0 = elsup_ind.size()-1;
  edge_ind.resize(nPoint0+1);
  edge_ind[0] = 0;
  for(int ip=0;ip<nPoint0;++ip){
    std::set<int> setIP;
    for(int ielsup=elsup_ind[ip];ielsup<elsup_ind[ip+1];++ielsup){
      int iq0 = elsup[ielsup];
      for(int ie=0;ie<neElm;++ie){
        int inoel0 = aNoelEdge[ie][0];
        int inoel1 = aNoelEdge[ie][1];
        int ip0 = aElm0[iq0*nnoelElm+inoel0];
        int ip1 = aElm0[iq0*nnoelElm+inoel1];
        if( ip0 != ip && ip1 != ip ) continue;
        if( ip0 == ip ){
          if( is_bidirectional || ip1 > ip ){ setIP.insert(ip1); }
        }
        else{
          if( is_bidirectional || ip0 > ip ){ setIP.insert(ip0); }
        }
      }
    }
    for(std::set<int>::iterator itr = setIP.begin();itr!=setIP.end();++itr){
      edge.push_back(*itr);
    }
    edge_ind[ip+1] = edge_ind[ip] + (int)setIP.size();
  }
}


void MeshLine_JArrayEdge
(std::vector<unsigned int>& aLine,
 ////
 const std::vector<int>& psup_ind,
 const std::vector<int>& psup)
{
  aLine.reserve(psup.size()*2);
  const int np = psup_ind.size()-1;
  for(int ip=0;ip<np;++ip){
    for(int ipsup=psup_ind[ip];ipsup<psup_ind[ip+1];++ipsup){
      int jp = psup[ipsup];
      aLine.push_back(ip);
      aLine.push_back(jp);
    }
  }
}

void MeshLine_MeshElem
(std::vector<unsigned int>& aLine,
 const unsigned int* aElm0,
 unsigned int nElem,
 MESHELEM_TYPE elem_type,
 unsigned int nPo)
{
  std::vector<int> elsup_ind,elsup;
  const unsigned int nPoEl = mapMeshElemType2NNodeElem[elem_type];
  JArrayElemSurPoint_MeshElem(elsup_ind, elsup,
                           aElm0, nElem, nPoEl, nPo);
  std::vector<int> edge_ind, edge;
  JArrayEdge_MeshElem(edge_ind, edge,
                      aElm0,
                      elem_type,
                      elsup_ind,elsup,false);
  MeshLine_JArrayEdge(aLine,
                      edge_ind,edge);
}


//////////////////////////////////////////////////////////////////////////////////////////////////

void JArray_AddMasterSlavePattern
(std::vector<int>& index,
 std::vector<int>& array,
 const int* aMSFlag,
 int ndim,
 const int* psup_ind0,
 int npsup_ind0,
 const int* psup0)
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
      int idim0 = idof0 - ino0*ndim;
      assert( ino0 < nno && idim0 < ndim );
//      std::cout << idim1 << " " << idim0 << " " << ndim << std::endl;
      assert( idim1 == idim0 );
      mapM2S[ino0].push_back(ino1);
    }
  }
  ////
  index.assign(nno+1,0);
  array.clear();
  std::vector<int> aflg(nno,-1);
  /////
  for(int ino0=0;ino0<nno;++ino0){
    aflg[ino0] = ino0;
    for(int icrs=psup_ind0[ino0];icrs<psup_ind0[ino0+1];++icrs){
      const int jno = psup0[icrs];
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
      for(int jcrs=psup_ind0[ino1];jcrs<psup_ind0[ino1+1];++jcrs){
        const int jno1 = psup0[jcrs];
        if( aflg[jno1] == ino0 ){ continue; }
        aflg[jno1] = ino0;
        index[ino0+1]++;
      }
    }
    for(int icrs=psup_ind0[ino0];icrs<psup_ind0[ino0+1];++icrs){
      const int jno = psup0[icrs];
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
    for(int icrs=psup_ind0[ino0];icrs<psup_ind0[ino0+1];++icrs){
      const int jno = psup0[icrs];
      if( aflg[jno] == ino0 ){ continue; }
      aflg[jno] = ino0;
      const int ind = index[ino0];
      array[ind] = jno;
      index[ino0]++;
    }
    for(unsigned int jjno=0;jjno<mapM2S[ino0].size();++jjno){
      const int jno = mapM2S[ino0][jjno];
      if( aflg[jno] != ino0 ){
        aflg[jno] = ino0;
        const int ind = index[ino0];
        array[ind] = jno;
        index[ino0]++;
      }
      for(int jcrs=psup_ind0[jno];jcrs<psup_ind0[jno+1];++jcrs){
        const int kno = psup0[jcrs];
        if( aflg[kno] == ino0 ){ continue; }
        aflg[kno] = ino0;
        const int ind = index[ino0];
        array[ind] = kno;
        index[ino0]++;
      }
    }
    for(int icrs=psup_ind0[ino0];icrs<psup_ind0[ino0+1];++icrs){
      const int jno = psup0[icrs];
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
  /////
  for(int ino=nno;ino>0;ino--){ index[ino] = index[ino-1]; }
  index[0] = 0;
}


/////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////

void MarkConnectedElements
(std::vector<int>& aIndGroup,
 int itri_ker,
 int igroup,
 const std::vector<int>& aTriSurRel,
 const int nfael)
{
  aIndGroup[itri_ker] = igroup;
  std::stack<int> next;
  next.push(itri_ker);
  while(!next.empty()){
    int itri0 = next.top();
    next.pop();
    for(int ie=0;ie<nfael;++ie){
      const int ita = aTriSurRel[(itri0*nfael+ie)*2+0];
      if( ita == -1 ) continue;
      if( aIndGroup[ita] != igroup ){
        aIndGroup[ita] = igroup;
        next.push(ita);
      }
    }
  }
}

void MarkConnectedElements
(std::vector<int>& aIndGroup,
 int itri_ker,
 int igroup,
 const std::vector<int>& aElemFaceInd,
 const std::vector<int>& aElemFaceRel)
{
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

void MakeGroupElem
(int& ngroup,
 std::vector<int>& aIndGroup,
 const std::vector<int>& aTri,
 const std::vector<int>& aTriSurRel,
 const int nfael,
 const int nnoel)
{
  ////
  const int nelem = aTri.size()/nnoel;
  aIndGroup.assign(nelem,-1);
  int igroup = -1;
  for(;;){
    int itri_ker = 0;
    for(;itri_ker<nelem;++itri_ker){
      if( aIndGroup[itri_ker]==-1) break;
    }
    if( itri_ker == nelem ) break;
    igroup++;
    MarkConnectedElements(aIndGroup, itri_ker, igroup, aTriSurRel,nfael);
  }
  ngroup = igroup+1;
}

void MakeGroupElem_Tri
(int& ngroup,
 std::vector<int>& aIndGroup,
 const std::vector<int>& aTri,
 const std::vector<int>& aTriSurRel)
{
  MakeGroupElem(ngroup,aIndGroup,
                aTri,aTriSurRel,3,3);
}

void MakeGroupElem
(int& ngroup,
 std::vector<int>& aIndGroup,
 ////
 const std::vector<int>& aElemInd,
 const std::vector<int>& aElem,
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


void MakeGroupElem
(int& ngroup,
 std::vector<int>& aIndGroup,
 /////
 const std::vector<int>& aElemInd,
 const std::vector<int>& aElem,
 const std::vector<MESHELEM_TYPE>& aElemType,
 int nPo)
{
  std::vector<int> elsup_ind, elsup;
  JArrayElemSurPoint_MeshMix(elsup_ind, elsup,
                           aElemInd,aElem,nPo);
  std::vector<int> aElemFaceInd, aElemFaceRel;
  makeSurroundingRelationship(aElemFaceInd, aElemFaceRel,
                              aElemInd,aElem,aElemType,
                              elsup_ind, elsup);
  MakeGroupElem(ngroup, aIndGroup,
                aElemInd,aElem,aElemFaceInd,aElemFaceRel);
}

void ClipGroup(std::vector<int>& aElemInd1,
               std::vector<int>& aElem1,
               std::vector<MESHELEM_TYPE>& aElemType1,
               ///
               const std::vector<int>& aElemInd,
               const std::vector<int>& aElem,
               const std::vector<MESHELEM_TYPE>& aElemType,
               int igroup,
               const std::vector<int>& aIndGroup)
{
  aElem1.clear();
  aElemType1.clear();
  aElemInd1.clear();
  aElemInd1.push_back(0);
  int nelem = aElemInd.size()-1;
  for(int ie=0;ie<nelem;++ie){
    if( aIndGroup[ie] != igroup ) continue;
    MESHELEM_TYPE type = aElemType[ie];
    aElemType1.push_back(type);
    aElemInd1.push_back( nNodeElem(type) );
    for(int iip=aElemInd[ie];iip<aElemInd[ie+1];++iip){
      int ip0 = aElem[iip];
      aElem1.push_back(ip0);
    }
  }
  const int ne = aElemInd1.size()-1;
  for(int ie=0;ie<ne;++ie){
    aElemInd1[ie+1] += aElemInd1[ie];
  }
}

////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////

int findEdge
(int ip0, int ip1,
 const std::vector<int>& psup_ind,
 const std::vector<int>& psup)
{
  if( ip1 > ip0 ){
    for(int ipsup=psup_ind[ip0];ipsup<psup_ind[ip0+1];++ipsup){
      int ip2 = psup[ipsup];
      if( ip2 == ip1 ){ return ipsup; }
    }
  }
  else{
    for(int ipsup=psup_ind[ip1];ipsup<psup_ind[ip1+1];++ipsup){
      int ip2 = psup[ipsup];
      if( ip2 == ip0 ){ return ipsup; }
    }
  }
  return -1;
}

int findFace
(int ip0, int ip1, int ip2, int ip3,
 const std::vector<unsigned int>& aQuad,
 const std::vector<int>& elsupInd,
 const std::vector<int>& elsup)
{
  if( ip0 < 0 || ip0 >= (int)elsupInd.size()-1 ) return -1;
  assert( ip0 >=0 && ip0 < (int)elsupInd.size()-1 );
  for(int ielsup=elsupInd[ip0];ielsup<elsupInd[ip0+1];++ielsup){
    int ie0 = elsup[ielsup];
    int iq0 = aQuad[ie0*4+0];
    int iq1 = aQuad[ie0*4+1];
    int iq2 = aQuad[ie0*4+2];
    int iq3 = aQuad[ie0*4+3];
    if( ip0!=iq0 && ip0!=iq1 && ip0!=iq2 && ip0!=iq3 ) continue;
    if( ip1!=iq0 && ip1!=iq1 && ip1!=iq2 && ip1!=iq3 ) continue;
    if( ip2!=iq0 && ip2!=iq1 && ip2!=iq2 && ip2!=iq3 ) continue;
    if( ip3!=iq0 && ip3!=iq1 && ip3!=iq2 && ip3!=iq3 ) continue;
    return ie0;
  }
  return -1;
}

// new points is in the order of [old points], [edge points], [face points]
void QuadSubdiv
(std::vector<unsigned int>& aQuad1,
 std::vector<int>& psup_ind,
 std::vector<int>& psup,
 std::vector<int>& aEdgeFace0, // two points on the edge and two quads touching the edge
 const unsigned int* aQuad0, int nQuad0,
 unsigned int nPoint0)
{
  const int nq0 = nQuad0;
  std::vector<int> elsup_ind, elsup;
  JArrayElemSurPoint_MeshElem(elsup_ind,elsup,
                           aQuad0,nQuad0,4,nPoint0);
  JArrayEdge_MeshElem(psup_ind,psup,
                       aQuad0, MESHELEM_QUAD, elsup_ind, elsup,
                       false); // is_bidirectional = false
  const unsigned int ne0 = (int)psup.size();
  aEdgeFace0.resize(0);
  aEdgeFace0.reserve(ne0*4);
  for(int ip=0;ip<(int)nPoint0;++ip){
    for(int ipsup=psup_ind[ip];ipsup<psup_ind[ip+1];++ipsup){
      int ip1 = psup[ipsup];
      aEdgeFace0.push_back(ip);
      aEdgeFace0.push_back(ip1);
      int iq0=-1, iq1=-1;
      for(int ielsup=elsup_ind[ip];ielsup<elsup_ind[ip+1];++ielsup){
        int jq0 = elsup[ielsup];
        int jp0 = aQuad0[jq0*4+0];
        int jp1 = aQuad0[jq0*4+1];
        int jp2 = aQuad0[jq0*4+2];
        int jp3 = aQuad0[jq0*4+3];
        if( (jp0!=ip) && (jp1!=ip) && (jp2!=ip) && (jp3!=ip) ){ continue; }
        if( (jp0!=ip1) && (jp1!=ip1) && (jp2!=ip1) && (jp3!=ip1) ){ continue; }
        //////
        if( iq0 == -1 ){ iq0 = jq0; }
        else{
          assert( iq1 == -1 );
          iq1 = jq0;
        }
      }
      aEdgeFace0.push_back(iq0);
      aEdgeFace0.push_back(iq1);
    }
  }
  aQuad1.resize(0);
  aQuad1.reserve(nQuad0*4);
  for(int iq=0;iq<nq0;++iq){
    int ip0 = aQuad0[iq*4+0];
    int ip1 = aQuad0[iq*4+1];
    int ip2 = aQuad0[iq*4+2];
    int ip3 = aQuad0[iq*4+3];
    int ie01 = findEdge(ip0,ip1, psup_ind,psup); assert( ie01 != -1 );
    int ie12 = findEdge(ip1,ip2, psup_ind,psup); assert( ie12 != -1 );
    int ie23 = findEdge(ip2,ip3, psup_ind,psup); assert( ie23 != -1 );
    int ie30 = findEdge(ip3,ip0, psup_ind,psup); assert( ie30 != -1 );
    int ip01 = ie01 + nPoint0;
    int ip12 = ie12 + nPoint0;
    int ip23 = ie23 + nPoint0;
    int ip30 = ie30 + nPoint0;
    int ip0123 = iq + nPoint0 + ne0;
    aQuad1.push_back(ip0);   aQuad1.push_back(ip01); aQuad1.push_back(ip0123); aQuad1.push_back(ip30);
    aQuad1.push_back(ip1);   aQuad1.push_back(ip12); aQuad1.push_back(ip0123); aQuad1.push_back(ip01);
    aQuad1.push_back(ip2);   aQuad1.push_back(ip23); aQuad1.push_back(ip0123); aQuad1.push_back(ip12);
    aQuad1.push_back(ip3);   aQuad1.push_back(ip30); aQuad1.push_back(ip0123); aQuad1.push_back(ip23);
  }
}


// new points is in the order of [old points], [edge points]
void TetSubdiv
(std::vector<unsigned int>& aTet1,
 std::vector<int>& psup_ind,
 std::vector<int>& psup,
 const unsigned int* aTet0, int nTet0,
 unsigned int nPoint0)
{
  const int nt0 = nTet0;
  std::vector<int> elsup_ind, elsup;
  JArrayElemSurPoint_MeshElem(elsup_ind,elsup,
                           aTet0,nTet0,4,nPoint0);
  JArrayEdge_MeshElem(psup_ind,psup,
                       aTet0, MESHELEM_TET, elsup_ind, elsup,
                       false);
  aTet1.resize(0);
  aTet1.reserve(nTet0*4);
  for(int it=0;it<nt0;++it){
    int ip0 = aTet0[it*4+0];
    int ip1 = aTet0[it*4+1];
    int ip2 = aTet0[it*4+2];
    int ip3 = aTet0[it*4+3];
    int ie01 = findEdge(ip0,ip1, psup_ind,psup); assert( ie01 != -1 );
    int ie02 = findEdge(ip0,ip2, psup_ind,psup); assert( ie02 != -1 );
    int ie03 = findEdge(ip0,ip3, psup_ind,psup); assert( ie03 != -1 );
    int ie12 = findEdge(ip1,ip2, psup_ind,psup); assert( ie12 != -1 );
    int ie13 = findEdge(ip1,ip3, psup_ind,psup); assert( ie13 != -1 );
    int ie23 = findEdge(ip2,ip3, psup_ind,psup); assert( ie23 != -1 );
    int ip01 = ie01 + nPoint0;
    int ip02 = ie02 + nPoint0;
    int ip03 = ie03 + nPoint0;
    int ip12 = ie12 + nPoint0;
    int ip13 = ie13 + nPoint0;
    int ip23 = ie23 + nPoint0;
    aTet1.push_back(ip0);  aTet1.push_back(ip01); aTet1.push_back(ip02); aTet1.push_back(ip03);
    aTet1.push_back(ip1);  aTet1.push_back(ip01); aTet1.push_back(ip13); aTet1.push_back(ip12);
    aTet1.push_back(ip2);  aTet1.push_back(ip02); aTet1.push_back(ip12); aTet1.push_back(ip23);
    aTet1.push_back(ip3);  aTet1.push_back(ip03); aTet1.push_back(ip23); aTet1.push_back(ip13);
    aTet1.push_back(ip01); aTet1.push_back(ip23); aTet1.push_back(ip13); aTet1.push_back(ip12);
    aTet1.push_back(ip01); aTet1.push_back(ip23); aTet1.push_back(ip12); aTet1.push_back(ip02);
    aTet1.push_back(ip01); aTet1.push_back(ip23); aTet1.push_back(ip02); aTet1.push_back(ip03);
    aTet1.push_back(ip01); aTet1.push_back(ip23); aTet1.push_back(ip03); aTet1.push_back(ip13);
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


/////////////////////////////////////////////////

void HexSubdiv
(std::vector<unsigned int>& aHex1,
 std::vector<int>& psupIndHex0,
 std::vector<int>& psupHex0,
 std::vector<unsigned int>& aQuadHex0,
 ///
 const unsigned int* aHex0, int nHex0,
 const int nhp0)
{
  //  int nhp0 = (int)aHexPoint0.size(); // hex point
  std::vector<int> elsupIndHex0, elsupHex0;
  JArrayElemSurPoint_MeshElem(elsupIndHex0, elsupHex0,
                           aHex0,nHex0,8,nhp0);
  
  //edge
  JArrayEdge_MeshElem(psupIndHex0, psupHex0,
                       aHex0, MESHELEM_HEX, elsupIndHex0,elsupHex0,
                       false); // is_directional = false
  
  //face
  aQuadHex0.clear();
  {
    std::vector<int> aHexSurRel0;
    makeSurroundingRelationship(aHexSurRel0,
                                aHex0,nHex0,8,
                                elsupIndHex0,elsupHex0,
                                nFaceElem(MESHELEM_HEX),
                                nNodeElemFace(MESHELEM_HEX, 0),
                                noelElemFace(MESHELEM_HEX));
    for(unsigned int ih=0;ih<(unsigned int)nHex0;++ih){
      for(int ifh=0;ifh<6;++ifh){
        int jh0 = aHexSurRel0[ih*6*2+ifh*2+0];
        if( jh0!=-1 && (int)ih>jh0 ) continue;
        for(int inofa=0;inofa<4;++inofa){
          int inoel0 = noelElemFace_Hex[ifh][inofa];
          int igp0 = aHex0[ih*8+inoel0];
          aQuadHex0.push_back(igp0);
        }
      }
    }
  }
  std::vector<int> elsupIndQuadHex0, elsupQuadHex0;
  JArrayElemSurPoint_MeshElem(elsupIndQuadHex0,elsupQuadHex0,
                           aQuadHex0.data(),aQuadHex0.size()/4,4,nhp0);
  
  const int neh0 = (int)psupHex0.size();
  const int nfh0 = (int)aQuadHex0.size()/4;
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
  for(int ih=0;ih<nHex0;++ih){
    int ihc0 = aHex0[ih*8+0];
    int ihc1 = aHex0[ih*8+1];
    int ihc2 = aHex0[ih*8+2];
    int ihc3 = aHex0[ih*8+3];
    int ihc4 = aHex0[ih*8+4];
    int ihc5 = aHex0[ih*8+5];
    int ihc6 = aHex0[ih*8+6];
    int ihc7 = aHex0[ih*8+7];
    int ihc01 = findEdge(ihc0,ihc1, psupIndHex0,psupHex0)+nhp0; assert(ihc01>=nhp0&&ihc01<nhp0+neh0);
    int ihc12 = findEdge(ihc1,ihc2, psupIndHex0,psupHex0)+nhp0; assert(ihc12>=nhp0&&ihc12<nhp0+neh0);
    int ihc23 = findEdge(ihc2,ihc3, psupIndHex0,psupHex0)+nhp0; assert(ihc23>=nhp0&&ihc23<nhp0+neh0);
    int ihc30 = findEdge(ihc3,ihc0, psupIndHex0,psupHex0)+nhp0; assert(ihc30>=nhp0&&ihc30<nhp0+neh0);
    int ihc45 = findEdge(ihc4,ihc5, psupIndHex0,psupHex0)+nhp0; assert(ihc45>=nhp0&&ihc45<nhp0+neh0);
    int ihc56 = findEdge(ihc5,ihc6, psupIndHex0,psupHex0)+nhp0; assert(ihc56>=nhp0&&ihc56<nhp0+neh0);
    int ihc67 = findEdge(ihc6,ihc7, psupIndHex0,psupHex0)+nhp0; assert(ihc67>=nhp0&&ihc67<nhp0+neh0);
    int ihc74 = findEdge(ihc7,ihc4, psupIndHex0,psupHex0)+nhp0; assert(ihc74>=nhp0&&ihc74<nhp0+neh0);
    int ihc04 = findEdge(ihc0,ihc4, psupIndHex0,psupHex0)+nhp0; assert(ihc04>=nhp0&&ihc04<nhp0+neh0);
    int ihc15 = findEdge(ihc1,ihc5, psupIndHex0,psupHex0)+nhp0; assert(ihc15>=nhp0&&ihc15<nhp0+neh0);
    int ihc26 = findEdge(ihc2,ihc6, psupIndHex0,psupHex0)+nhp0; assert(ihc26>=nhp0&&ihc26<nhp0+neh0);
    int ihc37 = findEdge(ihc3,ihc7, psupIndHex0,psupHex0)+nhp0; assert(ihc37>=nhp0&&ihc37<nhp0+neh0);
    int ihc0473 = findFace(ihc0,ihc4,ihc7,ihc3, aQuadHex0,elsupIndQuadHex0,elsupQuadHex0)+nhp0+neh0; assert(ihc0473>=nhp0+neh0&&ihc0473<nhp0+neh0+nfh0);
    int ihc1265 = findFace(ihc1,ihc2,ihc6,ihc5, aQuadHex0,elsupIndQuadHex0,elsupQuadHex0)+nhp0+neh0; assert(ihc1265>=nhp0+neh0&&ihc1265<nhp0+neh0+nfh0);
    int ihc0154 = findFace(ihc0,ihc1,ihc5,ihc4, aQuadHex0,elsupIndQuadHex0,elsupQuadHex0)+nhp0+neh0; assert(ihc0154>=nhp0+neh0&&ihc0154<nhp0+neh0+nfh0);
    int ihc3762 = findFace(ihc3,ihc7,ihc6,ihc2, aQuadHex0,elsupIndQuadHex0,elsupQuadHex0)+nhp0+neh0; assert(ihc3762>=nhp0+neh0&&ihc3762<nhp0+neh0+nfh0);
    int ihc0321 = findFace(ihc0,ihc3,ihc2,ihc1, aQuadHex0,elsupIndQuadHex0,elsupQuadHex0)+nhp0+neh0; assert(ihc0321>=nhp0+neh0&&ihc0321<nhp0+neh0+nfh0);
    int ihc4567 = findFace(ihc4,ihc5,ihc6,ihc7, aQuadHex0,elsupIndQuadHex0,elsupQuadHex0)+nhp0+neh0; assert(ihc4567>=nhp0+neh0&&ihc4567<nhp0+neh0+nfh0);
    int ihc01234567 = ih + nhp0 + neh0 + nfh0;
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

