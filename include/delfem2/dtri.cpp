/*
 * Copyright (c) 2019 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */


#include <algorithm>
#include <map>
#include <stack>
#include <iostream>
#include <cassert>

#include "delfem2/dtri.h"

// ---------------------------------------------


DFM2_INLINE void delfem2::MakeInnerRelationTri
(std::vector<CDynTri>& aTri, const unsigned int npoin,
 const std::vector<int>& elsup_ind,
 const std::vector<int>& elsup)
{
 
  std::vector<int> tmp_poin(npoin,0);
  unsigned int inpofa[2];
  
  const unsigned int nTri = (int)aTri.size();
  for(unsigned int itri=0;itri<nTri;itri++){
    for(unsigned int iedtri=0;iedtri<3;iedtri++){
      for(unsigned int ipoed=0;ipoed<2;ipoed++){
        inpofa[ipoed] = aTri[itri].v[ (iedtri+1+ipoed)%3 ];
        tmp_poin[ inpofa[ipoed] ] = 1;
      }
      const unsigned int ipoin0= inpofa[0];
      bool iflg = false;
      for(int ielsup=elsup_ind[ipoin0];ielsup<elsup_ind[ipoin0+1];ielsup++){
        const unsigned int jtri0 = elsup[ielsup];
        if( jtri0 == itri ) continue;
        for(unsigned int jedtri=0;jedtri<3;jedtri++){
          iflg = true;
          for(unsigned int jpoed=0;jpoed<2;jpoed++){
            const unsigned int jpoin0 =  aTri[jtri0].v[ (jedtri+1+jpoed)%3 ];
            if( tmp_poin[ jpoin0 ] == 0 ){ iflg = false; break; }
          }
          if( iflg ){
            //            aTri[itri].g2[iedtri] = -2;
            aTri[itri].s2[iedtri] = jtri0;
//            aTri[itri].r2[iedtri] = EdEd2Rel[iedtri][jedtri];
            break;
          }
        }
        if( iflg ) break;
      }
      if( !iflg ){
        aTri[itri].s2[iedtri] = -1;
      }
      for(unsigned int ipofa : inpofa){
        tmp_poin[ ipofa ] = 0;
      }
    }
  }
}

DFM2_INLINE bool delfem2::JArray_MakeElSuP
(std::vector<int>& elsup_ind, std::vector<int>& elsup,
 const std::vector<CDynTri>& aTri, const unsigned int npoin)
{
  elsup_ind.assign(npoin+1, 0);
  for(const auto & itri : aTri){
    for(int inotri : itri.v){
      elsup_ind[ inotri+1 ]++;
    }
  }
  for(unsigned int ipoin=0;ipoin<npoin;ipoin++){
    elsup_ind[ipoin+1] += elsup_ind[ipoin];
  }
  const int nelsup = elsup_ind[npoin];
  elsup.resize(nelsup);
  for(std::size_t itri=0;itri<aTri.size();itri++){
    for(unsigned int ipoin0 : aTri[itri].v){
      const unsigned int ielsup = elsup_ind[ipoin0];
      elsup[ielsup] = itri;
      elsup_ind[ipoin0]++;
    }
  }
  for(unsigned int ipoin=npoin;ipoin>0;ipoin--){
    elsup_ind[ipoin] = elsup_ind[ipoin-1];
  }
  elsup_ind[0] = 0;
  return true;
}

DFM2_INLINE void delfem2::JArray_PSuP
(std::vector<int>& psup_ind, std::vector<int>& psup,
 const std::vector<CDynTri>& aTri, const unsigned int npoin,
 const std::vector<int>& elsup_ind, const std::vector<int>& elsup)
{
  std::vector<unsigned int> aflg(npoin,0);
  psup_ind[0] = 0;
  for (unsigned int ino = 0; ino<npoin; ino++){
    psup_ind[ino+1] = psup_ind[ino];
    aflg[ino] = ino;
    for (int ielsup = elsup_ind[ino]; ielsup<elsup_ind[ino+1]; ielsup++){
      unsigned int itri1 = elsup[ielsup];
      for (unsigned int ino1 : aTri[itri1].v){
        if (aflg[ino1]==ino) continue;
        psup_ind[ino+1]++;
        aflg[ino1] = ino;
      }
    }
  }
  const int npsup = psup_ind[npoin];
  psup.resize(npsup);
  for (unsigned int ino = 0; ino<npoin; ino++){ aflg[ino] = 0; }
  unsigned int iedge = 0;
  for (unsigned int ino = 0; ino<npoin; ino++){
    assert(psup_ind[ino]==(int)iedge);
    aflg[ino] = ino;
    for (int ielsup = elsup_ind[ino]; ielsup<elsup_ind[ino+1]; ielsup++){
      unsigned int itri1 = elsup[ielsup];
      for (unsigned int ino1 : aTri[itri1].v){
        if (aflg[ino1]==ino) continue;
        psup[iedge] = ino1;
        iedge++;
        aflg[ino1] = ino;
      }
    }
  }
  assert((int)iedge==npsup);
}


// -----------------------------------------------------------------


DFM2_INLINE bool delfem2::InsertPoint_ElemEdge
(const int ipo_ins,    //the index of the new point
 const int itri_ins,  //triangle index
 const int ied_ins,  //edge index
 std::vector<CDynPntSur>& aPo,
 std::vector<CDynTri>& aTri )
{
  assert(itri_ins < (int)aTri.size());
  assert(ipo_ins < (int)aPo.size());
  
  if (aTri[itri_ins].s2[ied_ins]==-1){
    assert(0);
  }
 
  const int itri_adj = aTri[itri_ins].s2[ied_ins];
  const int ied_adj = FindAdjEdgeIndex(aTri[itri_ins],ied_ins,aTri);
  assert(itri_adj < (int)aTri.size());
  assert(ied_ins < 3);
  
  const int itri0 = itri_ins;
  const int itri1 = itri_adj;
  const int itri2 = (int)aTri.size();
  const int itri3 = (int)aTri.size()+1;
  
  aTri.resize(aTri.size()+2);
  
  CDynTri oldA = aTri[itri_ins];
  CDynTri oldB = aTri[itri_adj];
  
  const int inoA0 = ied_ins;
  const int inoA1 = (ied_ins+1)%3;
  const int inoA2 = (ied_ins+2)%3;
  
  const int inoB0 = ied_adj;
  const int inoB1 = (ied_adj+1)%3;
  const int inoB2 = (ied_adj+2)%3;
  
  assert(oldA.v[inoA1]==oldB.v[inoB2]);
  assert(oldA.v[inoA2]==oldB.v[inoB1]);
  assert(oldA.s2[inoA0]==itri1);
  assert(oldB.s2[inoB0]==itri0);
  
  aPo[ipo_ins].e = itri0;        aPo[ipo_ins].d = 0;
  aPo[oldA.v[inoA2]].e = itri0;  aPo[oldA.v[inoA2]].d = 1;
  aPo[oldA.v[inoA0]].e = itri1;  aPo[oldA.v[inoA0]].d = 1;
  aPo[oldB.v[inoB2]].e = itri2;  aPo[oldB.v[inoB2]].d = 1;
  aPo[oldB.v[inoB0]].e = itri3;  aPo[oldB.v[inoB0]].d = 1;
  
  {
    CDynTri& tri0 = aTri[itri0];
    tri0.v[0] = ipo_ins;
    tri0.v[1] = oldA.v[inoA2];
    tri0.v[2] = oldA.v[inoA0];
    tri0.s2[0] = oldA.s2[inoA1];
    tri0.s2[1] = itri1;
    tri0.s2[2] = itri3;
  }
  if (oldA.s2[inoA1]>=0&&oldA.s2[inoA1]<(int)aTri.size()){
    unsigned int jt0 = oldA.s2[inoA1];
    assert( jt0 < (int)aTri.size());
    unsigned int jno0 = FindAdjEdgeIndex(oldA, inoA1, aTri);
    aTri[jt0].s2[jno0] = itri0;
  }
  
  {
    CDynTri& tri1 = aTri[itri1];
    tri1.v[0] = ipo_ins;
    tri1.v[1] = oldA.v[inoA0];
    tri1.v[2] = oldA.v[inoA1];
    tri1.s2[0] = oldA.s2[inoA2];
    tri1.s2[1] = itri2;
    tri1.s2[2] = itri0;
  }
  if (oldA.s2[inoA2]>=0&&oldA.s2[inoA2]<(int)aTri.size()){
    unsigned int jt0 = oldA.s2[inoA2];
    assert( jt0 < (int)aTri.size());
    unsigned int jno0 = FindAdjEdgeIndex(oldA, inoA2, aTri);
    aTri[jt0].s2[jno0] = itri1;
  }
  
  {
    CDynTri& tri2 = aTri[itri2];
    tri2.v[0] = ipo_ins;
    tri2.v[1] = oldB.v[inoB2];
    tri2.v[2] = oldB.v[inoB0];
    tri2.s2[0] = oldB.s2[inoB1];
    tri2.s2[1] = itri3;
    tri2.s2[2] = itri1;
  }
  if(oldB.s2[inoB1]>=0&&oldB.s2[inoB1]<(int)aTri.size()){
    unsigned int jt0 = oldB.s2[inoB1];
    assert( jt0 < (int)aTri.size());
    unsigned int jno0 = FindAdjEdgeIndex(oldB, inoB1, aTri);
     aTri[jt0].s2[jno0] = itri2;
  }
  
  {
    CDynTri& tri3 = aTri[itri3];
    tri3.v[0] = ipo_ins;
    tri3.v[1] = oldB.v[inoB0];
    tri3.v[2] = oldB.v[inoB1];
    tri3.s2[0] = oldB.s2[inoB2];
    tri3.s2[1] = itri0;
    tri3.s2[2] = itri2;
  }
  if (oldB.s2[inoB2]>=0&&oldB.s2[inoB2]<(int)aTri.size()){
    unsigned int jt0 = oldB.s2[inoB2];
    assert( jt0 < (int)aTri.size());
    unsigned int jno0 = FindAdjEdgeIndex(oldB, inoB2, aTri);
    aTri[jt0].s2[jno0] = itri3;
  }
  
  return true;
}



DFM2_INLINE bool delfem2::InsertPoint_Elem
(const int ipo_ins,
 const int itri_ins,
 std::vector<CDynPntSur>& aPo,
 std::vector<CDynTri>& aTri)
{
  assert(itri_ins < (int)aTri.size());
  assert(ipo_ins < (int)aPo.size());
  
  const int itA = itri_ins;
  const int itB = (int)aTri.size();
  const int itC = (int)aTri.size()+1;
  
  aTri.resize(aTri.size()+2);
  const CDynTri old = aTri[itri_ins];
  
  aPo[ipo_ins].e = itA;       aPo[ipo_ins].d = 0;
  aPo[old.v[0]].e = itB;  aPo[old.v[0]].d = 2;
  aPo[old.v[1]].e = itC;  aPo[old.v[1]].d = 2;
  aPo[old.v[2]].e = itA;  aPo[old.v[2]].d = 2;
  
  {
    CDynTri& triA = aTri[itA];
    triA.v[0] = ipo_ins;
    triA.v[1] = old.v[1];
    triA.v[2] = old.v[2];
    triA.s2[0] = old.s2[0];
    triA.s2[1] = itB;
    triA.s2[2] = itC;
  }
  if (old.s2[0]>=0&&old.s2[0]<(int)aTri.size()){
    const unsigned int jt0 = old.s2[0];
    const unsigned int jno0 = FindAdjEdgeIndex(old, 0, aTri);
    aTri[jt0].s2[jno0] = itA;
  }
  
  {
    CDynTri& triB = aTri[itB];
    triB.v[0] = ipo_ins;
    triB.v[1] = old.v[2];
    triB.v[2] = old.v[0];
    triB.s2[0] = old.s2[1];
    triB.s2[1] = itC;
    triB.s2[2] = itA;
  }
  if (old.s2[1]>=0&&old.s2[1]<(int)aTri.size()){
    const unsigned int jt0 = old.s2[1];
    const unsigned int jno0 = FindAdjEdgeIndex(old, 1, aTri);
    aTri[jt0].s2[jno0] = itB;
  }

  {
    CDynTri& triC = aTri[itC];
    triC.v[0] = ipo_ins;
    triC.v[1] = old.v[0];
    triC.v[2] = old.v[1];
    triC.s2[0] = old.s2[2];
    triC.s2[1] = itA;
    triC.s2[2] = itB;
  }
  if (old.s2[2]>=0&&old.s2[2]<(int)aTri.size()){
    const unsigned int jt0 = old.s2[2];
    const unsigned int jno0 = FindAdjEdgeIndex(old, 2, aTri);
    aTri[jt0].s2[jno0] = itC;
    assert( jt0 < aTri.size() );
  }
  
  return true;
}


DFM2_INLINE bool delfem2::FlipEdge(
    unsigned int itriA,
    unsigned int ied0,
    std::vector<CDynPntSur>& aPo,
    std::vector<CDynTri>& aTri)
{
  assert(itriA < aTri.size());
  assert(ied0 < 3);
  assert(aTri[itriA].s2[ied0]>=0&&aTri[itriA].s2[ied0]<(int)aTri.size());
  
  const int itriB = aTri[itriA].s2[ied0];
  const unsigned int ied1 = FindAdjEdgeIndex(aTri[itriA], ied0, aTri);
  assert(itriB < (int)aTri.size());
  assert(ied1 < 3);
  assert(aTri[itriB].s2[ied1]>=0&&aTri[itriA].s2[ied0]<(int)aTri.size());
  
  const CDynTri oldA = aTri[itriA];
  const CDynTri oldB = aTri[itriB];
  
  const unsigned int noA0 = ied0;
  const unsigned int noA1 = (ied0+1)%3;
  const unsigned int noA2 = (ied0+2)%3;
  
  const unsigned int noB0 = ied1;
  const unsigned int noB1 = (ied1+1)%3;
  const unsigned int noB2 = (ied1+2)%3;

  assert(oldA.v[noA1]==oldB.v[noB2]);
  assert(oldA.v[noA2]==oldB.v[noB1]);
  
  aPo[oldA.v[noA1]].e = itriA;  aPo[oldA.v[noA1]].d = 0;
  aPo[oldA.v[noA0]].e = itriA;  aPo[oldA.v[noA0]].d = 2;
  aPo[oldB.v[noB1]].e = itriB;  aPo[oldB.v[noB1]].d = 0;
  aPo[oldB.v[noB0]].e = itriB;  aPo[oldB.v[noB0]].d = 2;
  
  {
    CDynTri& triA = aTri[itriA];
    triA.v[0] = oldA.v[noA1];
    triA.v[1] = oldB.v[noB0];
    triA.v[2] = oldA.v[noA0];
    triA.s2[0] = itriB;
    triA.s2[1] = oldA.s2[noA2];
    triA.s2[2] = oldB.s2[noB1];
  }
  if (oldA.s2[noA2]>=0){
    const unsigned int jt0 = oldA.s2[noA2];
    assert( jt0 < aTri.size() && jt0 != itriB && jt0 != itriA );
    const unsigned int jno0 = FindAdjEdgeIndex(oldA, noA2, aTri);
    aTri[jt0].s2[jno0] = (int)itriA;
  }
  if (oldB.s2[noB1]>=0){
    const unsigned int jt0 = oldB.s2[noB1];
    assert( jt0 < aTri.size() && jt0 != itriB && jt0 != itriA );
    const unsigned int jno0 = FindAdjEdgeIndex(oldB, noB1, aTri);
    aTri[jt0].s2[jno0] = (int)itriA;
  }
  
  {
    CDynTri& triB = aTri[itriB];
    triB.v[0] = oldB.v[noB1];
    triB.v[1] = oldA.v[noA0];
    triB.v[2] = oldB.v[noB0];
    triB.s2[0] = (int)itriA;
    triB.s2[1] = oldB.s2[noB2];
    triB.s2[2] = oldA.s2[noA1];
  }
  if (oldB.s2[noB2]>=0){
    const unsigned int jt0 = oldB.s2[noB2];
    assert( jt0 < aTri.size());
    const unsigned int jno0 = FindAdjEdgeIndex(oldB, noB2, aTri);
    aTri[jt0].s2[jno0] = itriB;
  }
  if (oldA.s2[noA1]>=0){
    const unsigned int jt0 = oldA.s2[noA1];
    assert( jt0 < aTri.size() );
    const unsigned int jno0 = FindAdjEdgeIndex(oldA, noA1, aTri);
    aTri[jt0].s2[jno0] = itriB;
  }
  return true;
}


DFM2_INLINE bool delfem2::FindEdge_LookAroundPoint
(unsigned int &itri0,unsigned int &inotri0, unsigned &inotri1,
 //
 const int ipo0, const int ipo1,
 const std::vector<CDynPntSur>& aPo,
 const std::vector<CDynTri>& aTri)
{
  const int itri_ini = aPo[ipo0].e;
  const unsigned int inotri_ini = aPo[ipo0].d;
  unsigned int inotri_cur = inotri_ini;
  int itri_cur = itri_ini;
  for (;;){  // serch clock-wise
    assert(aTri[itri_cur].v[inotri_cur]==ipo0);
    {
      const unsigned int inotri2 = (inotri_cur+1)%3; // indexRot3[1][inotri_cur];
      if (aTri[itri_cur].v[inotri2]==ipo1){
        itri0 = itri_cur;
        inotri0 = inotri_cur;
        inotri1 = inotri2;
        assert(aTri[itri0].v[inotri0]==ipo0);
        assert(aTri[itri0].v[inotri1]==ipo1);
        return true;
      }
    }
    {
      const unsigned int inotri2 = (inotri_cur+2)%3; // indexRot3[2][inotri_cur];
      if (aTri[itri_cur].s2[inotri2]==-1){ break; }
      const int itri_nex = aTri[itri_cur].s2[inotri2];
      const unsigned int jno = FindAdjEdgeIndex(aTri[itri_cur], inotri2, aTri);
      const unsigned int inotri3 = (jno+2)%3; // rel[inotri_cur];
      assert(aTri[itri_nex].v[inotri3]==ipo0);
      if (itri_nex==itri_ini) return false;
      itri_cur = itri_nex;
      inotri_cur = inotri3;
    }
  }
  
  inotri_cur = inotri_ini;
  itri_cur = itri_ini;
  for (;;){
    assert(aTri[itri_cur].v[inotri_cur]==ipo0);
    {
      const unsigned int inotri2 = (inotri_cur+1)%3; // indexRot3[1][inotri_cur];
      if (aTri[itri_cur].s2[inotri2]==-1){ break; }
      const int itri_nex = aTri[itri_cur].s2[inotri2];
      const unsigned int jno = FindAdjEdgeIndex(aTri[itri_cur], inotri2, aTri);
      const unsigned int inotri3 = (jno+1)%3; // rel[inotri_cur];
      assert(aTri[itri_nex].v[inotri3]==ipo0);
      if (itri_nex==itri_ini){  // end if it goes around
        itri0 = 0;
        inotri0 = 0; inotri1 = 0;
        return false;
      }
      itri_cur = itri_nex;
      inotri_cur = inotri3;
    }
    {
      const unsigned int inotri2 = (inotri_cur+1)%3; // indexRot3[1][inotri_cur];
      if (aTri[itri_cur].v[inotri2]==ipo1){
        itri0 = itri_cur;
        inotri0 = inotri_cur;
        inotri1 = inotri2;
        assert(aTri[itri0].v[inotri0]==ipo0);
        assert(aTri[itri0].v[inotri1]==ipo1);
        return true;
      }
    }
  }
  return false;
}



DFM2_INLINE bool delfem2::FindEdge_LookAllTriangles
(int& itri0, int& iedtri0,
 //
 const int ipo0, const int ipo1,
 const std::vector<CDynTri>& aTri)
{
  for(size_t itri=0;itri<aTri.size();++itri){
    for(int iedtri=0;iedtri<3;++iedtri){
      int jpo0 = aTri[itri].v[(iedtri+0)%3];
      int jpo1 = aTri[itri].v[(iedtri+1)%3];
      if( jpo0 == ipo0 && jpo1 == ipo1 ){
        itri0 = itri;
        iedtri0 = iedtri;
        return true;
      }
    }
  }
  return false;
}

DFM2_INLINE void delfem2::AssertDTri(
    const std::vector<CDynTri>& aTri )
{
	const unsigned int ntri = aTri.size();
	for(unsigned int itri=0;itri<ntri;itri++){
    const CDynTri& tri = aTri[itri];
    if( tri.v[0] == -1 ){
      assert(tri.v[1] == -1);
      assert(tri.v[2] == -1);
      continue;
    }
    assert( tri.v[0] != tri.v[1] );
    assert( tri.v[1] != tri.v[2] );
    assert( tri.v[2] != tri.v[0] );
    assert( (tri.s2[0] != tri.s2[1]) || tri.s2[0] == -1 );
    assert( (tri.s2[1] != tri.s2[2]) || tri.s2[1] == -1 );
    assert( (tri.s2[2] != tri.s2[0]) || tri.s2[0] == -1 );
		for(int iedtri=0;iedtri<3;iedtri++){
      if( tri.s2[iedtri] == -1 ) continue;
      assert( tri.s2[iedtri] >=0 && tri.s2[iedtri]<(int)aTri.size() );
      const unsigned int jtri = tri.s2[iedtri]; assert( jtri < ntri );
      const unsigned int jno = FindAdjEdgeIndex(aTri[itri], iedtri, aTri);
      assert( aTri[jtri].s2[jno] == itri );
      assert( aTri[itri].v[(iedtri+1)%3] == aTri[jtri].v[(jno+2)%3] );
      assert( aTri[itri].v[(iedtri+2)%3] == aTri[jtri].v[(jno+1)%3] );
		}
	}
}

DFM2_INLINE void delfem2::AssertMeshDTri
(const std::vector<CDynPntSur>& aPo3D,
 const std::vector<CDynTri>& aSTri)
{
  const unsigned int npo = aPo3D.size();
  const unsigned int ntri = aSTri.size();
  for (unsigned int itri = 0; itri<ntri; itri++){
    const CDynTri& tri0 = aSTri[itri];
    assert( tri0.v[0] < npo);
    assert( tri0.v[0] < npo);
    assert( tri0.v[0] < npo);
  }
  for (int ipoin = 0; ipoin<npo; ++ipoin){
    const int itri0 = aPo3D[ipoin].e;
    const unsigned int inoel0 = aPo3D[ipoin].d;
    if (aPo3D[ipoin].e>=0){
      assert(aPo3D[ipoin].d < 3);
      if (aSTri[itri0].v[inoel0]!=ipoin){}
      assert(aSTri[itri0].v[inoel0]==ipoin);
    }
  }
}


DFM2_INLINE void delfem2::InitializeMesh
(std::vector<CDynPntSur>& aPo3D,
 std::vector<CDynTri>& aSTri,
 //
 const unsigned int* aTri, int nTri,
 int nXYZ)
{
  aPo3D.resize(nXYZ);
  for (int ipo = 0; ipo<nXYZ; ++ipo){
    aPo3D[ipo].e = -1; // for unreffered point
    aPo3D[ipo].d = 0;
  }
  aSTri.resize(nTri);
  for (int itri = 0; itri<nTri; itri++){
    aSTri[itri].v[0] = aTri[itri*3+0];
    aSTri[itri].v[1] = aTri[itri*3+1];
    aSTri[itri].v[2] = aTri[itri*3+2];
  }
  for (int itri = 0; itri<nTri; itri++){
    unsigned int i1 = aSTri[itri].v[0];
    unsigned int i2 = aSTri[itri].v[1];
    unsigned int i3 = aSTri[itri].v[2];
    aPo3D[i1].e = itri; aPo3D[i1].d = 0;
    aPo3D[i2].e = itri; aPo3D[i2].d = 1;
    aPo3D[i3].e = itri; aPo3D[i3].d = 2;
  }
  {
    std::vector<int> elsup_ind, elsup;
    JArray_MakeElSuP(elsup_ind, elsup,
                     aSTri, (int)aPo3D.size());
    MakeInnerRelationTri(aSTri, (int)aPo3D.size(),
                         elsup_ind,elsup);
  }
}


// -----------------------------------------------------------------

DFM2_INLINE void delfem2::MoveCCW
(int& itri_cur,
 unsigned int &inotri_cur,
 bool& flag_is_wall,
 //
 const std::vector<CDynTri>& aTri)
{
  const unsigned int inotri1 = (inotri_cur+1)%3; // indexRot3[1][inotri_cur];
  if (aTri[itri_cur].s2[inotri1]==-1){ flag_is_wall = true; return; }
  flag_is_wall = false;
  const int itri_nex = aTri[itri_cur].s2[inotri1];
  unsigned int ino2 = FindAdjEdgeIndex(aTri[itri_cur],inotri1,aTri);
  const unsigned int inotri_nex = (ino2+1)%3;
  itri_cur = itri_nex;
  inotri_cur = inotri_nex;
}


DFM2_INLINE bool delfem2::DeleteTri
(int itri_to,
 std::vector<CDynPntSur>& aPo,
 std::vector<CDynTri>& aTri)
{
  if (itri_to>=(int)aTri.size()) return true;
  {
    assert(aTri[itri_to].s2[0]==-1);
    assert(aTri[itri_to].s2[1]==-1);
    assert(aTri[itri_to].s2[2]==-1);
  }
  const int itri_from = (int)aTri.size()-1;
  if (itri_to==itri_from){
    aTri.resize(aTri.size()-1);
    return true;
  }
  aTri[itri_to] = aTri[itri_from];
  aTri.resize(aTri.size()-1);
  for (int iedtri = 0; iedtri<3; iedtri++){
    if (aTri[itri_to].s2[iedtri]==-1) continue;
    const unsigned int itri_adj = aTri[itri_to].s2[iedtri];
    const unsigned int iedtri_adj = FindAdjEdgeIndex(aTri[itri_to],iedtri,aTri);
    assert(itri_adj < aTri.size());
    assert(aTri[itri_adj].s2[iedtri_adj]==itri_from);
    aTri[itri_adj].s2[iedtri_adj] = itri_to;
  }
  for (unsigned int inotri = 0; inotri<3; inotri++){
    const unsigned int ipo0 = aTri[itri_to].v[inotri];
    aPo[ipo0].e = itri_to;
    aPo[ipo0].d = inotri;
  }
  return true;
}

DFM2_INLINE unsigned int delfem2::FindAdjEdgeIndex
 (const CDynTri& t0,
  unsigned int ied0,
  const std::vector<delfem2::CDynTri>& aTri)
{
  int iv0 = t0.v[(ied0+1)%3];
  int iv1 = t0.v[(ied0+2)%3];
  assert( iv0 >=0 && iv1 >= 0 && iv0 != iv1 );
  assert( t0.s2[ied0] >= 0 );
  unsigned int it1 = t0.s2[ied0];
  if( aTri[it1].v[1] == iv1 && aTri[it1].v[2] == iv0 ){ return 0; }
  if( aTri[it1].v[2] == iv1 && aTri[it1].v[0] == iv0 ){ return 1; }
  if( aTri[it1].v[0] == iv1 && aTri[it1].v[1] == iv0 ){ return 2; }
  assert(0);
  return 0;
}


DFM2_INLINE bool delfem2::CollapseEdge_MeshDTri
(const unsigned int itri_del,
 const unsigned int ied_del,
 std::vector<CDynPntSur>& aPo,
 std::vector<CDynTri>& aTri)
{
  assert( itri_del < aTri.size() );
  assert( ied_del < 3 );
  if (aTri[itri_del].s2[ied_del]==-1){
    std::cout<<"Error!-->Not Implemented: Mesh with hole"<<std::endl;
    assert(0);
  }
  
  const unsigned int itri_adj = aTri[itri_del].s2[ied_del];
  const unsigned int ied_adj = FindAdjEdgeIndex(aTri[itri_del], ied_del, aTri);
  assert(itri_adj < aTri.size());
  assert(ied_adj < 3);
  assert(aTri[itri_adj].s2[ied_adj]==itri_del);
  
  // do nothing and return false if the collapsing edge is on the boundary
  if( aTri[itri_del].s2[(ied_del+1)%3] == -1 ) return false;
  if( aTri[itri_del].s2[(ied_del+2)%3] == -1 ) return false;
  if( aTri[itri_adj].s2[(ied_adj+1)%3] == -1 ) return false;
  if( aTri[itri_adj].s2[(ied_adj+2)%3] == -1 ) return false;
  
  const unsigned int itA = itri_del;
  const unsigned int itB = itri_adj;
  const unsigned int itC = aTri[itA].s2[(ied_del+1)%3];
  const unsigned int itD = aTri[itA].s2[(ied_del+2)%3];
  const unsigned int itE = aTri[itB].s2[(ied_adj+1)%3];
  const unsigned int itF = aTri[itB].s2[(ied_adj+2)%3];
  
  const unsigned int inoA0 = ied_del;
  const unsigned int inoA1 = (inoA0+1)%3;
  const unsigned int inoA2 = (inoA0+2)%3; // point to be deleted
  
  const unsigned int inoB0 = ied_adj;
  const unsigned int inoB1 = (inoB0+1)%3; // point to be deleted
  const unsigned int inoB2 = (inoB0+2)%3;
  
  const unsigned int inoC0 = FindAdjEdgeIndex(aTri[itA], inoA1, aTri);
  const unsigned int inoC1 = (inoC0+1)%3;
  const unsigned int inoC2 = (inoC0+2)%3;
  
  const unsigned int inoD0 = FindAdjEdgeIndex(aTri[itA], inoA2, aTri);
  const unsigned int inoD1 = (inoD0+1)%3;
  
  const unsigned int inoE0 = FindAdjEdgeIndex(aTri[itB], inoB1, aTri);
  const unsigned int inoE1 = (inoE0+1)%3;
  const unsigned int inoE2 = (inoE0+2)%3;
  
  const unsigned int inoF0 = FindAdjEdgeIndex(aTri[itB], inoB2, aTri);
  const unsigned int inoF1 = (inoF0+1)%3;
  
  if (aTri[itC].s2[inoC2]==itD){ // additinoal two triangles to be deleted
    assert( aTri[itD].s2[inoD1]==itC );
    // TODO: implement this collapse
    return false;
  }
  
  if (aTri[itE].s2[inoE2]==itF){ // additinoal two triangles to be deleted
    assert( aTri[itF].s2[inoF1]==itE );
    // TODO: implement this collapse
    return false;
  }
  
  if (itC==itF && itD==itE) return false;
  
  const CDynTri oldA = aTri[itA];
  const CDynTri oldB = aTri[itB];
  
  int ipoW = oldA.v[inoA0];
  int ipoX = oldA.v[inoA1];
  int ipoY = oldB.v[inoB0];
  int ipoZ = oldB.v[inoB1];  // point to be deleted
  
  assert(aTri[itD].v[inoD1]==ipoX);
  assert(aTri[itF].v[inoF1]==ipoZ);
  
  { // check if elements around ipX and elements around ipZ share common triangles
    std::vector<int> ring1;
    { // set triangle index from point 0 to point 1
      int jtri = itF;
      int jnoel_c = inoF1;
      int jnoel_b = (jnoel_c+1)%3;
      for (;;){
        assert(jtri < (int)aTri.size());
        assert(jnoel_c < 3);
        assert(aTri[jtri].v[jnoel_c]==ipoZ);
        ring1.push_back(aTri[jtri].v[(jnoel_c+2)%3]);
        const int ktri = aTri[jtri].s2[jnoel_b];
        if( ktri == -1 ) return false; // boundary
        assert(ktri>=0&&ktri<(int)aTri.size());
        const unsigned int knoel_f = FindAdjEdgeIndex(aTri[jtri], jnoel_b, aTri);
        const unsigned int knoel_c = (knoel_f+1)%3;
        const unsigned int knoel_b = (knoel_f+2)%3;
        assert(aTri[ktri].s2[knoel_f]==jtri);
        assert(aTri[ktri].v[knoel_c]==ipoZ);
        if (ktri==itC) break;
        jtri = ktri;
        jnoel_c = knoel_c;
        jnoel_b = knoel_b;
      }
    }
    std::vector<int> ring2;
    { // set triangle index from point 0 to point 1
      int jtri = itD;
      int jnoel_c = inoD1;
      int jnoel_b = (jnoel_c+1)%3;
      for (;;){
        assert(jtri < (int)aTri.size());
        assert(jnoel_c < 3);
        assert(aTri[jtri].v[jnoel_c]==ipoX);
        ring2.push_back(aTri[jtri].v[(jnoel_c+2)%3]);
        const int ktri = aTri[jtri].s2[jnoel_b];
        if( ktri == -1 ) return false; // boundary
        assert( ktri>=0 && ktri<(int)aTri.size());
        const unsigned int knoel_f = FindAdjEdgeIndex(aTri[jtri], jnoel_b, aTri);
        const unsigned int knoel_c = (knoel_f+1)%3;
        const unsigned int knoel_b = (knoel_f+2)%3;
        assert(aTri[ktri].s2[knoel_f]==jtri);
        assert(aTri[ktri].v[knoel_c]==ipoX);
        if (ktri==itE) break;
        jtri = ktri;
        jnoel_c = knoel_c;
        jnoel_b = knoel_b;
      }
    }
    sort(ring1.begin(), ring1.end());
    sort(ring2.begin(), ring2.end());
    std::vector<int> insc(ring1.size());
    auto it = set_intersection(ring1.begin(), ring1.end(),
                               ring2.begin(), ring2.end(),
                               insc.begin());
    if (it!=insc.begin()){ return false; } // ring1 and ring2 have intersection
  }
 
  assert(oldA.v[inoA1]==oldB.v[inoB2]);
  assert(oldA.v[inoA2]==oldB.v[inoB1]);
  assert(oldA.s2[inoA0]==itB);
  assert(oldB.s2[inoB0]==itA);
  
  // ---------------------------------
  // change from there
  
  aPo[ipoW].e = itC;  aPo[ipoW].d = inoC1;
  aPo[ipoY].e = itE;  aPo[ipoY].d = inoE1;
  aPo[ipoX].e = itD;  aPo[ipoX].d = inoD1;
  aPo[ipoZ].e = -1;
  
  aTri[itC].s2[inoC0] = oldA.s2[inoA2];
  if (oldA.s2[inoA2]>=0&&oldA.s2[inoA2]<(int)aTri.size()){
    assert(oldA.s2[inoA2] < (int)aTri.size());
    aTri[itD].s2[inoD0] = itC;
  }
  
  aTri[itD].s2[inoD0] = oldA.s2[inoA1];
  if (oldA.s2[inoA1]>=0&&oldA.s2[inoA1]<(int)aTri.size()){
    assert(oldA.s2[inoA1] < (int)aTri.size());
    aTri[itC].s2[inoC0] = itD;
  }
  
  aTri[itE].s2[inoE0] = oldB.s2[inoB2];
  if (oldB.s2[inoB2]>=0&&oldB.s2[inoB2] < (int)aTri.size()){
    assert(oldB.s2[inoB2] < (int)aTri.size());
    aTri[itF].s2[inoF0] = itE;
  }
  
  aTri[itF].s2[inoF0] = oldB.s2[inoB1];
  if (oldB.s2[inoB1]>=0&&oldB.s2[inoB1] < (int)aTri.size()){
    assert(oldB.s2[inoB1] < (int)aTri.size());
    aTri[itE].s2[inoE0] = itF;
  }
  
  { // set triangle vtx index from ipoZ to ipoX
    std::vector< std::pair<unsigned int, unsigned int> > aItIn; // itri, inode
    int jtri = itF;
    int jnoel_c = inoF1;
    int jnoel_b = (jnoel_c+1)%3; // noelTriEdge[jnoel_c][0];
    for (;;){
      aItIn.push_back( std::make_pair(jtri,jnoel_c) );
      assert(jtri < (int)aTri.size());
      assert(jnoel_c < 3);
      assert(aTri[jtri].v[jnoel_c]==ipoZ);
      assert(aTri[jtri].s2[jnoel_b]>=0&&aTri[jtri].s2[jnoel_b]<(int)aTri.size());
      int ktri = aTri[jtri].s2[jnoel_b];
      if (ktri==itD||ktri==itE) break;
      const unsigned int knoel_f = FindAdjEdgeIndex(aTri[jtri],jnoel_b,aTri);
      const unsigned int knoel_c = (knoel_f+1)%3;
      const unsigned int knoel_b = (knoel_f+2)%3;
      assert(aTri[ktri].v[knoel_c]==ipoZ);
      assert(itB < (int)aTri.size());
      jtri = ktri;
      jnoel_c = knoel_c;
      jnoel_b = knoel_b;
    }
    for(auto itr=aItIn.begin();itr!=aItIn.end();++itr) {
      unsigned int it0 = itr->first;
      unsigned int in0 = itr->second;
      assert(aTri[it0].v[in0] == ipoZ);
      aTri[it0].v[in0] = ipoX;
    }
  }
  
  {  // isolate two triangles to be deleted
    aTri[itA].s2[0] = -1;  aTri[itA].s2[1] = -1;  aTri[itA].s2[2] = -1;
    aTri[itB].s2[0] = -1;  aTri[itB].s2[1] = -1;  aTri[itB].s2[2] = -1;
    const int itri_1st = (itA > itB) ? itA : itB;
    const int itri_2nd = (itA < itB) ? itA : itB;
    DeleteTri(itri_1st, aPo, aTri);
    DeleteTri(itri_2nd, aPo, aTri);
  }
  return true;
}

DFM2_INLINE void delfem2::GetTriArrayAroundPoint
(std::vector< std::pair<unsigned int,unsigned int> >& aTriSurPo,
 unsigned int ipoin,
 const std::vector<CDynPntSur>& aPo,
 const std::vector<CDynTri>& aTri)
{
  const int itri_ini = aPo[ipoin].e;
  const unsigned int inoel_c_ini = aPo[ipoin].d;
  assert(itri_ini < (int)aTri.size());
  assert(inoel_c_ini < 3);
  assert(aTri[itri_ini].v[inoel_c_ini]==ipoin);
  int itri0 = itri_ini;
  unsigned int inoel_c0 = inoel_c_ini;
  unsigned int inoel_b0 = (inoel_c0+1)%3;
  for (;;){
    assert(itri0 < (int)aTri.size());
    assert(inoel_c0 < 3);
    assert(aTri[itri0].v[inoel_c0]==ipoin);
    aTriSurPo.emplace_back(itri0, inoel_c0);
    // ----------------------
    if (aTri[itri0].s2[inoel_b0]==-1){ break; }
    assert(aTri[itri0].s2[inoel_b0]>=0&&aTri[itri0].s2[inoel_b0] < (int)aTri.size());
    int itri1 = aTri[itri0].s2[inoel_b0];
    const unsigned int inoel_f1 = FindAdjEdgeIndex(aTri[itri0],inoel_b0,aTri);
    const unsigned int inoel_c1 = (inoel_f1+1)%3;
    const unsigned int inoel_b1 = (inoel_f1+2)%3;
    assert(itri1 < (int)aTri.size());
    assert(aTri[itri1].s2[inoel_f1]==(int)itri0);
    assert(aTri[itri1].v[inoel_c1]==ipoin);
    if (itri1==itri_ini) return;
    itri0 = itri1;
    inoel_c0 = inoel_c1;
    inoel_b0 = inoel_b1;
  }
}


DFM2_INLINE void delfem2::extractHoles
(std::vector< std::vector<int> >& aIndP_Hole,
 const int npo,
 const std::vector<CDynTri>& aETri)
{
  aIndP_Hole.clear();
  std::multimap<int,int> mapConnection;
  std::vector<int> aFlg(npo,0);
  for(const auto & itri : aETri){
    for(int inotri=0;inotri<3;++inotri){
      int itris0 = itri.s2[inotri];
      if( itris0 != -1 ) continue;
      const int ip0 = itri.v[(inotri+1)%3];
      const int ip1 = itri.v[(inotri+2)%3];
      mapConnection.insert( std::make_pair(ip0,ip1) );
//      mapConnection.insert( std::make_pair(ip1,ip0) ); // to make the hole ccw
      aFlg[ip0] = 1;
      aFlg[ip1] = 1;
    }
  }
  if( mapConnection.empty() ) return;
  for(int itr=0;itr<npo;++itr){
    int ip_ker0 = -1;
    for(int ipo=0;ipo<npo;++ipo){
      if( aFlg[ipo] == 0 ) continue;
      if( aFlg[ipo] == 1 ){
        ip_ker0 = ipo;
        break;
      }
    }
    if( ip_ker0 == -1 ) break;
    aIndP_Hole.resize(aIndP_Hole.size()+1);
    std::vector<int>& hole = aIndP_Hole[aIndP_Hole.size()-1];
    std::stack<int> stackNext;
    stackNext.push(ip_ker0);
    while(!stackNext.empty()){
      int ip0 = stackNext.top();
      stackNext.pop();
      if( aFlg[ip0] != 1 ) continue;
      aFlg[ip0] = 2;
      hole.push_back(ip0);
      std::pair<std::multimap<int,int>::iterator, std::multimap<int,int>::iterator > its;
      its = mapConnection.equal_range(ip0);
      for(auto it=its.first;it!=its.second;it++){
        assert( it->first == ip0 );
        int ip1 = it->second;
        stackNext.push(ip1);
      }
    }
  }
}

