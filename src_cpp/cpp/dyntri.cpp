/*
 * Copyright (c) 2019 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */


#include <set>
#include <algorithm>
#include <map>
#include <stack>
#include <iostream>
#include <cassert>

#include "delfem2/dyntri.h"

////////////////////////////////////////////////////////////////////////////////////////////////////


void MakeInnerRelationTri
(std::vector<ETri>& aTri, const unsigned int npoin,
 const std::vector<int>& elsup_ind,
 const std::vector<int>& elsup)
{
  const unsigned int EdEd2Rel[3][3] = {
    { 0, 2, 1 },
    { 2, 1, 0 },
    { 1, 0, 2 } };
  
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
            aTri[itri].r2[iedtri] = EdEd2Rel[iedtri][jedtri];
            break;
          }
        }
        if( iflg ) break;
      }
      if( !iflg ){
        aTri[itri].s2[iedtri] = -1;
      }
      for(unsigned int ipofa=0;ipofa<2;ipofa++){
        tmp_poin[ inpofa[ipofa] ] = 0;
      }
    }
  }
}

bool JArray_MakeElSuP
(std::vector<int>& elsup_ind, std::vector<int>& elsup,
 const std::vector<ETri>& aTri, const unsigned int npoin)
{
  const unsigned int nnotri = 3;
  
  elsup_ind.assign(npoin+1, 0);
  for(unsigned int itri=0;itri<aTri.size();itri++){
    for(unsigned int inotri=0;inotri<nnotri;inotri++){
      elsup_ind[ aTri[itri].v[inotri]+1 ]++;
    }
  }
  for(unsigned int ipoin=0;ipoin<npoin;ipoin++){
    elsup_ind[ipoin+1] += elsup_ind[ipoin];
  }
  const int nelsup = elsup_ind[npoin];
  elsup.resize(nelsup);
  for(unsigned int itri=0;itri<aTri.size();itri++){
    for(unsigned int inotri=0;inotri<nnotri;inotri++){
      const unsigned int ipoin0 = aTri[itri].v[inotri];
      const unsigned int ielsup = elsup_ind[ipoin0];
      elsup[ielsup] = itri;
      elsup_ind[ipoin0]++;
    }
  }
  for(int ipoin=npoin;ipoin>0;ipoin--){
    elsup_ind[ipoin] = elsup_ind[ipoin-1];
  }
  elsup_ind[0] = 0;
  return true;
}

void JArray_PSuP
(std::vector<int>& psup_ind, std::vector<int>& psup,
 const std::vector<ETri>& aTri, const unsigned int npoin,
 const std::vector<int>& elsup_ind, const std::vector<int>& elsup)
{
  std::vector<unsigned int> aflg(npoin,0);
  psup_ind[0] = 0;
  for (unsigned int ino = 0; ino<npoin; ino++){
    psup_ind[ino+1] = psup_ind[ino];
    aflg[ino] = ino;
    for (int ielsup = elsup_ind[ino]; ielsup<elsup_ind[ino+1]; ielsup++){
      unsigned int itri1 = elsup[ielsup];
      for (unsigned int inotri = 0; inotri<3; inotri++){
        unsigned int ino1 = aTri[itri1].v[inotri];
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
      for (unsigned int inotri = 0; inotri<3; inotri++){
        unsigned int ino1 = aTri[itri1].v[inotri];
        if (aflg[ino1]==ino) continue;
        psup[iedge] = ino1;
        iedge++;
        aflg[ino1] = ino;
      }
    }
  }
  assert((int)iedge==npsup);
}


////////////////////////////////////////////////////////////////////////////////////////////////////


bool InsertPoint_ElemEdge
(const int ipo_ins,    //the index of the new point
 const int itri_ins,  //triangle index
 const int ied_ins,  //edge index
 std::vector<CEPo2>& aPo,
 std::vector<ETri>& aTri )
{
  assert(itri_ins < (int)aTri.size());
  assert(ipo_ins < (int)aPo.size());
  
  if (aTri[itri_ins].s2[ied_ins]==-1){
    assert(0);
  }
  
  // (node index opp to 0)*3+(node index opp to 1) -> relation index
  const int noel2RelTriTri[9] = {
    -1,  // 0 00
    -1,  // 1 01
    0,  // 2 02
    2, // 3 10
    -1, // 4 11
    -1,  // 5 12
    -1,  // 6 20
    1,  // 7 21
    -1, // 8 22
  };
  
  const int itri_adj = aTri[itri_ins].s2[ied_ins];
  const int ied_adj = (int)relTriTri[(int)aTri[itri_ins].r2[ied_ins]][ied_ins];
  assert(itri_adj < (int)aTri.size());
  assert(ied_ins < 3);
  
  const int itri0 = itri_ins;
  const int itri1 = itri_adj;
  const int itri2 = (int)aTri.size();
  const int itri3 = (int)aTri.size()+1;
  
  aTri.resize(aTri.size()+2);
  
  ETri old0 = aTri[itri_ins];
  ETri old1 = aTri[itri_adj];
  
  const int ino0_0 = ied_ins;
  const int ino1_0 = (ied_ins+1)%3; //noelTriEdge[ied_ins][0];
  const int ino2_0 = (ied_ins+2)%3; //noelTriEdge[ied_ins][1];
  
  const int ino0_1 = ied_adj;
  const int ino1_1 = (ied_adj+1)%3; //noelTriEdge[ied_adj][0];
  const int ino2_1 = (ied_adj+2)%3; //noelTriEdge[ied_adj][1];
  
  assert(old0.v[ino1_0]==old1.v[ino2_1]);
  assert(old0.v[ino2_0]==old1.v[ino1_1]);
  assert(old0.s2[ino0_0]==itri1);
  assert(old1.s2[ino0_1]==itri0);
  
  aPo[ipo_ins].e = itri0;         aPo[ipo_ins].d = 0;
  aPo[old0.v[ino2_0]].e = itri0;  aPo[old0.v[ino2_0]].d = 1;
  aPo[old0.v[ino0_0]].e = itri1;  aPo[old0.v[ino0_0]].d = 1;
  aPo[old1.v[ino2_1]].e = itri2;  aPo[old1.v[ino2_1]].d = 1;
  aPo[old1.v[ino0_1]].e = itri3;  aPo[old1.v[ino0_1]].d = 1;
  
  {
    ETri& ref_tri = aTri[itri0];
    ////////////////
    ref_tri.v[0] = ipo_ins;          ref_tri.v[1] = old0.v[ino2_0];  ref_tri.v[2] = old0.v[ino0_0];
    ref_tri.s2[0] = old0.s2[ino1_0];  ref_tri.s2[1] = itri1;          ref_tri.s2[2] = itri3;
    ////////////////
    if (old0.s2[ino1_0]>=0&&old0.s2[ino1_0]<(int)aTri.size()){
      assert(old0.r2[ino1_0] < 3);
      const unsigned int* rel = relTriTri[old0.r2[ino1_0]];
      ref_tri.r2[0] = noel2RelTriTri[rel[ino1_0]*3+rel[ino2_0]];
      assert(ref_tri.r2[0]>=0&&ref_tri.r2[0] < 3);
      assert(old0.s2[ino1_0] < (int)aTri.size());
      aTri[old0.s2[ino1_0]].s2[rel[ino1_0]] = itri0;
      aTri[old0.s2[ino1_0]].r2[rel[ino1_0]] = invRelTriTri[ref_tri.r2[0]];
    }
    ref_tri.r2[1] = 0;
    ref_tri.r2[2] = 0;
  }
  {
    ETri& ref_tri = aTri[itri1];
    ////////////////
    ref_tri.v[0] = ipo_ins;      ref_tri.v[1] = old0.v[ino0_0];  ref_tri.v[2] = old0.v[ino1_0];
    ref_tri.s2[0] = old0.s2[ino2_0];  ref_tri.s2[1] = itri2;      ref_tri.s2[2] = itri0;
    ////////////////
    if (old0.s2[ino2_0]>=0&&old0.s2[ino2_0]<(int)aTri.size()){
      assert(old0.r2[ino2_0] < 3);
      const unsigned int* rel = relTriTri[old0.r2[ino2_0]];
      ref_tri.r2[0] = noel2RelTriTri[rel[ino2_0]*3+rel[ino0_0]];
      assert(ref_tri.r2[0]>=0&&ref_tri.r2[0] < 3);
      assert(old0.s2[ino2_0] < (int)aTri.size());
      aTri[old0.s2[ino2_0]].s2[rel[ino2_0]] = itri1;
      aTri[old0.s2[ino2_0]].r2[rel[ino2_0]] = invRelTriTri[ref_tri.r2[0]];
    }
    ref_tri.r2[1] = 0;
    ref_tri.r2[2] = 0;
  }
  {
    ETri& ref_tri = aTri[itri2];
    ////////////////
    ref_tri.v[0] = ipo_ins;      ref_tri.v[1] = old1.v[ino2_1];  ref_tri.v[2] = old1.v[ino0_1];
    ref_tri.s2[0] = old1.s2[ino1_1];  ref_tri.s2[1] = itri3;      ref_tri.s2[2] = itri1;
    ////////////////
    if (old1.s2[ino1_1]>=0&&old1.s2[ino1_1]<(int)aTri.size()){
      assert(old1.r2[ino1_1] < 3);
      const unsigned int* rel = relTriTri[old1.r2[ino1_1]];
      ref_tri.r2[0] = noel2RelTriTri[rel[ino1_1]*3+rel[ino2_1]];
      assert(ref_tri.r2[0]>=0&&ref_tri.r2[0] < 3);
      assert(old1.s2[ino1_1] < (int)aTri.size());
      aTri[old1.s2[ino1_1]].s2[rel[ino1_1]] = itri2;
      aTri[old1.s2[ino1_1]].r2[rel[ino1_1]] = invRelTriTri[ref_tri.r2[0]];
    }
    ref_tri.r2[1] = 0;
    ref_tri.r2[2] = 0;
  }
  {
    ETri& ref_tri = aTri[itri3];
    ref_tri.v[0] = ipo_ins;      ref_tri.v[1] = old1.v[ino0_1];  ref_tri.v[2] = old1.v[ino1_1];
    ref_tri.s2[0] = old1.s2[ino2_1];  ref_tri.s2[1] = itri0;      ref_tri.s2[2] = itri2;
    if (old1.s2[ino2_1]>=0&&old1.s2[ino2_1]<(int)aTri.size()){
      assert(old1.r2[ino2_1] < 3);
      const unsigned int* rel = relTriTri[old1.r2[ino2_1]];
      ref_tri.r2[0] = noel2RelTriTri[rel[ino2_1]*3+rel[ino0_1]];
      assert(ref_tri.r2[0]>=0&&ref_tri.r2[0] < 3);
      assert(old1.s2[ino2_1] < (int)aTri.size());
      aTri[old1.s2[ino2_1]].s2[rel[ino2_1]] = itri3;
      aTri[old1.s2[ino2_1]].r2[rel[ino2_1]] = invRelTriTri[ref_tri.r2[0]];
    }
    ref_tri.r2[1] = 0;
    ref_tri.r2[2] = 0;
  }
  return true;
}



bool InsertPoint_Elem
(const int ipo_ins,
 const int itri_ins,
 std::vector<CEPo2>& aPo,
 std::vector<ETri>& aTri)
{
  assert(itri_ins < (int)aTri.size());
  assert(ipo_ins < (int)aPo.size());
  
  // (node index opp to 0)*3+(node index opp to 1) -> relation index
  const int noel2RelTriTri[9] = {
    -1,  // 0 00
    -1,  // 1 01
    0,  // 2 02
    2, // 3 10
    -1, // 4 11
    -1,  // 5 12
    -1,  // 6 20
    1,  // 7 21
    -1, // 8 22
  };
  
  const int itri0 = itri_ins;
  const int itri1 = (int)aTri.size();
  const int itri2 = (int)aTri.size()+1;
  
  aTri.resize(aTri.size()+2);
  const ETri old_tri = aTri[itri_ins];
  
  aPo[ipo_ins].e = itri0;        aPo[ipo_ins].d = 0;
  aPo[old_tri.v[0]].e = itri1;  aPo[old_tri.v[0]].d = 2;
  aPo[old_tri.v[1]].e = itri2;  aPo[old_tri.v[1]].d = 2;
  aPo[old_tri.v[2]].e = itri0;  aPo[old_tri.v[2]].d = 2;
  
  {
    ETri& ref_tri = aTri[itri0];
    ref_tri.v[0] = ipo_ins;      ref_tri.v[1] = old_tri.v[1];  ref_tri.v[2] = old_tri.v[2];
    ref_tri.s2[0] = old_tri.s2[0];  ref_tri.s2[1] = itri1;      ref_tri.s2[2] = itri2;
    if (old_tri.s2[0]>=0&&old_tri.s2[0]<(int)aTri.size()){
      assert(old_tri.r2[0] < 3);
      const unsigned int* rel = &relTriTri[old_tri.r2[0]][0];
      ref_tri.r2[0] = noel2RelTriTri[rel[0]*3+rel[1]];
      assert(ref_tri.r2[0]>=0&&ref_tri.r2[0] < 3);
      assert(old_tri.s2[0] < (int)aTri.size());
      aTri[old_tri.s2[0]].s2[rel[0]] = itri0;
      aTri[old_tri.s2[0]].r2[rel[0]] = invRelTriTri[ref_tri.r2[0]];
    }
    ref_tri.r2[1] = 0;
    ref_tri.r2[2] = 0;
  }
  {
    ETri& ref_tri = aTri[itri1];
    ref_tri.v[0] = ipo_ins;      ref_tri.v[1] = old_tri.v[2];  ref_tri.v[2] = old_tri.v[0];
    ref_tri.s2[0] = old_tri.s2[1];  ref_tri.s2[1] = itri2;      ref_tri.s2[2] = itri0;
    if (old_tri.s2[1]>=0&&old_tri.s2[1]<(int)aTri.size()){
      assert(old_tri.r2[1] < 3);
      const unsigned int* rel = &relTriTri[old_tri.r2[1]][0];
      ref_tri.r2[0] = noel2RelTriTri[rel[1]*3+rel[2]];
      assert(ref_tri.r2[0]>=0&&ref_tri.r2[0] < 3);
      assert(old_tri.s2[1] < (int)aTri.size());
      aTri[old_tri.s2[1]].s2[rel[1]] = itri1;
      aTri[old_tri.s2[1]].r2[rel[1]] = invRelTriTri[ref_tri.r2[0]];
    }
    ref_tri.r2[1] = 0;
    ref_tri.r2[2] = 0;
  }
  {
    ETri& ref_tri = aTri[itri2];
    ref_tri.v[0] = ipo_ins;     ref_tri.v[1] = old_tri.v[0];  ref_tri.v[2] = old_tri.v[1];
    ref_tri.s2[0] = old_tri.s2[2];  ref_tri.s2[1] = itri0;      ref_tri.s2[2] = itri1;
    if (old_tri.s2[2]>=0&&old_tri.s2[2]<(int)aTri.size()){
      assert(old_tri.r2[2] < 3);
      const unsigned int* rel = &relTriTri[old_tri.r2[2]][0];
      ref_tri.r2[0] = noel2RelTriTri[rel[2]*3+rel[0]];
      assert(ref_tri.r2[0]>=0&&ref_tri.r2[0] < 3);
      assert(old_tri.s2[2] < (int)aTri.size());
      aTri[old_tri.s2[2]].s2[rel[2]] = itri2;
      aTri[old_tri.s2[2]].r2[rel[2]] = invRelTriTri[ref_tri.r2[0]];
    }
    ref_tri.r2[1] = 0;
    ref_tri.r2[2] = 0;
  }
  return true;
}


bool FlipEdge
(int itri0, int ied0,
 std::vector<CEPo2>& aPo,
 std::vector<ETri>& aTri)
{
  assert(itri0 < (int)aTri.size());
  assert(ied0 < 3);
  assert(aTri[itri0].s2[ied0]>=0&&aTri[itri0].s2[ied0]<(int)aTri.size());
  
  const int itri1 = aTri[itri0].s2[ied0];
  const int ied1 = relTriTri[aTri[itri0].r2[ied0]][ied0];
  assert(itri1 < (int)aTri.size());
  assert(ied1 < 3);
  assert(aTri[itri1].s2[ied1]>=0&&aTri[itri0].s2[ied0]<(int)aTri.size());
  
  // (node index opp to 0)*3+(node index opp to 1) -> relation index
  const int noel2RelTriTri[9] = {
    -1,  // 0 00
    -1,  // 1 01
    0,  // 2 02
    2, // 3 10
    -1, // 4 11
    -1,  // 5 12
    -1,  // 6 20
    1,  // 7 21
    -1, // 8 22
  };
  
  //  std::cout << itri0 << "-" << ied0 << "    " << itri1 << "-" << ied1 << std::endl;
  
  ETri old0 = aTri[itri0];
  ETri old1 = aTri[itri1];
  
  const int no0_0 = ied0;
  const int no1_0 = (ied0+1)%3; //noelTriEdge[ied0][0];
  const int no2_0 = (ied0+2)%3; //noelTriEdge[ied0][1];
  
  const int no0_1 = ied1;
  const int no1_1 = (ied1+1)%3; //noelTriEdge[ied1][0];
  const int no2_1 = (ied1+2)%3; //noelTriEdge[ied1][1];
  
  if(old0.s2[no1_0]==old1.s2[no2_1]) return false; // inverted mesh
  if(old0.s2[no2_0]==old1.s2[no1_1]) return false; // inverted mesh
  
  assert(old0.v[no1_0]==old1.v[no2_1]);
  assert(old0.v[no2_0]==old1.v[no1_1]);
  
  aPo[old0.v[no1_0]].e = itri0;  aPo[old0.v[no1_0]].d = 0;
  aPo[old0.v[no0_0]].e = itri0;  aPo[old0.v[no0_0]].d = 2;
  aPo[old1.v[no1_1]].e = itri1;  aPo[old1.v[no1_1]].d = 0;
  aPo[old1.v[no0_1]].e = itri1;  aPo[old1.v[no0_1]].d = 2;
  
  {
    ETri& ref_tri = aTri[itri0];
    ////////////////
    ref_tri.v[0] = old0.v[no1_0];  ref_tri.v[1] = old1.v[no0_1];  ref_tri.v[2] = old0.v[no0_0];
    ref_tri.s2[0] = itri1;      ref_tri.s2[1] = old0.s2[no2_0];  ref_tri.s2[2] = old1.s2[no1_1];
    ////////////////
    ref_tri.r2[0] = 0;
    if (old0.s2[no2_0]>=0){
      assert(old0.r2[no2_0] < 3);
      const unsigned int* rel = relTriTri[old0.r2[no2_0]];
      assert(old0.s2[no2_0] < (int)aTri.size());
      assert(old0.s2[no2_0]!=itri0);
      assert(old0.s2[no2_0]!=itri1);
      ref_tri.r2[1] = noel2RelTriTri[rel[no1_0]*3+rel[no2_0]];
      assert(ref_tri.r2[1]>=0&&ref_tri.r2[1] < 3);
      aTri[old0.s2[no2_0]].s2[rel[no2_0]] = itri0;
      aTri[old0.s2[no2_0]].r2[rel[no2_0]] = invRelTriTri[ref_tri.r2[1]];
    }
    if (old1.s2[no1_1]>=0){
      assert(old1.r2[no1_1] < 3);
      const unsigned int* rel = relTriTri[old1.r2[no1_1]];
      assert(old1.s2[no1_1] < (int)aTri.size());
      ref_tri.r2[2] = noel2RelTriTri[rel[no2_1]*3+rel[no0_1]];
      assert(ref_tri.r2[2]>=0&&ref_tri.r2[2] < 3);
      aTri[old1.s2[no1_1]].s2[rel[no1_1]] = itri0;
      aTri[old1.s2[no1_1]].r2[rel[no1_1]] = invRelTriTri[ref_tri.r2[2]];
    }
  }
  
  {
    ETri& ref_tri = aTri[itri1];
    ////////////////
    ref_tri.v[0] = old1.v[no1_1];  ref_tri.v[1] = old0.v[no0_0];  ref_tri.v[2] = old1.v[no0_1];
    ref_tri.s2[0] = itri0;      ref_tri.s2[1] = old1.s2[no2_1];  ref_tri.s2[2] = old0.s2[no1_0];
    ////////////////
    ref_tri.r2[0] = 0;
    if (old1.s2[no2_1]>=0){
      assert(old1.r2[no2_1] < 3);
      const unsigned int* rel = relTriTri[old1.r2[no2_1]];
      assert(old1.s2[no2_1] < (int)aTri.size());
      ref_tri.r2[1] = noel2RelTriTri[rel[no1_1]*3+rel[no2_1]];
      assert(ref_tri.r2[1]>=0&&ref_tri.r2[1] < 3);
      aTri[old1.s2[no2_1]].s2[rel[no2_1]] = itri1;
      aTri[old1.s2[no2_1]].r2[rel[no2_1]] = invRelTriTri[ref_tri.r2[1]];
    }
    if (old0.s2[no1_0]>=0){
      assert(old0.r2[no1_0] < 3);
      const unsigned int* rel = relTriTri[old0.r2[no1_0]];
      assert(old0.s2[no1_0] < (int)aTri.size());
      ref_tri.r2[2] = noel2RelTriTri[rel[no2_0]*3+rel[no0_0]];
      assert(ref_tri.r2[2]>=0&&ref_tri.r2[2] < 3);
      aTri[old0.s2[no1_0]].s2[rel[no1_0]] = itri1;
      aTri[old0.s2[no1_0]].r2[rel[no1_0]] = invRelTriTri[ref_tri.r2[2]];
    }
  }
  
  return true;
}

bool FindEdge_LookAroundPoint
(int& itri0, int& inotri0, int& inotri1,
 ///
 const int ipo0, const int ipo1,
 const std::vector<CEPo2>& aPo,
 const std::vector<ETri>& aTri)
{
  const int itri_ini = aPo[ipo0].e;
  const int inotri_ini = aPo[ipo0].d;
  int inotri_cur = inotri_ini;
  int itri_cur = itri_ini;
  for (;;){  // serch clock-wise
    assert(aTri[itri_cur].v[inotri_cur]==ipo0);
    {
      const int inotri2 = (inotri_cur+1)%3; // indexRot3[1][inotri_cur];
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
      const int inotri2 = (inotri_cur+2)%3; // indexRot3[2][inotri_cur];
      if (aTri[itri_cur].s2[inotri2]==-1){ break; }
      const int itri_nex = aTri[itri_cur].s2[inotri2];
      const unsigned int* rel = relTriTri[aTri[itri_cur].r2[inotri2]];
      const int inotri3 = rel[inotri_cur];
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
      const int inotri2 = (inotri_cur+1)%3; // indexRot3[1][inotri_cur];
      if (aTri[itri_cur].s2[inotri2]==-1){ break; }
      const int itri_nex = aTri[itri_cur].s2[inotri2];
      const unsigned int* rel = relTriTri[aTri[itri_cur].r2[inotri2]];
      const int inotri3 = rel[inotri_cur];
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
      const int inotri2 = (inotri_cur+1)%3; // indexRot3[1][inotri_cur];
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



bool FindEdge_LookAllTriangles
(int& itri0, int& iedtri0,
 ///
 const int ipo0, const int ipo1,
 const std::vector<ETri>& aTri)
{
  for(unsigned int itri=0;itri<aTri.size();++itri){
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

bool CheckTri( const std::vector<ETri>& aTri )
{
	const int ntri = (int)aTri.size();
	for(int itri=0;itri<ntri;itri++){
		const ETri& ref_tri = aTri[itri];
    /*		for(int inotri=0;inotri<nNoTri;inotri++){
     assert( ref_tri.v[inotri] >= 0 );
     }*/
		for(int iedtri=0;iedtri<3;iedtri++){
			if( ref_tri.s2[iedtri] >=0 && ref_tri.s2[iedtri]<(int)aTri.size() ){
				const int itri_s = ref_tri.s2[iedtri];
				const int irel = ref_tri.r2[iedtri];
				assert( itri_s < ntri );
				assert( irel < 3 );
				// check sorounding
				{
					const int noel_dia = relTriTri[irel][iedtri];
					assert( aTri[itri_s].s2[noel_dia] == itri );
          //					std::cout << itri << " " << itri_s << std::endl;
				}
				// check relation 
				for(int inoed=0;inoed<2;inoed++){
					const int inoel = (iedtri+1+inoed)%3;
          if( ref_tri.v[inoel] != aTri[itri_s].v[ (int)relTriTri[irel][inoel] ] ){
						std::cout << itri << " " << iedtri << " " << itri_s << " " << ref_tri.v[inoel] << " " << aTri[itri_s].v[ (int)relTriTri[irel][inoel] ] << std::endl;
					}
          assert( ref_tri.v[inoel] == aTri[itri_s].v[ (int)relTriTri[irel][inoel] ] );
				}
			}
		}
	}
  
	return true;
}

bool CheckTri
(const std::vector<CEPo2>& aPo3D,
 const std::vector<ETri>& aSTri,
 bool is_assert)
{
  const int npo = (int)aPo3D.size();
  const int ntri = (int)aSTri.size();
  
  for (int itri = 0; itri<ntri; itri++){
    const ETri& tri0 = aSTri[itri];
    if( tri0.v[0] == -1 ){
      assert(tri0.v[1] == -1);
      assert(tri0.v[2] == -1);
      continue;
    }
    if (tri0.v[0]==tri0.v[1]){ //     assert(tri0.v[2]!=tri0.v[0]);
      if( is_assert ){ assert(0); }
      return false;
    }
    if (tri0.v[1]==tri0.v[2]){ //     assert(tri0.v[2]!=tri0.v[0]);
      if( is_assert ){ assert(0); }
      return false;
    }
    if (tri0.v[2]==tri0.v[0]){ //     assert(tri0.v[2]!=tri0.v[0]);
      if( is_assert ){ assert(0); }
      return false;
    }
    ////
    if( (tri0.s2[0]==tri0.s2[1]) && tri0.s2[0] >= 0 ){ // assert(tri0.s2[0]!=tri0.s2[2]);
      //      std::cout << tri0.s2[0] << " " << tri0.s2[1] << std::endl;
      if( is_assert ){ assert(0); }
      return false;
    }
    if( (tri0.s2[1]==tri0.s2[2]) && tri0.s2[1] >= 0 ){ // assert(tri0.s2[1]!=tri0.s2[2]);
      if( is_assert ){ assert(0); }
      return false;
    }
    if( (tri0.s2[2]==tri0.s2[0]) && tri0.s2[0] >= 0 ){ // assert(tri0.s2[2]!=tri0.s2[0]);
      if( is_assert ){ assert(0); }
      return false;
    }
    /////
    for (unsigned int inotri = 0; inotri<3; inotri++){
      assert(tri0.v[inotri] < npo);
    }
    for (int iedtri = 0; iedtri<3; iedtri++){
      if (tri0.s2[iedtri]>=0&&tri0.s2[iedtri]<ntri){
        const int itri_s = tri0.s2[iedtri];
        const int irel = tri0.r2[iedtri];
        assert(itri_s < ntri);
        assert(irel < 3);
        {
          const int noel_dia = relTriTri[irel][iedtri];
          assert(noel_dia < 3);
          if (aSTri[itri_s].s2[noel_dia]!=itri){ // neibough relation broken
            //            std::cout << itri << " " << itri_s << " " << noel_dia << std::endl;
          }
          assert(aSTri[itri_s].s2[noel_dia]==itri);
        }
        // check relation
        for (int inoed = 0; inoed<2; inoed++){
          const int inoel = (iedtri+1+inoed)%3;//noelTriEdge[iedtri][inoed];
          if (tri0.v[inoel]!=aSTri[itri_s].v[(int)relTriTri[irel][inoel]]){
          }
          assert(tri0.v[inoel]==aSTri[itri_s].v[(int)relTriTri[irel][inoel]]);
        }
      }
    }
  }
  for (int ipoin = 0; ipoin<npo; ++ipoin){
    const int itri0 = aPo3D[ipoin].e;
    const int inoel0 = aPo3D[ipoin].d;
    if (aPo3D[ipoin].e>=0){
      assert(aPo3D[ipoin].d>=0&&aPo3D[ipoin].d < 3);
      if (aSTri[itri0].v[inoel0]!=ipoin){}
      assert(aSTri[itri0].v[inoel0]==ipoin);
    }
  }
  return true;
}



void InitializeMesh
(std::vector<CEPo2>& aPo3D,
 std::vector<ETri>& aSTri,
 ////
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


//////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////

void MoveCCW
(int& itri_cur,
 int& inotri_cur,
 bool& flag_is_wall,
 ////
 std::vector<ETri>& aTri)
{
  const int inotri1 = (inotri_cur+1)%3; // indexRot3[1][inotri_cur];
  if (aTri[itri_cur].s2[inotri1]==-1){ flag_is_wall = true; return; }
  flag_is_wall = false;
  const int itri_nex = aTri[itri_cur].s2[inotri1];
  const unsigned int* rel_nex = relTriTri[aTri[itri_cur].r2[inotri1]];
  const int inotri_nex = rel_nex[inotri_cur];
  itri_cur = itri_nex;
  inotri_cur = inotri_nex;
}


bool DeleteTri
(int itri_to,
 std::vector<CEPo2>& aPo,
 std::vector<ETri>& aTri)
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
    assert(itri_adj < aTri.size());
    const unsigned int* rel = relTriTri[(int)aTri[itri_to].r2[iedtri]];
    const unsigned int iedtri_adj = rel[iedtri];
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


bool Collapse_ElemEdge
(const int itri_del,
 const int ied_del,
 std::vector<CEPo2>& aPo,
 std::vector<ETri>& aTri)
{
  assert( itri_del >= 0 && itri_del < (int)aTri.size() );
  assert( ied_del >= 0 && ied_del < 3 );
  if (aTri[itri_del].s2[ied_del]==-1){
    std::cout<<"Error!-->Not Implemented: Mesh with hole"<<std::endl;
    assert(0);
  }
  
  const int itri_adj = aTri[itri_del].s2[ied_del];
  const int ied_adj = (int)relTriTri[(int)aTri[itri_del].r2[ied_del]][ied_del];
  assert(itri_adj < (int)aTri.size());
  assert(ied_adj < 3);
  assert(aTri[itri_adj].s2[ied_adj]==itri_del);
  
  const int itri0 = itri_del;
  const int itri1 = itri_adj;
  const int itri2 = aTri[itri0].s2[(ied_del+1)%3];
  const int itri3 = aTri[itri0].s2[(ied_del+2)%3];
  const int itri4 = aTri[itri1].s2[(ied_adj+1)%3];
  const int itri5 = aTri[itri1].s2[(ied_adj+2)%3];
  if( itri2 == -1 ) return false;
  if( itri3 == -1 ) return false;
  if( itri4 == -1 ) return false;
  if( itri5 == -1 ) return false;
  
  const int ino0_0 = ied_del;
  const int ino1_0 = (ino0_0+1)%3;
  const int ino2_0 = (ino0_0+2)%3;
  
  const int ino0_1 = ied_adj;
  const int ino1_1 = (ino0_1+1)%3; // point to be deleted
  const int ino2_1 = (ino0_1+2)%3;
  
  const int ino0_2 = relTriTri[(int)aTri[itri0].r2[ino1_0]][ino1_0];
  const int ino1_2 = (ino0_2+1)%3;
  const int ino2_2 = (ino0_2+2)%3;
  
  const int ino0_3 = relTriTri[(int)aTri[itri0].r2[ino2_0]][ino2_0];
  const int ino1_3 = (ino0_3+1)%3;
  
  const int ino0_4 = relTriTri[(int)aTri[itri1].r2[ino1_1]][ino1_1];
  const int ino1_4 = (ino0_4+1)%3;
  const int ino2_4 = (ino0_4+2)%3;
  
  const int ino0_5 = relTriTri[(int)aTri[itri1].r2[ino2_1]][ino2_1];
  const int ino1_5 = (ino0_5+1)%3;
  
  if (aTri[itri2].s2[ino2_2]==itri3) return false;
  if (aTri[itri3].s2[ino1_3]==itri2) return false;
  if (aTri[itri4].s2[ino2_4]==itri5) return false;
  if (aTri[itri5].s2[ino1_5]==itri4) return false;
  if (itri2==itri5 && itri3==itri4) return false;
  
  const ETri old0 = aTri[itri0];
  const ETri old1 = aTri[itri1];
  
  int ipo0 = old0.v[ino0_0];
  int ipo1 = old0.v[ino1_0];
  int ipo2 = old1.v[ino0_1];
  int ipo3 = old1.v[ino1_1];  // point to be deleted
  
  assert(aTri[itri3].v[ino1_3]==ipo1);
  assert(aTri[itri5].v[ino1_5]==ipo3);
  
  {
    std::vector<int> ring1;
    { // set triangle index from point 0 to point 1
      int jtri = itri5;
      int jnoel_c = ino1_5;
      int jnoel_b = (jnoel_c+1)%3; // noelTriEdge[jnoel_c][0];
      for (;;){
        assert(jtri < (int)aTri.size());
        assert(jnoel_c < 3);
        assert(aTri[jtri].v[jnoel_c]==ipo3);
        ring1.push_back(aTri[jtri].v[(jnoel_c+2)%3]);
        const int ktri = aTri[jtri].s2[jnoel_b];
        if( ktri == -1 ) return false;
        assert(ktri>=0&&ktri<(int)aTri.size());
        const int rel01 = aTri[jtri].r2[jnoel_b];
        const int knoel_c = relTriTri[rel01][jnoel_c];
        const int knoel_b = relTriTri[rel01][(jnoel_c+2)%3];
        assert(aTri[ktri].s2[relTriTri[rel01][jnoel_b]]==jtri);
        assert(aTri[ktri].v[knoel_c]==ipo3);
        if (ktri==itri2) break;
        jtri = ktri;
        jnoel_c = knoel_c;
        jnoel_b = knoel_b;
      }
    }
    std::vector<int> ring2;
    { // set triangle index from point 0 to point 1
      int jtri = itri3;
      int jnoel_c = ino1_3;
      int jnoel_b = (jnoel_c+1)%3; // noelTriEdge[jnoel_c][0];
      for (;;){
        assert(jtri < (int)aTri.size());
        assert(jnoel_c < 3);
        assert(aTri[jtri].v[jnoel_c]==ipo1);
        ring2.push_back(aTri[jtri].v[(jnoel_c+2)%3]);
        const int ktri = aTri[jtri].s2[jnoel_b];
        if( ktri == -1 ) return false;
        assert( ktri>=0 && ktri<(int)aTri.size());
        const int rel01 = aTri[jtri].r2[jnoel_b];
        const int knoel_c = relTriTri[rel01][jnoel_c];
        const int knoel_b = relTriTri[rel01][(jnoel_c+2)%3];
        assert(aTri[ktri].s2[relTriTri[rel01][jnoel_b]]==jtri);
        assert(aTri[ktri].v[knoel_c]==ipo1);
        if (ktri==itri4) break;
        jtri = ktri;
        jnoel_c = knoel_c;
        jnoel_b = knoel_b;
      }
    }
    sort(ring1.begin(), ring1.end());
    sort(ring2.begin(), ring2.end());
    std::vector<int> insc(ring1.size());
    std::vector<int>::iterator it = set_intersection(ring1.begin(), ring1.end(),
                                                     ring2.begin(), ring2.end(),
                                                     insc.begin());
    if (it!=insc.begin()){ return  false; }
  }
  
  /*
   std::cout<<std::endl;
   std::cout<<"stt"<<std::endl;
   std::cout<<"tris : "<<itri0<<" "<<itri1<<" "<<itri2<<" "<<itri3<<" "<<itri4<<" "<<itri5<<std::endl;
   std::cout<<"vtxs : "<<ipo0<<" "<<ipo1<<" "<<ipo2<<" "<<ipo3<<std::endl;
   */
  
  assert(old0.v[ino1_0]==old1.v[ino2_1]);
  assert(old0.v[ino2_0]==old1.v[ino1_1]);
  assert(old0.s2[ino0_0]==itri1);
  assert(old1.s2[ino0_1]==itri0);
  
  {
    aPo[ipo0].e = itri2;  aPo[ipo0].d = ino1_2;
    aPo[ipo2].e = itri4;  aPo[ipo2].d = ino1_4;
    aPo[ipo1].e = itri3;  aPo[ipo1].d = ino1_3;
    aPo[ipo3].e = -1;
  }
  
  // (edge index)*3+(opp edge index) -> relation index
  const int ed2RelTriTri[9] = {
    0,  // 0 00
    2,  // 1 01
    1,  // 2 02
    2,  // 3 10
    1,  // 4 11
    0,  // 5 12
    1,  // 6 20
    0,  // 7 21
    2,  // 8 22
  };
  
  { // change itri2
    ETri& tri = aTri[itri2];
    //    tri.g2[ino0_2] = old0.g2[ino2_0];
    tri.s2[ino0_2] = old0.s2[ino2_0];
    if (old0.s2[ino2_0]>=0&&old0.s2[ino2_0]<(int)aTri.size()){
      assert(old0.r2[ino2_0] < 3);
      assert(old0.s2[ino2_0] < (int)aTri.size());
      tri.r2[ino0_2] = ed2RelTriTri[ino0_2*3+ino0_3];
      assert(tri.r2[ino0_2]>=0&&tri.r2[ino0_2] < 3);
      aTri[itri3].s2[ino0_3] = itri2;
      aTri[itri3].r2[ino0_3] = invRelTriTri[tri.r2[ino0_2]];
    }
  }
  
  { // change itri3
    ETri& tri = aTri[itri3];
    //    tri.g2[ino0_3] = old0.g2[ino1_0];
    tri.s2[ino0_3] = old0.s2[ino1_0];
    if (old0.s2[ino1_0]>=0&&old0.s2[ino1_0]<(int)aTri.size()){
      //    if( old0.g2[ino1_0] == -2 || old0.g2[ino1_0] == -3 ){
      assert(old0.r2[ino1_0] < 3);
      assert(old0.s2[ino1_0] < (int)aTri.size());
      tri.r2[ino0_3] = ed2RelTriTri[ino0_3*3+ino0_2];
      assert(tri.r2[ino0_3]>=0&&tri.r2[ino0_3] < 3);
      aTri[itri2].s2[ino0_2] = itri3;
      aTri[itri2].r2[ino0_2] = invRelTriTri[tri.r2[ino0_3]];
    }
  }
  
  { // change itri4
    ETri& tri = aTri[itri4];
    //    tri.g2[ino0_4] = old1.g2[ino2_1];
    tri.s2[ino0_4] = old1.s2[ino2_1];
    //    if( old1.g2[ino2_1] == -2 || old1.g2[ino2_1] == -3 ){
    if (old1.s2[ino2_1]>=0&&old1.s2[ino2_1] < (int)aTri.size()){
      assert(old1.r2[ino2_1] < 3);
      assert(old1.s2[ino2_1] < (int)aTri.size());
      tri.r2[ino0_4] = ed2RelTriTri[ino0_4*3+ino0_5];
      assert(tri.r2[ino0_4]>=0&&tri.r2[ino0_4] < 3);
      aTri[itri5].s2[ino0_5] = itri4;
      aTri[itri5].r2[ino0_5] = invRelTriTri[tri.r2[ino0_4]];
      assert((int)relTriTri[(int)aTri[itri4].r2[ino0_4]][ino0_4]==ino0_5);
      assert((int)relTriTri[(int)aTri[itri5].r2[ino0_5]][ino0_5]==ino0_4);
    }
  }
  
  { // change itri5
    ETri& tri = aTri[itri5];
    //    tri.g2[ino0_5] = old1.g2[ino1_1];
    tri.s2[ino0_5] = old1.s2[ino1_1];
    //    if( old1.g2[ino1_1] == -2 || old1.g2[ino1_1] == -3 ){
    if (old1.s2[ino1_1]>=0&&old1.s2[ino1_1] < (int)aTri.size()){
      assert(old1.r2[ino1_1] < 3);
      assert(old1.s2[ino1_1] < (int)aTri.size());
      tri.r2[ino0_5] = ed2RelTriTri[ino0_5*3+ino0_4];
      assert(tri.r2[ino0_5]>=0&&tri.r2[ino0_5] < 3);
      aTri[itri4].s2[ino0_4] = itri5;
      aTri[itri4].r2[ino0_4] = invRelTriTri[tri.r2[ino0_5]];
      assert((int)relTriTri[aTri[itri5].r2[ino0_5]][ino0_5]==ino0_4);
      assert((int)relTriTri[aTri[itri4].r2[ino0_4]][ino0_4]==ino0_5);
    }
  }
  
  { // set triangle index from point 0 to point 1
    int jtri = itri5;
    int jnoel_c = ino1_5;
    int jnoel_b = (jnoel_c+1)%3; // noelTriEdge[jnoel_c][0];
    for (;;){
      assert(jtri < (int)aTri.size());
      assert(jnoel_c < 3);
      assert(aTri[jtri].v[jnoel_c]==ipo3);
      //      std::cout << " fan : " << jtri << "     " << aTri[jtri].v[0] << " " << aTri[jtri].v[1] << " " << aTri[jtri].v[2] << std::endl;
      aTri[jtri].v[jnoel_c] = ipo1;
      //      assert( aTri[jtri].g2[jnoel_b] == -2 );
      assert(aTri[jtri].s2[jnoel_b]>=0&&aTri[jtri].s2[jnoel_b]<(int)aTri.size());
      int ktri = aTri[jtri].s2[jnoel_b];
      const int rel01 = aTri[jtri].r2[jnoel_b];
      const int knoel_c = relTriTri[rel01][jnoel_c];
      const int knoel_b = relTriTri[rel01][(jnoel_c+2)%3];
      assert(itri1 < (int)aTri.size());
      assert(aTri[ktri].s2[relTriTri[rel01][jnoel_b]]==jtri);
      if (ktri==itri3||ktri==itri4) break;
      jtri = ktri;
      jnoel_c = knoel_c;
      jnoel_b = knoel_b;
    }
  }
  
  {  // isolate two triangles to be deleted
    aTri[itri0].s2[0] = -1;  aTri[itri0].s2[1] = -1;  aTri[itri0].s2[2] = -1;
    aTri[itri1].s2[0] = -1;  aTri[itri1].s2[1] = -1;  aTri[itri1].s2[2] = -1;
    const int itri_1st = (itri0 > itri1) ? itri0 : itri1;
    const int itri_2nd = (itri0 < itri1) ? itri0 : itri1;
    DeleteTri(itri_1st, aPo, aTri);
    DeleteTri(itri_2nd, aPo, aTri);
  }
  return true;
}

void GetTriArrayAroundPoint
(std::vector< std::pair<int,int> >& aTriSurPo,
 int ipoin,
 const std::vector<CEPo2>& aPo,
 const std::vector<ETri>& aTri)
{
  const int itri_ini = aPo[ipoin].e;
  const int inoel_c_ini = aPo[ipoin].d;
  assert(itri_ini < (int)aTri.size());
  assert(inoel_c_ini < 3);
  assert(aTri[itri_ini].v[inoel_c_ini]==ipoin);
  int itri0 = itri_ini;
  int inoel_c0 = inoel_c_ini;
  int inoel_b0 = (inoel_c0+1)%3;
  for (;;){
    assert(itri0 < (int)aTri.size());
    assert(inoel_c0 < 3);
    assert(aTri[itri0].v[inoel_c0]==ipoin);
    aTriSurPo.push_back(std::make_pair(itri0, inoel_c0));
    /////
    if (aTri[itri0].s2[inoel_b0]==-1){ break; }
    assert(aTri[itri0].s2[inoel_b0]>=0&&aTri[itri0].s2[inoel_b0] < (int)aTri.size());
    int itri1 = aTri[itri0].s2[inoel_b0];
    const int rel01 = aTri[itri0].r2[inoel_b0];
    const int inoel_c1 = relTriTri[rel01][inoel_c0];
    const int inoel_b1 = relTriTri[rel01][(inoel_c0+2)%3];
    assert(itri1 < (int)aTri.size());
    assert(aTri[itri1].s2[relTriTri[rel01][inoel_b0]]==(int)itri0);
    assert(aTri[itri1].v[inoel_c1]==ipoin);
    if (itri1==itri_ini) return;
    itri0 = itri1;
    inoel_c0 = inoel_c1;
    inoel_b0 = inoel_b1;
  }
}






void extractHoles
(std::vector< std::vector<int> >& aIndP_Hole,
 const int npo,
 const std::vector<ETri>& aETri)
{
  aIndP_Hole.clear();
  std::multimap<int,int> mapConnection;
  std::vector<int> aFlg(npo,0);
  for(int itri=0;itri<(int)aETri.size();itri++){
    for(int inotri=0;inotri<3;++inotri){
      int itris0 = aETri[itri].s2[inotri];
      if( itris0 != -1 ) continue;
      const int ip0 = aETri[itri].v[(inotri+1)%3];
      const int ip1 = aETri[itri].v[(inotri+2)%3];
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
      for(std::multimap<int,int>::iterator it=its.first;it!=its.second;it++){
        assert( it->first == ip0 );
        int ip1 = it->second;
        stackNext.push(ip1);
      }
    }
  }
}

