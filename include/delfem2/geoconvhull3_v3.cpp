/*
 * Copyright (c) 2019 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */

#include <cmath>
#include <stack>
#include "delfem2/geoconvhull3_v3.h"

#ifndef M_PI
#  define M_PI 3.14159265358979323846
#endif

namespace delfem2 {
namespace convhull3 {

template <typename T>
bool IsOut(
    int itri,
    const CVec3 <T> &v,
    const std::vector<CVec3 < T>> & aXYZ,
    const std::vector<int> &aTri)
{
    int i0 = aTri[itri * 3 + 0];
    int i1 = aTri[itri * 3 + 1];
    int i2 = aTri[itri * 3 + 2];
    const CVec3 <T> &v0 = aXYZ[i0];
    const CVec3 <T> &v1 = aXYZ[i1];
    const CVec3 <T> &v2 = aXYZ[i2];
    CVec3 <T> n;
    Normal(n, v0, v1, v2);
    double dot = Dot(v - v0, n);
    return dot > 0;
}

//! Volume of a tetrahedra
template <typename T>
T Volume_Tet(
    const CVec3<T>& v0,
    const CVec3<T>& v1,
    const CVec3<T>& v2,
    const CVec3<T>& v3 )
{
//  return delfem2::Volume_Tet3(v0.p, v1.p, v2.p, v3.p);
  double v = (v1.p[0]-v0.p[0])*( (v2.p[1]-v0.p[1])*(v3.p[2]-v0.p[2]) - (v3.p[1]-v0.p[1])*(v2.p[2]-v0.p[2]) )
             + (v1.p[1]-v0.p[1])*( (v2.p[2]-v0.p[2])*(v3.p[0]-v0.p[0]) - (v3.p[2]-v0.p[2])*(v2.p[0]-v0.p[0]) )
             + (v1.p[2]-v0.p[2])*( (v2.p[0]-v0.p[0])*(v3.p[1]-v0.p[1]) - (v3.p[0]-v0.p[0])*(v2.p[1]-v0.p[1]) );
  return v*0.16666666666666666666666666666667;
}

}
}

// -------------------------------------

template <typename T>
void delfem2::ConvexHull(
    std::vector<int>& aTri,
    const std::vector<CVec3<T>>& aXYZ)
{
  std::vector<int> aBflg( aXYZ.size(), -1 );
  aTri.reserve(aXYZ.size()*6);
  aTri.resize(4*3);
  aTri[ 0] = 1;  aTri[ 1] = 2;  aTri[ 2] = 3; // 0
  aTri[ 3] = 0;  aTri[ 4] = 3;  aTri[ 5] = 2; // 1
  aTri[ 6] = 0;  aTri[ 7] = 1;  aTri[ 8] = 3; // 2
  aTri[ 9] = 0;  aTri[10] = 2;  aTri[11] = 1; // 3
  std::vector< std::pair<int,int> > aTriSur;
  aTriSur.resize(4*3);
  aTriSur[ 0] = std::make_pair(1,0);
  aTriSur[ 1] = std::make_pair(2,0);
  aTriSur[ 2] = std::make_pair(3,0);
  aTriSur[ 3] = std::make_pair(0,0);
  aTriSur[ 4] = std::make_pair(3,2);
  aTriSur[ 5] = std::make_pair(2,1);
  aTriSur[ 6] = std::make_pair(0,1);
  aTriSur[ 7] = std::make_pair(1,2);
  aTriSur[ 8] = std::make_pair(3,1);
  aTriSur[ 9] = std::make_pair(0,2);
  aTriSur[10] = std::make_pair(2,2);
  aTriSur[11] = std::make_pair(1,1);
  {
    double vol = delfem2::convhull3::Volume_Tet(
        aXYZ[0],
        aXYZ[1],
        aXYZ[2],
        aXYZ[3]);
    if( vol < 0 ){
      aTri[ 0] = 3;  aTri[ 1] = 2;  aTri[ 2] = 1; // 0
      aTri[ 3] = 2;  aTri[ 4] = 3;  aTri[ 5] = 0; // 1
      aTri[ 6] = 3;  aTri[ 7] = 1;  aTri[ 8] = 0; // 2
      aTri[ 9] = 1;  aTri[10] = 2;  aTri[11] = 0; // 3
      aTriSur[ 0] = std::make_pair(3,2);
      aTriSur[ 1] = std::make_pair(2,2);
      aTriSur[ 2] = std::make_pair(1,2);
      aTriSur[ 3] = std::make_pair(2,1);
      aTriSur[ 4] = std::make_pair(3,0);
      aTriSur[ 5] = std::make_pair(0,2);
      aTriSur[ 6] = std::make_pair(3,1);
      aTriSur[ 7] = std::make_pair(1,0);
      aTriSur[ 8] = std::make_pair(0,1);
      aTriSur[ 9] = std::make_pair(1,1);
      aTriSur[10] = std::make_pair(2,0);
      aTriSur[11] = std::make_pair(0,0);
    }
  }
  const int triEd[3][2] = {
    { 1, 2 },
    { 2, 0 },
    { 0, 1 } };
  for(int iv=4;iv<(int)aXYZ.size();iv++){
    CVec3<T> v = aXYZ[iv];
    int itri_ker = -1;
    for(int itri=0;itri<(int)aTri.size()/3;itri++){
      if( delfem2::convhull3::IsOut(itri,v,aXYZ,aTri) ){ itri_ker = itri; break; }
    }
#ifndef NDEBUG
    {
      for(unsigned int itri=0;itri<aTri.size()/3;itri++){
        for(int ied=0;ied<3;ied++){
          int ied1 = triEd[ied][0];
          int ied2 = triEd[ied][1];
          int itri_s = aTriSur[itri*3+ied].first;
          int ied_s0 = aTriSur[itri*3+ied].second;
          assert( aTriSur[itri_s*3+ied_s0].first  == (int)itri );
          assert( aTriSur[itri_s*3+ied_s0].second == ied );
          int ied_s1 = triEd[ied_s0][0];
          int ied_s2 = triEd[ied_s0][1];
          assert( aTri[itri*3+ied1] == aTri[itri_s*3+ied_s2] );
          assert( aTri[itri*3+ied2] == aTri[itri_s*3+ied_s1] );
        }
      }
    }
#endif
    if( itri_ker == -1 ) continue; // inside
    std::vector< std::pair<int,int> > aB;
    std::vector<int> isDelTri( aTri.size()/3, -1 );
    {
      std::vector<int> isLookedEdge( aTri.size(), -1 );
      std::stack< std::pair<int,int> > sBound;
      { // initialize
        sBound.push( aTriSur[itri_ker*3+0] );
        sBound.push( aTriSur[itri_ker*3+1] );
        sBound.push( aTriSur[itri_ker*3+2] );
        isDelTri[itri_ker] = 1;
      }
      for(;;){
        if( sBound.empty() ) break;
        int itri0 = sBound.top().first;
        int ied0  = sBound.top().second;
        sBound.pop();
        if( isLookedEdge[itri0*3+ied0] == 1 ) continue;
        isLookedEdge[itri0*3+ied0] = 1;
        {
          const std::pair<int,int>& s0 = aTriSur[itri0*3+ied0];
          isLookedEdge[s0.first*3+s0.second] = 1;
        }
        isDelTri[itri0] = ( ::delfem2::convhull3::IsOut(itri0,v,aXYZ,aTri) ) ? 1 : 0;
        if( isDelTri[itri0] == 1 ){ // expand this boundary
          int ied1 = triEd[ied0][0];
          int ied2 = triEd[ied0][1];
          const std::pair<int,int>& s1 = aTriSur[itri0*3+ied1];
          const std::pair<int,int>& s2 = aTriSur[itri0*3+ied2];
          assert( aTriSur[s1.first*3+s1.second].first == itri0 );
          assert( aTriSur[s2.first*3+s2.second].first == itri0 );
          sBound.push( s1 );
          sBound.push( s2 );
        }
        else{ // this is a actuall boundary
          aB.emplace_back(itri0,ied0 );
        }
      }
    }
    std::vector<int> aBSur( aB.size()*2, -1);
    {
      for(auto & ib : aB){
        int itri0 = ib.first;
        int itn0  = ib.second;
        int itn1 = triEd[itn0][0];
        int itn2 = triEd[itn0][1];
        int iv1 = aTri[itri0*3+itn1];
        int iv2 = aTri[itri0*3+itn2];
        aBflg[iv1] = -1;
        aBflg[iv2] = -1;
      }
      for(int ib=0;ib<(int)aB.size();ib++){
        int itri0 = aB[ib].first;
        int itn0  = aB[ib].second;
        int itn1 = triEd[itn0][0];
        int itn2 = triEd[itn0][1];
        int iv1 = aTri[itri0*3+itn1];
        int iv2 = aTri[itri0*3+itn2];
        if(      aBflg[iv1] == -2 ){}
        else if( aBflg[iv1] == -1 ){ aBflg[iv1] = ib; }
        else{
          assert( aBflg[iv1] >= 0 );
          int ib0 = aBflg[iv1];
          aBSur[ib *2+1] = ib0;
          aBSur[ib0*2+0] = ib;
          aBflg[iv1] = -2;
        }
        if(      aBflg[iv2] == -2 ){}
        else if( aBflg[iv2] == -1 ){ aBflg[iv2] = ib; }
        else{
          assert( aBflg[iv2] >= 0 );
          int ib0 = aBflg[iv2];
          aBSur[ib *2+0] = ib0;
          aBSur[ib0*2+1] = ib;
          aBflg[iv2] = -2;
        }
      }
    }
#ifndef NDEBUG
    for(std::size_t ib=0;ib<aB.size();ib++){
      for(int inb=0;inb<2;inb++){
        int itri0 = aB[ib].first;
        int itn0  = aB[ib].second;
        int iv1 = aTri[itri0*3+triEd[itn0][0]];
        int iv2 = aTri[itri0*3+triEd[itn0][1]];
        int ib_s0 = aBSur[ib*2+inb];
        int itri_s0 = aB[ib_s0].first;
        int itn_s0  = aB[ib_s0].second;
        int iv_s1 = aTri[itri_s0*3+triEd[itn_s0][0]];
        int iv_s2 = aTri[itri_s0*3+triEd[itn_s0][1]];
        if( inb == 0 ){ assert( iv2 == iv_s1 ); }
        else{           assert( iv1 == iv_s2 ); }
      }
    }
#endif
    std::vector<int> mapOld2New( aTri.size()/3, -1 );
    std::vector<int> aTri1; aTri1.reserve(aTri.size()+60);
    std::vector< std::pair<int,int> > aTriSur1; aTriSur1.reserve(aTriSur.size()+60);
    for(int itri=0;itri<(int)aTri.size()/3;itri++){
      if( isDelTri[itri] ==  1) continue;
      assert( !::delfem2::convhull3::IsOut(itri,v,aXYZ,aTri) );
      // itri is inside
      mapOld2New[itri] = (int)aTri1.size()/3;
      aTri1.push_back( aTri[itri*3+0] );
      aTri1.push_back( aTri[itri*3+1] );
      aTri1.push_back( aTri[itri*3+2] );
      aTriSur1.emplace_back(-1,0 );
      aTriSur1.emplace_back(-1,0 );
      aTriSur1.emplace_back(-1,0 );
    }
    for(int itri=0;itri<(int)aTri.size()/3;itri++){ // set old relation
      if( isDelTri[itri] ==  1) continue;
      int jtri0 = mapOld2New[itri];
      assert( jtri0 >= 0 && (int)aTri1.size()/3 );
      for(int iet=0;iet<3;iet++){
        int itri_s = aTriSur[itri*3+iet].first;
        if( mapOld2New[itri_s] == -1 ) continue;
        aTriSur1[jtri0*3+iet].first = mapOld2New[itri_s];
        aTriSur1[jtri0*3+iet].second = aTriSur[itri*3+iet].second;
      }
    }
    const int ntri_old = (int)aTri1.size()/3;
    for(int ib=0;ib<(int)aB.size();ib++){
      int itri0 = aB[ib].first;
      int itn0  = aB[ib].second;
      int itn1 = triEd[itn0][0];
      int itn2 = triEd[itn0][1];
      assert( !::delfem2::convhull3::IsOut(itri0,v,aXYZ,aTri) );
#ifndef NDEBUG
      {
        int itri_s = aTriSur[itri0*3+itn0].first;
        assert( ::delfem2::convhull3::IsOut(itri_s,v,aXYZ,aTri) );
        int ied_s0 = aTriSur[itri0*3+itn0].second;
        assert( aTriSur[itri_s*3+ied_s0].first == itri0 );
        assert( aTriSur[itri_s*3+ied_s0].second == itn0 );
        int ied_s1 = triEd[ied_s0][0];
        int ied_s2 = triEd[ied_s0][1];
        assert( aTri[itri0*3+itn1] == aTri[itri_s*3+ied_s2] );
        assert( aTri[itri0*3+itn2] == aTri[itri_s*3+ied_s1] );
      }
#endif
      assert( isDelTri[itri0] == 0 );
      int jtri0 = mapOld2New[itri0]; assert( jtri0 != -1 );
      int jtri1 = (int)aTri1.size()/3;
      assert( jtri1 == ntri_old + ib );
      aTri1.push_back( iv );
      aTri1.push_back( aTri[itri0*3+itn2] );
      aTri1.push_back( aTri[itri0*3+itn1] );
      aTriSur1[jtri0*3+itn0] = std::make_pair(jtri1,0);
      ////
      int jtri2 = aBSur[ib*2+0] + ntri_old;
      int jtri3 = aBSur[ib*2+1] + ntri_old;
      aTriSur1.emplace_back(jtri0,itn0 );
      aTriSur1.emplace_back(jtri3,2 );
      aTriSur1.emplace_back(jtri2,1 );
    }
    aTri    = aTri1;
    aTriSur = aTriSur1;
  }
}
#ifndef DFM2_HEADER_ONLY
template void delfem2::ConvexHull(
    std::vector<int>& aTri,
    const std::vector<CVec3d>& aXYZ);
#endif

