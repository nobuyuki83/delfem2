/*
 * Copyright (c) 2019 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */


#ifndef DFM2_SLICE_H
#define DFM2_SLICE_H

#include <stdio.h>
#include <set>
#include <stack>
#include <vector>

#include "delfem2/dfm2_inline.h"

namespace delfem2 {
  
class CSegInfo{
public:
  void Initialize(
      int jtri0,
      const unsigned int* aTri,
      size_t nTri,
      const double* aLevelVtx,
      double height);
  void Pos3D(
      double pA[3],
      double pB[3],
      const std::vector<double>& aXYZ,
      const std::vector<unsigned int>& aTri) const ;
  void Pos2D(
      double pA[2],
      double pB[2],
      const double* aXY,
      const unsigned int* aTri) const ;
public:
  int itri;
  int iedA, iedB;
  double r0A, r0B;
};

void AddContour(
    std::vector<CSegInfo>& aSeg,
    double thres,
    const unsigned int* aTri,
    size_t nTri,
    const double* aVal);

// ---------------------------------------
  
class CSliceTriMesh
{
  public:
    CSliceTriMesh(unsigned int ih): iHeight(ih) {}
    
    int IndHeight() const { return iHeight; } // requirement for MakeReebGraph
    size_t NumSeg() const { return aTriInfo.size(); } // requirement for MakeReebGraph
    unsigned int IndTri_Seg(unsigned int iseg) const { return aTriInfo[iseg].itri; } // requirement for MakeReebGraph
    
  public:
  public:
    unsigned int iHeight;
    std::vector<CSegInfo> aTriInfo;
};

void Slice_MeshTri3D_Heights(
    std::vector<CSliceTriMesh>& aCS,
    //
    const std::vector<double>& aHeight,
    const std::vector<double>& aHeightVtx,
    const std::vector<unsigned int>& aTri,
    const std::vector<unsigned int>& aTriSuTri);

// T must have following functions
//     int IndHeight() const;
//     unsigned int NumSeg() const;
//     unsigned int IndTri_Seg(unsigned int iseg) const;
template <typename T>
void MakeReebGraph
(std::vector< std::set<unsigned int> >& aConnectivity,
 const std::vector<T>& aCrossSection,
 const std::vector<unsigned int>& aTri,
 const std::vector<unsigned int>& aTriSuTri)
{
  const unsigned int ntri = (unsigned int)aTri.size()/3;
  const size_t nCS = aCrossSection.size();
  aConnectivity.clear();
  aConnectivity.resize(nCS);
  std::vector<int> Tri2HCS;
  Tri2HCS.resize(ntri,-1);
  for(unsigned int ics=0;ics<nCS;ics++){
    const T& cs0 = aCrossSection[ics];
    for(unsigned int iseg=0;iseg<cs0.NumSeg();iseg++){
      unsigned int itri0 = cs0.IndTri_Seg(iseg);
      if( Tri2HCS[itri0] != -1 ){
        unsigned int jcs0 = Tri2HCS[itri0];
        const T& cs1 = aCrossSection[jcs0];
        if( abs( cs0.IndHeight() - cs1.IndHeight()) == 1 ){
          aConnectivity[ics].insert(jcs0);
          aConnectivity[jcs0].insert(ics);
        }
      }
      else{
        Tri2HCS[itri0] = ics;
      }
    }
  }
  for(unsigned int ics=0;ics<nCS;ics++){
    const T& cs0 = aCrossSection[ics];
    for(unsigned int iseg=0;iseg<cs0.NumSeg();iseg++){
      unsigned int itri0 = cs0.IndTri_Seg(iseg);
      for(unsigned int iedtri=0;iedtri<3;iedtri++){
        unsigned int jtri0 = aTriSuTri[itri0*3+iedtri];
        if( Tri2HCS[jtri0] == -1 ){ continue; }
        unsigned int jhcs = Tri2HCS[jtri0];
        const T& cs1 = aCrossSection[jhcs];
        if( abs( cs0.IndHeight() - cs1.IndHeight() ) != 1 ) continue;
        aConnectivity[ics].insert(jhcs);
        aConnectivity[jhcs].insert(ics);
      }
    }
  }
  unsigned int itri_ker = 0;
  for(;;){
    std::stack<unsigned int> stackTri;
    for(;itri_ker<ntri;itri_ker++){
      if( Tri2HCS[itri_ker] == -1 ){ stackTri.push(itri_ker); break; }
    }
    if( stackTri.empty() )break;
    std::set<unsigned int> setAdjCS;
    std::vector<int> aFlg; aFlg.resize(nCS,-1);
    for(;;){
      if( stackTri.empty() ) break;
      unsigned int jtri0 = stackTri.top();
      stackTri.pop();
      if( Tri2HCS[jtri0] != -1 ) continue;
      Tri2HCS[jtri0] = -2;
      for(unsigned int jedtri=0;jedtri<3;jedtri++){
        unsigned int ktri0 = aTriSuTri[jtri0*3+jedtri];
        if(      Tri2HCS[ktri0] == -2 ) continue; // already studied
        else if( Tri2HCS[ktri0] == -1 ){ stackTri.push(ktri0); }
        else{
          unsigned int ihcs = Tri2HCS[ktri0];
          if( aFlg[ihcs] != -1 ){ continue; }
          else{
            setAdjCS.insert(ihcs);
            aFlg[ihcs] = 1;
          }
        }
      }
    }
    unsigned int nAdjCS = (unsigned int)setAdjCS.size();
    std::vector<unsigned int> aAdjCS;
    aAdjCS.reserve(nAdjCS);
    int is_single_height = -1;
    for(std::set<unsigned int>::iterator itr=setAdjCS.begin();itr!=setAdjCS.end();itr++){
      unsigned int ics = *itr;
      if( is_single_height == -1 ){
        is_single_height = aCrossSection[ics].IndHeight();
      }
      else if( is_single_height != aCrossSection[ics].IndHeight() ){
        is_single_height = -2;
      }
      aAdjCS.push_back(ics);
    }
    if( is_single_height == -2 ){ is_single_height = 0; }
    else{                         is_single_height = 1; }
    for(unsigned int iac=0;iac<nAdjCS;iac++){
      unsigned int ihcs = aAdjCS[iac];
      const T& cs0 = aCrossSection[ihcs];
      for(unsigned int jac=0;jac<nAdjCS;jac++){
        if( iac == jac ) continue;
        unsigned int jhcs = aAdjCS[jac];
        const T& cs1 = aCrossSection[jhcs];
        if( !is_single_height && cs0.IndHeight() == cs1.IndHeight() ){ continue; }
        if( abs( cs0.IndHeight() - cs1.IndHeight() ) > 1 ) continue;
        aConnectivity[ihcs].insert(jhcs);
        aConnectivity[jhcs].insert(ihcs);
      }
    }
  }
}

} // namespace delfem2

#ifndef DFM2_STATIC_LIBRARY
# include "delfem2/slice.cpp"
#endif

#endif /* SLICE_H */
