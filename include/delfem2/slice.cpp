/*
 * Copyright (c) 2019 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */

#include <cassert>
#include <climits>

#include "delfem2/slice.h"

// ---------------------------------------------

namespace delfem2 {
namespace slice {

DFM2_INLINE void IndexElement_OverlapLevels_MeshTri3D
(std::vector< std::vector<unsigned int> >& aCST,
 //
 const std::vector<double>& aH,
 const std::vector<double>& aLevelVtx,
 const std::vector<unsigned int>& aTri)
{
  const size_t ntri = aTri.size()/3;
  const size_t nH = aH.size();
  aCST.resize(nH);
  for(unsigned int itri=0;itri<ntri;itri++){
    const double ah[3] = {
      aLevelVtx[ aTri[itri*3+0] ],
      aLevelVtx[ aTri[itri*3+1] ],
      aLevelVtx[ aTri[itri*3+2] ],
    };
    for(unsigned int ih=0;ih<nH;ih++){
      unsigned int icnt = 0;
      for(double inotri : ah){
        if( inotri-aH[ih] > 0 ){ icnt++; }
      }
      if( icnt == 1 || icnt == 2 ){
        aCST[ih].push_back(itri);
      }
    }
  }
}

DFM2_INLINE bool TraverseBoundaryLoop(
	CSliceTriMesh& cs,
	std::vector<int>& aFlgSeg,
	int iseg_ker, 
	[[maybe_unused]] int ih,
	const std::vector<int>& Tri2Seg,
	const std::vector<unsigned int>& aCST,
	double height,
	const std::vector<double>& aLevelVtx,
	const std::vector<unsigned int>& aTri,
	const std::vector<unsigned int>& aTriSuTri)
{
  cs.aTriInfo.clear();
  unsigned int iseg_next = iseg_ker;
  for(;;){ // seg in cs loop
    assert( aFlgSeg[iseg_next] == 0 );
    int jtri0 = aCST[iseg_next];
    aFlgSeg[iseg_next] = 1;
    CSegInfo info;
    info.Initialize(
        jtri0,
        aTri.data(),aTri.size()/3,
        aLevelVtx.data(),height);
    cs.aTriInfo.push_back(info);
    unsigned int iedge_next = info.iedB;
    int itri_next1 = aTriSuTri[jtri0*3+iedge_next];
    if( itri_next1 == -1 ){ break; } // open loop discard
    int iseg_next1 = Tri2Seg[itri_next1];
    assert( iseg_next1 != -1 );
    if( iseg_next1 == iseg_ker ){ return true; } // is_open == false
    iseg_next = iseg_next1;
  }
  // reach here if the loop is open
  cs.aTriInfo.clear();
  //
  iseg_next = iseg_ker;
  for(;;){ // seg in cs loop
    int jtri0 = aCST[iseg_next];
    unsigned int iflg = 0;
    double aH[3];
    for(unsigned int inotri=0;inotri<3;inotri++){
      aH[inotri] = aLevelVtx[ aTri[jtri0*3+inotri] ] -height;
      if( aH[inotri] < 0 ){ continue; }
      if( inotri == 0 ){ iflg += 1; }
      if( inotri == 1 ){ iflg += 2; }
      if( inotri == 2 ){ iflg += 4; }
    }
    int iedge_next = -1;
    if( iflg == 1 ){ iedge_next = 2; } // 0
    if( iflg == 2 ){ iedge_next = 0; } // 1
    if( iflg == 4 ){ iedge_next = 1; } // 2
    if( iflg == 3 ){ iedge_next = 0; } // 01
    if( iflg == 5 ){ iedge_next = 2; } // 02
    if( iflg == 6 ){ iedge_next = 1; } // 12
    assert( iedge_next != -1 );
    unsigned int itri_next1 = aTriSuTri[jtri0*3+iedge_next];
    if( itri_next1 == UINT_MAX ){ break; } // reached boundary
    const int iseg_next1 = Tri2Seg[itri_next1];
    assert( iseg_next1 != -1 );
    if( iseg_next1 == iseg_ker ) break;
    iseg_next = iseg_next1;
    assert( aFlgSeg[iseg_next] == 0 );
    jtri0 = aCST[iseg_next];
    aFlgSeg[iseg_next] = 1;
  }
  return false;
}

}
}

// ----------------------------------------------

void delfem2::CSegInfo::Initialize(
    int jtri0,
    const unsigned int* aTri,
    [[maybe_unused]] size_t nTri,
    const double* aLevelVtx,
    double height)
{
  this->itri = jtri0;
  unsigned int iflg = 0;
  double aH[3];
  for(unsigned int inotri=0;inotri<3;inotri++){
    aH[inotri] = aLevelVtx[ aTri[jtri0*3+inotri] ] - height;
    if( aH[inotri] < 0 ){ continue; }
    if( inotri == 0 ){ iflg += 1; }
    if( inotri == 1 ){ iflg += 2; }
    if( inotri == 2 ){ iflg += 4; }
  }
  if( iflg == 1 ){
    this->iedA = 2; this->r0A = +aH[0]/(aH[0]-aH[1]);
    this->iedB = 1; this->r0B = -aH[2]/(aH[0]-aH[2]);
  }
  if( iflg == 2 ){
    this->iedA = 0; this->r0A = +aH[1]/(aH[1]-aH[2]);
    this->iedB = 2; this->r0B = -aH[0]/(aH[1]-aH[0]);
  }
  if( iflg == 4 ){
    this->iedA = 1; this->r0A = +aH[2]/(aH[2]-aH[0]);
    this->iedB = 0; this->r0B = -aH[1]/(aH[2]-aH[1]);
  }
  if( iflg == 3 ){
    this->iedA = 0; this->r0A = -aH[1]/(aH[2]-aH[1]);
    this->iedB = 1; this->r0B = +aH[2]/(aH[2]-aH[0]);
  }
  if( iflg == 5 ){
    this->iedA = 2; this->r0A = -aH[0]/(aH[1]-aH[0]);
    this->iedB = 0; this->r0B = +aH[1]/(aH[1]-aH[2]);
  }
  if( iflg == 6 ){
    this->iedA = 1; this->r0A = -aH[2]/(aH[0]-aH[2]);
    this->iedB = 2; this->r0B = +aH[0]/(aH[0]-aH[1]);
  }
}

void delfem2::CSegInfo::Pos3D(
    double pA[3],
    double pB[3],
    //
    const std::vector<double>& aXYZ,
    const std::vector<unsigned int>& aTri) const
{
  const unsigned int i0 = aTri[this->itri*3+0];
  const unsigned int i1 = aTri[this->itri*3+1];
  const unsigned int i2 = aTri[this->itri*3+2];
  const double aP[3][3] = {
    { aXYZ[i0*3+0], aXYZ[i0*3+1], aXYZ[i0*3+2] },
    { aXYZ[i1*3+0], aXYZ[i1*3+1], aXYZ[i1*3+2] },
    { aXYZ[i2*3+0], aXYZ[i2*3+1], aXYZ[i2*3+2] } };
  pA[0] = (1.0-this->r0A)*aP[(this->iedA+1)%3][0] + (this->r0A)*aP[(this->iedA+2)%3][0];
  pA[1] = (1.0-this->r0A)*aP[(this->iedA+1)%3][1] + (this->r0A)*aP[(this->iedA+2)%3][1];
  pA[2] = (1.0-this->r0A)*aP[(this->iedA+1)%3][2] + (this->r0A)*aP[(this->iedA+2)%3][2];
  pB[0] = (1.0-this->r0B)*aP[(this->iedB+1)%3][0] + (this->r0B)*aP[(this->iedB+2)%3][0];
  pB[1] = (1.0-this->r0B)*aP[(this->iedB+1)%3][1] + (this->r0B)*aP[(this->iedB+2)%3][1];
  pB[2] = (1.0-this->r0B)*aP[(this->iedB+1)%3][2] + (this->r0B)*aP[(this->iedB+2)%3][2];
}

void delfem2::CSegInfo::Pos2D(
    double pA[2],
    double pB[2],
    //
    const double* aXY,
    const unsigned int* aTri) const
{
  const unsigned int i0 = aTri[this->itri*3+0];
  const unsigned int i1 = aTri[this->itri*3+1];
  const unsigned int i2 = aTri[this->itri*3+2];
  const double aP[3][2] = {
    { aXY[i0*2+0], aXY[i0*2+1] },
    { aXY[i1*2+0], aXY[i1*2+1] },
    { aXY[i2*2+0], aXY[i2*2+1] } };
  pA[0] = (1.0-this->r0A)*aP[(this->iedA+1)%3][0] + (this->r0A)*aP[(this->iedA+2)%3][0];
  pA[1] = (1.0-this->r0A)*aP[(this->iedA+1)%3][1] + (this->r0A)*aP[(this->iedA+2)%3][1];
  pB[0] = (1.0-this->r0B)*aP[(this->iedB+1)%3][0] + (this->r0B)*aP[(this->iedB+2)%3][0];
  pB[1] = (1.0-this->r0B)*aP[(this->iedB+1)%3][1] + (this->r0B)*aP[(this->iedB+2)%3][1];
}

void delfem2::AddContour(
    std::vector<CSegInfo>& aSeg,
    double thres,
    const unsigned int* vec_tri,
    size_t num_tri,
    const double* value_tri)
{
  for(unsigned int itri=0;itri<num_tri;++itri){
    const double v0 = value_tri[ vec_tri[itri*3+0] ];
    const double v1 = value_tri[ vec_tri[itri*3+1] ];
    const double v2 = value_tri[ vec_tri[itri*3+2] ];
    if(   (v0-thres)*(v1-thres) >= 0
       && (v1-thres)*(v2-thres) >= 0
       && (v2-thres)*(v0-thres) >= 0 ){
      continue;
    }
    delfem2::CSegInfo seg;
    seg.Initialize(itri,
                   vec_tri,num_tri,value_tri,
                   thres);
    aSeg.push_back(seg);
  }
}

// ----------------------------------------------




void delfem2::Slice_MeshTri3D_Heights
(std::vector<CSliceTriMesh>& aCS,
 ////
 const std::vector<double>& aLevel,
 const std::vector<double>& aLevelVtx,
 const std::vector<unsigned int>& aTri,
 const std::vector<unsigned int>& aTriSuTri)
{
  const std::size_t ntri = aTri.size()/3;
  const std::size_t nH = aLevel.size();
  //
  std::vector< std::vector<unsigned int> > aCST;
  slice::IndexElement_OverlapLevels_MeshTri3D(aCST,
                                              aLevel,
                                              aLevelVtx,
                                              aTri);
  aCS.clear();
  std::vector<int> Tri2Seg;
  Tri2Seg.resize(ntri,-1);
  for(unsigned int ih=0;ih<nH;ih++){ // h loop
    for(unsigned int isg=0;isg<aCST[ih].size();isg++){
      unsigned int itri = aCST[ih][isg];
      Tri2Seg[itri] = isg;
    }
    std::vector<int> aFlgSeg;
    aFlgSeg.resize(aCST[ih].size(),0);
    unsigned int iseg_ker = 0;
    for(;;){ // cs loop
      for(;iseg_ker<aCST[ih].size();iseg_ker++){
        if( aFlgSeg[iseg_ker] == 0 ){ break; }
      }
      if( iseg_ker == aCST[ih].size() ) break;
      CSliceTriMesh cs(ih);
      const bool is_closed = slice::TraverseBoundaryLoop(cs, aFlgSeg,
                                                         iseg_ker, ih, Tri2Seg,
                                                         aCST[ih], aLevel[ih],
                                                         aLevelVtx,
                                                         aTri, aTriSuTri);
      if( !is_closed ){ continue; }
      aCS.push_back(cs);
    }
    //
    for(unsigned int itri : aCST[ih]){
      Tri2Seg[itri] = -1;
    }
  }
}


/*
void Slice_MeshTri2D_Contour
 (std::vector<CSliceTriMesh>& aCS,
  //
  const std::vector<double>& aLevel,
  const std::vector<double>& aXY,
  const std::vector<double>& aVal,
  const std::vector<unsigned int>& aTri,
  const std::vector<int>& aTriSur)
{
  const unsigned int ntri = (unsigned int)aTri.size()/3;
  const unsigned int nH = (unsigned int)aLevel.size();
    //
  std::vector< std::vector<unsigned int> > aCST;
  IndexElement_OverlapLevels_MeshTri3D(aCST,
                                       aLevel,norm,origin,
                                       aXYZ,aTri);
    //
  aCS.clear();
  std::vector<int> Tri2Seg;
  Tri2Seg.resize(ntri,-1);
  for(unsigned int ih=0;ih<nH;ih++){ // h loop
    for(unsigned int isg=0;isg<aCST[ih].size();isg++){
      unsigned int itri = aCST[ih][isg];
      Tri2Seg[itri] = isg;
    }
      //
    std::vector<int> aFlgSeg;
    aFlgSeg.resize(aCST[ih].size(),0);
    unsigned int iseg_ker = 0;
    for(;;){ // cs loop
      for(;iseg_ker<aCST[ih].size();iseg_ker++){
        if( aFlgSeg[iseg_ker] == 0 ){ break; }
      }
      if( iseg_ker == aCST[ih].size() ) break;
        //
      CSliceTriMesh cs(ih);
      const bool is_closed = TraverseBoundaryLoop(cs, aFlgSeg,
                                                  iseg_ker, ih, Tri2Seg,
                                                  aCST[ih], norm, origin, aLevel[ih],
                                                  aXYZ, aTri, aTriSur);
      if( !is_closed ){ continue; }
      aCS.push_back(cs);
    }
      //
    for(unsigned int isg=0;isg<aCST[ih].size();isg++){
      unsigned int itri = aCST[ih][isg];
      Tri2Seg[itri] = -1;
    }
  }
}
 */
