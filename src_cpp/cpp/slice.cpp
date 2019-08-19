/*
 * Copyright (c) 2019 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */


#include "delfem2/slice.h"




void IndexElement_OverlapLevels_MeshTri3D
(std::vector< std::vector<unsigned int> >& aCST,
 ////
 const std::vector<double>& aH,
 const double norm[3],
 const double centerOfGravity[3],
 const std::vector<double>& aXYZ,
 const std::vector<unsigned int>& aTri)
{
  const unsigned int ntri = aTri.size()/3;
  const unsigned int nH = aH.size();
  aCST.resize(nH);
  for(unsigned int itri=0;itri<ntri;itri++){
    double ah[3];
    for(unsigned int inotri=0;inotri<3;inotri++){
      unsigned int ino0 = aTri[itri*3+inotri];
      const double p0[3] = {
        aXYZ[ino0*3+0]-centerOfGravity[0],
        aXYZ[ino0*3+1]-centerOfGravity[1],
        aXYZ[ino0*3+2]-centerOfGravity[2] };
      ah[inotri] = p0[0]*norm[0] +  p0[1]*norm[1] + p0[2]*norm[2];
    }
    for(unsigned int ih=0;ih<nH;ih++){
      unsigned int icnt = 0;
      for(unsigned int inotri=0;inotri<3;inotri++){
        if( ah[inotri]-aH[ih] > 0 ){ icnt++; }
      }
      ////
      if( icnt == 1 || icnt == 2 ){
        aCST[ih].push_back(itri);
      }
    }
  }
}


bool TraverseBoundaryLoop
(CSliceTriMesh& cs,
 std::vector<int>& aFlgSeg,
 int iseg_ker, int ih,
 const std::vector<int>& Tri2Seg,
 const std::vector<unsigned int>& aCST,
 const double norm[3],
 const double origin[3],
 double height,
 const std::vector<double>& aXYZ,
 const std::vector<unsigned int>& aTri,
 const std::vector<int>& aTriSur)
{
  cs.aTriInfo.clear();
  unsigned int iseg_next = iseg_ker;
  for(;;){ // seg in cs loop
    assert( aFlgSeg[iseg_next] == 0 );
    int jtri0 = aCST[iseg_next];
    aFlgSeg[iseg_next] = 1;
    CSliceTriMesh::CSegInfo info;
    {
      info.itri = jtri0;
      unsigned int iflg = 0;
      double aP[3][3], aH[3];
      for(unsigned int inotri=0;inotri<3;inotri++){
        unsigned int jno0 = aTri[jtri0*3+inotri];
        aP[inotri][0] = aXYZ[jno0*3+0];
        aP[inotri][1] = aXYZ[jno0*3+1];
        aP[inotri][2] = aXYZ[jno0*3+2];
        aH[inotri]
        = (aP[inotri][0]-origin[0])*norm[0]
        + (aP[inotri][1]-origin[1])*norm[1]
        + (aP[inotri][2]-origin[2])*norm[2] -height;
        if( aH[inotri] < 0 ){ continue; }
        if( inotri == 0 ){ iflg += 1; }
        if( inotri == 1 ){ iflg += 2; }
        if( inotri == 2 ){ iflg += 4; }
      }
      if( iflg == 1 ){
        info.iedA = 2; info.r0A = +aH[0]/(aH[0]-aH[1]);
        info.iedB = 1; info.r0B = -aH[2]/(aH[0]-aH[2]);
      }
      if( iflg == 2 ){
        info.iedA = 0; info.r0A = +aH[1]/(aH[1]-aH[2]);
        info.iedB = 2; info.r0B = -aH[0]/(aH[1]-aH[0]);
      }
      if( iflg == 4 ){
        info.iedA = 1; info.r0A = +aH[2]/(aH[2]-aH[0]);
        info.iedB = 0; info.r0B = -aH[1]/(aH[2]-aH[1]);
      }
      if( iflg == 3 ){
        info.iedA = 0; info.r0A = -aH[1]/(aH[2]-aH[1]);
        info.iedB = 1; info.r0B = +aH[2]/(aH[2]-aH[0]);
      }
      if( iflg == 5 ){
        info.iedA = 2; info.r0A = -aH[0]/(aH[1]-aH[0]);
        info.iedB = 0; info.r0B = +aH[1]/(aH[1]-aH[2]);
      }
      if( iflg == 6 ){
        info.iedA = 1; info.r0A = -aH[2]/(aH[0]-aH[2]);
        info.iedB = 2; info.r0B = +aH[0]/(aH[0]-aH[1]);
      }
      info.pA[0] = (1.0-info.r0A)*aP[(info.iedA+1)%3][0] + (info.r0A)*aP[(info.iedA+2)%3][0];
      info.pA[1] = (1.0-info.r0A)*aP[(info.iedA+1)%3][1] + (info.r0A)*aP[(info.iedA+2)%3][1];
      info.pA[2] = (1.0-info.r0A)*aP[(info.iedA+1)%3][2] + (info.r0A)*aP[(info.iedA+2)%3][2];
      info.pB[0] = (1.0-info.r0B)*aP[(info.iedB+1)%3][0] + (info.r0B)*aP[(info.iedB+2)%3][0];
      info.pB[1] = (1.0-info.r0B)*aP[(info.iedB+1)%3][1] + (info.r0B)*aP[(info.iedB+2)%3][1];
      info.pB[2] = (1.0-info.r0B)*aP[(info.iedB+1)%3][2] + (info.r0B)*aP[(info.iedB+2)%3][2];
    }
    cs.aTriInfo.push_back(info);
    unsigned int iedge_next = info.iedB;
    ////////////////
    int itri_next1 = aTriSur[jtri0*6+iedge_next*2+0];
    if( itri_next1 == -1 ){ break; } // open loop discard
    int iseg_next1 = Tri2Seg[itri_next1];
    assert( iseg_next1 != -1 );
    if( iseg_next1 == iseg_ker ){ return true; } // is_open == false
    iseg_next = iseg_next1;
  }
  // reach here if the loop is open
  cs.aTriInfo.clear();
  /////
  iseg_next = iseg_ker;
  for(;;){ // seg in cs loop
    int jtri0 = aCST[iseg_next];
    unsigned int iflg = 0;
    double aP[3][3], aH[3];
    for(unsigned int inotri=0;inotri<3;inotri++){
      unsigned int jno0 = aTri[jtri0*3+inotri];
      aP[inotri][0] = aXYZ[jno0*3+0];
      aP[inotri][1] = aXYZ[jno0*3+1];
      aP[inotri][2] = aXYZ[jno0*3+2];
      aH[inotri]
      = (aP[inotri][0]-origin[0])*norm[0]
      + (aP[inotri][1]-origin[1])*norm[1]
      + (aP[inotri][2]-origin[2])*norm[2] -height;
      if( aH[inotri] < 0 ){ continue; }
      if( inotri == 0 ){ iflg += 1; }
      if( inotri == 1 ){ iflg += 2; }
      if( inotri == 2 ){ iflg += 4; }
    }
    unsigned int iedge_next;
    if( iflg == 1 ){ iedge_next = 2; } // 0
    if( iflg == 2 ){ iedge_next = 0; } // 1
    if( iflg == 4 ){ iedge_next = 1; } // 2
    if( iflg == 3 ){ iedge_next = 0; } // 01
    if( iflg == 5 ){ iedge_next = 2; } // 02
    if( iflg == 6 ){ iedge_next = 1; } // 12
    int itri_next1 = aTriSur[jtri0*6+iedge_next*2+0];
    if( itri_next1 == -1 ){ break; } // reached boundary
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


void Slice_MeshTri3D_Heights
(std::vector<CSliceTriMesh>& aCS,
 ////
 const std::vector<double>& aHeight,
 const double norm[3],
 const double origin[3],
 const std::vector<double>& aXYZ,
 const std::vector<unsigned int>& aTri,
 const std::vector<int>& aTriSur)
{
  const unsigned int ntri = (unsigned int)aTri.size()/3;
  const unsigned int nH = (unsigned int)aHeight.size();
  ////
  std::vector< std::vector<unsigned int> > aCST;
  IndexElement_OverlapLevels_MeshTri3D(aCST,
                                       aHeight,norm,origin,
                                       aXYZ,aTri);
  /////
  aCS.clear();
  std::vector<int> Tri2Seg;
  Tri2Seg.resize(ntri,-1);
  for(unsigned int ih=0;ih<nH;ih++){ // h loop
    for(unsigned int isg=0;isg<aCST[ih].size();isg++){
      unsigned int itri = aCST[ih][isg];
      Tri2Seg[itri] = isg;
    }
    ////
    std::vector<int> aFlgSeg;
    aFlgSeg.resize(aCST[ih].size(),0);
    unsigned int iseg_ker = 0;
    for(;;){ // cs loop
      for(;iseg_ker<aCST[ih].size();iseg_ker++){
        if( aFlgSeg[iseg_ker] == 0 ){ break; }
      }
      if( iseg_ker == aCST[ih].size() ) break;
      /////
      CSliceTriMesh cs(ih);
      const bool is_closed = TraverseBoundaryLoop(cs, aFlgSeg,
                                                  iseg_ker, ih, Tri2Seg,
                                                  aCST[ih], norm, origin, aHeight[ih],
                                                  aXYZ, aTri, aTriSur);
      if( !is_closed ){ continue; }
      aCS.push_back(cs);
    }
    ////
    for(unsigned int isg=0;isg<aCST[ih].size();isg++){
      unsigned int itri = aCST[ih][isg];
      Tri2Seg[itri] = -1;
    }
  }
}
