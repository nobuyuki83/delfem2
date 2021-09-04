/*
 * Copyright (c) 2019 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */

#include "delfem2/cad3d.h"

#include <deque>
#include <set>

#include "delfem2/cagedef.h"
#include "delfem2/opengl/old/v3q.h"
#include "delfem2/gizmo_geo3.h"
#include "delfem2/mshuni.h"
#include "delfem2/mshmisc.h"
#include "delfem2/geoproximity3_v3.h"
#include "delfem2/geo3_v23m34q.h"
#include "delfem2/geo_polyline2.h"
#include "delfem2/dtri2_v2dtri.h"
#include "delfem2/pgeo.h"

// =========================================

namespace delfem2{
namespace cad3d{

DFM2_INLINE void FaceCenterNormal(
    CVec3d& cg,
    CVec3d& nf,
    const std::vector< std::pair<unsigned int,bool> >& aIE,
    const std::vector<CCad3D_Edge>& aEdge)
{
  const std::size_t nIE = aIE.size();
  cg.setZero();
  double len_tot = 0.0;
  for (unsigned int iie = 0; iie<nIE; ++iie){
    int ie0 = aIE[(iie+0)%nIE].first;
    bool dir0 = aIE[(iie+0)%nIE].second;
    CVec3d pA = dir0 ? aEdge[ie0].p0 : aEdge[ie0].p1;
    CVec3d pB = dir0 ? aEdge[ie0].p1 : aEdge[ie0].p0;
    double lenAB = Distance(pA, pB);
    cg += (pA+pB)*(0.5*lenAB);
    len_tot += lenAB;
  }
  cg /= len_tot;
  // ---
  nf.setZero();
  for (unsigned int iie = 0; iie<nIE; ++iie){
    int ie0 = aIE[(iie+0)%nIE].first;
//    int ie1 = aIE[(iie+1)%nIE].first;
    bool dir0 = aIE[(iie+0)%nIE].second;
    CVec3d pA = dir0 ? aEdge[ie0].p0 : aEdge[ie0].p1;
    CVec3d pB = dir0 ? aEdge[ie0].p1 : aEdge[ie0].p0;
    nf += ((pA-cg)^(pB-cg));
  }
  nf.normalize();
}

DFM2_INLINE void AddSphere_ZSym(
    std::vector<CCad3D_Vertex>& aVertex,
    std::vector<CCad3D_Edge>& aEdge,
    std::vector<CCad3D_Face>& aFace,
    double elen)
{
  int icp0 = (int)aVertex.size();
  aVertex.emplace_back(CVec3d(-1, 0, 0) ); // icp0+0
  aVertex.emplace_back(CVec3d( 0,+1, 0) ); // icp0+1
  aVertex.emplace_back(CVec3d(+1, 0, 0) ); // icp0+2
  aVertex.emplace_back(CVec3d( 0,-1, 0) ); // icp0+3
  aVertex.emplace_back(CVec3d( 0, 0,+1) ); // icp0+4
  ////
  ////
  int ie0 = (int)aEdge.size();
  CCad3D_Edge e0(icp0+0,icp0+1,true,2);
  CCad3D_Edge e1(icp0+1,icp0+2,true,2);
  CCad3D_Edge e2(icp0+2,icp0+3,true,2);
  CCad3D_Edge e3(icp0+3,icp0+0,true,2);
  CCad3D_Edge e4(icp0+4,icp0+0,false,1);
  CCad3D_Edge e5(icp0+4,icp0+1,false,0);
  CCad3D_Edge e6(icp0+2,icp0+4,false,1);
  CCad3D_Edge e7(icp0+3,icp0+4,false,0);
  aEdge.push_back(e0);
  aEdge.push_back(e1);
  aEdge.push_back(e2);
  aEdge.push_back(e3);
  aEdge.push_back(e4);
  aEdge.push_back(e5);
  aEdge.push_back(e6);
  aEdge.push_back(e7);
  /////
  int ifc0 = (int)aFace.size();
  { // face041
    std::vector< std::pair<unsigned int,bool> > aIE;
    aIE.emplace_back(ie0+0,false );
    aIE.emplace_back(ie0+4,false );
    aIE.emplace_back(ie0+5,true );
    aFace.emplace_back(aIE);
  }
  { // face041
    std::vector< std::pair<unsigned int,bool> > aIE;
    aIE.emplace_back(ie0+1,false );
    aIE.emplace_back(ie0+5,false );
    aIE.emplace_back(ie0+6,false );
    aFace.emplace_back(aIE);
  }
  { // face041
    std::vector< std::pair<unsigned int,bool> > aIE;
    aIE.emplace_back(ie0+2,false );
    aIE.emplace_back(ie0+6,true );
    aIE.emplace_back(ie0+7,false );
    aFace.emplace_back(aIE);
  }
  { // face041
    std::vector< std::pair<unsigned int,bool> > aIE;
    aIE.emplace_back(ie0+3,false );
    aIE.emplace_back(ie0+7,true );
    aIE.emplace_back(ie0+4,true );
    aFace.emplace_back(aIE);
  }
  {
    for(auto & iv : aVertex){
      iv.isConst[0] = false;
      iv.isConst[1] = false;
      iv.isConst[2] = false;
    }
    for(auto & ie : aEdge){
      int iv0 = ie.iv0;
      int iv1 = ie.iv1;
      int inorm = ie.inorm;
      if( inorm < 0 || inorm >= 3 ){ continue; }
      aVertex[iv0].isConst[inorm] = true;
      aVertex[iv1].isConst[inorm] = true;
    }
  }
  elen = 0.1;
  for(int ie=ie0;ie<ie0+8;ie++){
    aEdge[ie].Initialize(aVertex,elen); // ie0+0
  }
  for(int ifc=ifc0;ifc<ifc0+4;ifc++){
    aFace[ifc].Initialize(aVertex,aEdge, elen); // ie0+0
  }
  MakeItSmooth(aVertex,aEdge,aFace);
}

DFM2_INLINE void AddTorus_XSym(
    std::vector<CCad3D_Vertex>& aVertex,
    std::vector<CCad3D_Edge>& aEdge,
    std::vector<CCad3D_Face>& aFace,
    double elen)
{
  int icp0 = (int)aVertex.size();
  aVertex.emplace_back(CVec3d(0,-1.0, 0.0) ); // icp0+0
  aVertex.emplace_back(CVec3d(0,-0.3, 0.0) ); // icp0+0
  aVertex.emplace_back(CVec3d(0, 0.0,+1.0) ); // icp0+1
  aVertex.emplace_back(CVec3d(0, 0.0,+0.3) ); // icp0+1
  aVertex.emplace_back(CVec3d(0,+1.0, 0.0) ); // icp0+2
  aVertex.emplace_back(CVec3d(0,+0.3, 0.0) ); // icp0+2
  aVertex.emplace_back(CVec3d(0, 0.0,-1.0) ); // icp0+3
  aVertex.emplace_back(CVec3d(0, 0.0,-0.3) ); // icp0+3
  aVertex[icp0+0].norm = CVec3d(0,-1,0);
  aVertex[icp0+1].norm = CVec3d(0,+1,0);
  aVertex[icp0+2].norm = CVec3d(0,0,+1);
  aVertex[icp0+3].norm = CVec3d(0,0,-1);
  aVertex[icp0+4].norm = CVec3d(0,+1,0);
  aVertex[icp0+5].norm = CVec3d(0,-1,0);
  aVertex[icp0+6].norm = CVec3d(0,0,-1);
  aVertex[icp0+7].norm = CVec3d(0,0,+1);
  /////
  int ie0 = (int)aEdge.size();
  aEdge.emplace_back(icp0+0,icp0+2,true,0 ); // 0
  aEdge.emplace_back(icp0+2,icp0+4,true,0 ); // 1
  aEdge.emplace_back(icp0+4,icp0+6,true,0 ); // 2
  aEdge.emplace_back(icp0+6,icp0+0,true,0 ); // 3
  aEdge.emplace_back(icp0+3,icp0+1,true,0 ); // 4
  aEdge.emplace_back(icp0+5,icp0+3,true,0 ); // 5
  aEdge.emplace_back(icp0+7,icp0+5,true,0 ); // 6
  aEdge.emplace_back(icp0+1,icp0+7,true,0 ); // 7
  aEdge.emplace_back(icp0+1,icp0+0,false,2 ); // 8
  aEdge.emplace_back(icp0+3,icp0+2,false,1 ); // 9
  aEdge.emplace_back(icp0+4,icp0+5,false,2 ); // 10
  aEdge.emplace_back(icp0+6,icp0+7,false,1 ); // 11
  for(int ie=ie0+8;ie<ie0+12;++ie){
    aEdge[ie].r0 = 1;
    aEdge[ie].r1 = 1;
  }
  for(int ie=ie0;ie<ie0+12;ie++){
    aEdge[ie].Initialize(aVertex,elen); // ie0+0
  }
  //
  int ifc0 = (int)aFace.size();
  { // face041
    std::vector< std::pair<unsigned int,bool> > aIE;
    aIE.emplace_back(ie0+0,false );
    aIE.emplace_back(ie0+8,false );
    aIE.emplace_back(ie0+4,false );
    aIE.emplace_back(ie0+9,true );
    aFace.emplace_back(aIE);
  }
  { // face041
    std::vector< std::pair<unsigned int,bool> > aIE;
    aIE.emplace_back(ie0+1,false );
    aIE.emplace_back(ie0+9,false );
    aIE.emplace_back(ie0+5,false );
    aIE.emplace_back(ie0+10,false );
    aFace.emplace_back(aIE);
  }
  { // face041
    std::vector< std::pair<unsigned int,bool> > aIE;
    aIE.emplace_back(ie0+2,false );
    aIE.emplace_back(ie0+10,true );
    aIE.emplace_back(ie0+6,false );
    aIE.emplace_back(ie0+11,false );
    aFace.emplace_back(aIE);
  }
  { // face041
    std::vector< std::pair<unsigned int,bool> > aIE;
    aIE.emplace_back(ie0+3,false );
    aIE.emplace_back(ie0+11,true );
    aIE.emplace_back(ie0+7,false );
    aIE.emplace_back(ie0+8,true );
    aFace.emplace_back(aIE);
  }
  {
    for(auto & iv : aVertex){
      iv.isConst[0] = false;
      iv.isConst[1] = false;
      iv.isConst[2] = false;
    }
    for(auto & ie : aEdge){
      int iv0 = ie.iv0;
      int iv1 = ie.iv1;
      int inorm = ie.inorm;
      if( inorm < 0 || inorm >= 3 ){ continue; }
      aVertex[iv0].isConst[inorm] = true;
      aVertex[iv1].isConst[inorm] = true;
    }
  }
  for(int ifc=ifc0;ifc<ifc0+4;ifc++){
    aFace[ifc].Initialize(aVertex,aEdge, elen); // ie0+0
  }
//  MakeItSmooth(aVertex,aEdge,aFace);
}

DFM2_INLINE void AddCube(
    std::vector<CCad3D_Vertex>& aVertex,
    std::vector<CCad3D_Edge>& aEdge,
    std::vector<CCad3D_Face>& aFace,
    double elen)
{
  int iv0 = (int)aVertex.size();
  aVertex.emplace_back(CVec3d(-1, -1, -1)); // icp0+0
  aVertex.emplace_back(CVec3d(-1, -1, +1)); // icp0+1
  aVertex.emplace_back(CVec3d(-1, +1, -1)); // icp0+2
  aVertex.emplace_back(CVec3d(-1, +1, +1)); // icp0+3
  aVertex.emplace_back(CVec3d(+1, -1, -1)); // icp0+4
  aVertex.emplace_back(CVec3d(+1, -1, +1)); // icp0+5
  aVertex.emplace_back(CVec3d(+1, +1, -1)); // icp0+6
  aVertex.emplace_back(CVec3d(+1, +1, +1)); // icp0+7
  ////
  int ie0 = (int)aEdge.size();
  aEdge.emplace_back(iv0+0, iv0+1, false, 0); // 0
  aEdge.emplace_back(iv0+1, iv0+3, false, 0); // 1
  aEdge.emplace_back(iv0+3, iv0+2, false, 0); // 2
  aEdge.emplace_back(iv0+2, iv0+0, false, 0); // 3
  /////
  aEdge.emplace_back(iv0+6, iv0+4, false, 0); // 4
  aEdge.emplace_back(iv0+7, iv0+6, false, 0); // 5
  aEdge.emplace_back(iv0+5, iv0+7, false, 0); // 6
  aEdge.emplace_back(iv0+4, iv0+5, false, 0); // 7
  /////
  aEdge.emplace_back(iv0+0, iv0+4, false, 1); // 8
  aEdge.emplace_back(iv0+5, iv0+1, false, 1); // 9
  aEdge.emplace_back(iv0+2, iv0+6, false, 1); // 10
  aEdge.emplace_back(iv0+7, iv0+3, false, 1); // 11
  /////
  int ifc0 = (int)aFace.size();
  { // face0132
    std::vector< std::pair<unsigned int, bool> > aIE;
    aIE.emplace_back(ie0+0, true);
    aIE.emplace_back(ie0+1, true);
    aIE.emplace_back(ie0+2, true);
    aIE.emplace_back(ie0+3, true);
    aFace.emplace_back(aIE);
  }
  { // face4567
    std::vector< std::pair<unsigned int, bool> > aIE;
    aIE.emplace_back(ie0+4, false);
    aIE.emplace_back(ie0+5, false);
    aIE.emplace_back(ie0+6, false);
    aIE.emplace_back(ie0+7, false);
    aFace.emplace_back(aIE);
  }
  { // face0451
    std::vector< std::pair<unsigned int, bool> > aIE;
    aIE.emplace_back(ie0+0, false);
    aIE.emplace_back(ie0+8, true);
    aIE.emplace_back(ie0+7, true);
    aIE.emplace_back(ie0+9, true);
    aFace.emplace_back(aIE);
  }
  { // face041
    std::vector< std::pair<unsigned int, bool> > aIE;
    aIE.emplace_back(ie0+2,  false);
    aIE.emplace_back(ie0+11, false);
    aIE.emplace_back(ie0+5,  true);
    aIE.emplace_back(ie0+10, false);
    aFace.emplace_back(aIE);
  }
  { // face041
    std::vector< std::pair<unsigned int, bool> > aIE;
    aIE.emplace_back(ie0+3,  false);
    aIE.emplace_back(ie0+10, true);
    aIE.emplace_back(ie0+4,  true);
    aIE.emplace_back(ie0+8,  false);
    aFace.emplace_back(aIE);
  }
  { // face041
    std::vector< std::pair<unsigned int, bool> > aIE;
    aIE.emplace_back(ie0+1,  false);
    aIE.emplace_back(ie0+9,  false);
    aIE.emplace_back(ie0+6,  true);
    aIE.emplace_back(ie0+11, true);
    aFace.emplace_back(aIE);
  }
  {
    for (int iv = iv0; iv<iv0+8; ++iv){
      aVertex[iv].isConst[0] = false;
      aVertex[iv].isConst[1] = false;
      aVertex[iv].isConst[2] = false;
    }
    for (int ie = ie0; ie<ie0+12; ++ie){
      int iv0e = aEdge[ie].iv0;
      int iv1e = aEdge[ie].iv1;
      int inorm = aEdge[ie].inorm;
      if (inorm < 0||inorm>=3){ continue; }
      aVertex[iv0e].isConst[inorm] = true;
      aVertex[iv1e].isConst[inorm] = true;
    }
  }

  elen = 0.1;
  for (int ie = ie0; ie<ie0+12; ie++){
    aEdge[ie].Initialize(aVertex, elen); // ie0+0
  }
  for (int ifc = ifc0; ifc<ifc0+6; ifc++){
    aFace[ifc].Initialize(aVertex, aEdge, elen); // ie0+0
  }
  MakeItSmooth(aVertex, aEdge, aFace);
}


} // cad3d
} // delfem2



// ------------------------------------


bool delfem2::CCad3D_Edge::isPick(
	double& ratio,
	const CVec2d& sp0,
	const float mMV[16], 
	const float mPj[16]) const
{
  const int np = (int)aP.size();
  for(int ie=0;ie<np-1;ie++){
    int ip0 = ie;
    int ip1 = ie+1;
    CVec3d p0a = aP[ip0];
    CVec3d p1a = aP[ip1];
    CVec2d s0 = screenXYProjection(p0a, mMV, mPj);
    CVec2d s1 = screenXYProjection(p1a, mMV, mPj);
    double dist = GetDist_LineSeg_Point(sp0,s0,s1);
    if( dist < 0.03 ){
      ratio = (ip0+0.5)/(np-1.0);
      return true;
    }
  }
  return false;
}

// -----------------------------------------------

void delfem2::CCad3D_Face::Initialize(
    const std::vector<CCad3D_Vertex>& aVertex,
    const std::vector<CCad3D_Edge>& aEdge,
    [[maybe_unused]] double elen_)
{
  aPInfo.resize(0);
  std::vector<double> aXYZ_B0;
  std::vector<double> aXYZ_B1;
  const unsigned int ne = static_cast<unsigned int>(aIE.size());
  for(std::size_t iie=0;iie<aIE.size();++iie){
    unsigned int ie0 = aIE[iie].first;
    assert( ie0<aEdge.size() );
    const CCad3D_Edge& e0 = aEdge[ie0];
    const bool dir0 = aIE[iie].second;
    int iv0 = (dir0) ? e0.iv0 : e0.iv1;
    {
      CVec3d p0 = (dir0) ? e0.p0 : e0.p1;
      aXYZ_B1.push_back(p0.x);
      aXYZ_B1.push_back(p0.y);
      aXYZ_B1.push_back(p0.z);
    }
    const unsigned nep = static_cast<unsigned int>(e0.aP.size());
    assert( nep > 0 );
    for(unsigned int iep=0;iep<nep-1;++iep){
      unsigned int iep0 = (dir0) ? iep : nep-1-iep;
      double ratio = (double)iep0/(nep-1.0);
      CVec3d pep = (1-ratio)*e0.p0 + ratio*e0.p1;
      aXYZ_B0.push_back(pep.x);
      aXYZ_B0.push_back(pep.y);
      aXYZ_B0.push_back(pep.z);
      CFacePointInfo pinfo;
      if( iep==0 ){
        pinfo.itype = 0;
        pinfo.iv = iv0;
      }
      else{
        pinfo.itype = 1;
        pinfo.iv = -1;
      }
      pinfo.ie = ie0;
      pinfo.iep = iep0;
      aPInfo.push_back(pinfo);
    }
#ifndef NDEBUG
	{
      unsigned int iie1 = (iie+ne-1)%ne; // back
      unsigned int ie1 = aIE[iie1].first;
      assert( ie1<aEdge.size() );
      const CCad3D_Edge& e1 = aEdge[ie1];
      bool dir1 = aIE[iie1].second;
      int iv1 = (dir1) ? e1.iv1 : e1.iv0;
      assert( iv0 == iv1 );
    }
#endif
  }
  ////
  CVec3d cg, norm; cad3d::FaceCenterNormal(cg,norm, aIE, aEdge);
  CVec3d axis_x, axis_y; GetVertical2Vector(norm, axis_x, axis_y);
  std::vector<double> aXY_B0;
  for(unsigned int ixyz=0;ixyz<aXYZ_B0.size()/3;++ixyz){
    CVec3d p(aXYZ_B0[ixyz*3+0],aXYZ_B0[ixyz*3+1],aXYZ_B0[ixyz*3+2]);
    aXY_B0.push_back((p-cg).dot(axis_x));
    aXY_B0.push_back((p-cg).dot(axis_y));
  }
  std::vector<double> aXY_B1;
  for(unsigned int ixyz=0;ixyz<aXYZ_B1.size()/3;++ixyz){
    CVec3d p(aXYZ_B1[ixyz*3+0],aXYZ_B1[ixyz*3+1],aXYZ_B1[ixyz*3+2]);
    aXY_B1.push_back((p-cg).dot(axis_x));
    aXY_B1.push_back((p-cg).dot(axis_y));
  }
  std::vector<double> aXY_out;
  {
    std::vector< std::vector<double> > aaXY;
    aaXY.push_back(aXY_B0);
    /////
    std::vector<int> loopIP_ind,loopIP;
    std::vector<CVec2d> aVec2;
    double elen0 = 0.05;
    {
      JArray_FromVecVec_XY(loopIP_ind,loopIP, aVec2,
                           aaXY);
      if( !CheckInputBoundaryForTriangulation(loopIP_ind,aVec2) ){
        return;
      }
      FixLoopOrientation(loopIP,
                         loopIP_ind,aVec2);
    }
    {
      std::vector<CDynPntSur> aPo2D;
      std::vector<CDynTri> aETri;
      Meshing_SingleConnectedShape2D(aPo2D, aVec2, aETri,
                                     loopIP_ind,loopIP);
      if( elen0 > 1.0e-10 ){
        CInputTriangulation_Uniform param(1.0);
        std::vector<int> aFlgPnt(aVec2.size());
        std::vector<unsigned int> aFlgTri(aETri.size());
        MeshingInside(aPo2D,aETri,aVec2, aFlgPnt,aFlgTri,
                      aVec2.size(), 0, elen0, param);
      }
      MeshTri2D_Export(aXY_out,aTri, aVec2,aETri);
    }
  }
  const std::size_t nxy_bound = aXY_B0.size()/2;
  for(std::size_t ip=nxy_bound;ip<aXY_out.size()/2;++ip){
    double x0 = aXY_out[ip*2+0];
    double y0 = aXY_out[ip*2+1];
    CFacePointInfo pinfo;
    pinfo.itype = 2;
    pinfo.aW0.resize(aXY_B0.size()/2);
    pinfo.aW1.resize(aXY_B1.size()/2);
    MeanValueCoordinate_Polygon2<CVec2d>(
        pinfo.aW0.data(),
        x0,y0,aXY_B0.data(),(unsigned int)(aXY_B0.size()/2));
    MeanValueCoordinate_Polygon2<CVec2d>(
        pinfo.aW1.data(),
        x0,y0,aXY_B1.data(),(unsigned int)(aXY_B1.size()/2));
    aPInfo.push_back(pinfo);
  }
  MovePoints(aVertex,aEdge);
}

void delfem2::CCad3D_Face::MovePoints(
    const std::vector<CCad3D_Vertex>& aVertex,
    const std::vector<CCad3D_Edge>& aEdge)
{
  aXYZ.resize(aPInfo.size()*3);
  for(unsigned int ip=0;ip<aPInfo.size();++ip){
    if( aPInfo[ip].itype == 0 ){
      int iv0 = aPInfo[ip].iv;
      aPInfo[ip].n = aVertex[iv0].norm;
      aXYZ[ip*3+0] = aVertex[iv0].pos.x;
      aXYZ[ip*3+1] = aVertex[iv0].pos.y;
      aXYZ[ip*3+2] = aVertex[iv0].pos.z;
    }
    else if( aPInfo[ip].itype == 1 ){
      int ie0 = aPInfo[ip].ie;
      int iep0 = aPInfo[ip].iep;
      CVec3d ne = aEdge[ie0].getNorm();
      CVec3d te = aEdge[ie0].GetTangentInEdge((double)iep0/(aEdge[ie0].aP.size()-1));
      CVec3d nep = ne^te;
      nep.normalize();
      aPInfo[ip].n = nep;
      aXYZ[ip*3+0] = aEdge[ie0].aP[iep0].x;
      aXYZ[ip*3+1] = aEdge[ie0].aP[iep0].y;
      aXYZ[ip*3+2] = aEdge[ie0].aP[iep0].z;
    }
  }
  //////
  if( aIE.size() == 3 ){
    CVec3d aP[9] = {
      aEdge[aIE[0].first].getVtxPos(aIE[0].second,0),
      aEdge[aIE[0].first].getVtxPos(aIE[0].second,1),
      aEdge[aIE[0].first].getVtxPos(aIE[0].second,2),
      ////
      aEdge[aIE[1].first].getVtxPos(aIE[1].second,0),
      aEdge[aIE[1].first].getVtxPos(aIE[1].second,1),
      aEdge[aIE[1].first].getVtxPos(aIE[1].second,2),
      ////
      aEdge[aIE[2].first].getVtxPos(aIE[2].second,0),
      aEdge[aIE[2].first].getVtxPos(aIE[2].second,1),
      aEdge[aIE[2].first].getVtxPos(aIE[2].second,2),
    };
    for(unsigned int ip=0;ip<aPInfo.size();++ip){
      if( aPInfo[ip].itype != 2 ){ continue; }
      const std::vector<double>& aW1 = aPInfo[ip].aW1;
      assert( aW1.size() == 3 );
      CVec3d p =  delfem2::getPointCoonsTri_CubicBezierEdge(aW1[0],aW1[1],aW1[2],aP);
      aXYZ[ip*3+0] = p.x;
      aXYZ[ip*3+1] = p.y;
      aXYZ[ip*3+2] = p.z;
    }
  }
  else if( aIE.size() == 4 ){
    CVec3d aP[12] = {
      aEdge[aIE[0].first].getVtxPos(aIE[0].second,0),
      aEdge[aIE[0].first].getVtxPos(aIE[0].second,1),
      aEdge[aIE[0].first].getVtxPos(aIE[0].second,2),
      ////
      aEdge[aIE[1].first].getVtxPos(aIE[1].second,0),
      aEdge[aIE[1].first].getVtxPos(aIE[1].second,1),
      aEdge[aIE[1].first].getVtxPos(aIE[1].second,2),
      ////
      aEdge[aIE[2].first].getVtxPos(aIE[2].second,0),
      aEdge[aIE[2].first].getVtxPos(aIE[2].second,1),
      aEdge[aIE[2].first].getVtxPos(aIE[2].second,2),
      ////
      aEdge[aIE[3].first].getVtxPos(aIE[3].second,0),
      aEdge[aIE[3].first].getVtxPos(aIE[3].second,1),
      aEdge[aIE[3].first].getVtxPos(aIE[3].second,2),
    };
    for(unsigned int ip=0;ip<aPInfo.size();++ip){
      if( aPInfo[ip].itype != 2 ){ continue; }
      const std::vector<double>& aW1 = aPInfo[ip].aW1;
      assert( aW1.size() == 4 );
      const double u = aW1[1] + aW1[2];
      const double v = aW1[2] + aW1[3];
//      CVector3 p =  getPointCoons_CubicBezier(u,v,aP);
      CVec3d p = getPointHermetianQuad(u,v,aP);
      aXYZ[ip*3+0] = p.x;
      aXYZ[ip*3+1] = p.y;
      aXYZ[ip*3+2] = p.z;
    }
  }
  else{
    for(std::size_t ip=0;ip<aPInfo.size();++ip){
      if( aPInfo[ip].itype != 2 ){ continue; }
      aXYZ[ip*3+0] = 0;
      aXYZ[ip*3+1] = 0;
      aXYZ[ip*3+2] = 0;
      const std::vector<double>& aW = aPInfo[ip].aW0;
      for(std::size_t jp=0;jp<aW.size();++jp){
        aXYZ[ip*3+0] += aW[jp]*aXYZ[jp*3+0];
        aXYZ[ip*3+1] += aW[jp]*aXYZ[jp*3+1];
        aXYZ[ip*3+2] += aW[jp]*aXYZ[jp*3+2];
      }
      /*
      CVector3 pi(aXYZ[ip*3+0],aXYZ[ip*3+1],aXYZ[ip*3+2]);
      for(int jp=0;jp<aW.size();++jp){
        CVector3 pj(aXYZ[jp*3+0],aXYZ[jp*3+1],aXYZ[jp*3+2]);
        const CVector3 nj = aPInfo[jp].n;
        CVector3 dp = ((pj-pi)*nj)*nj*0.8; // control per edge ?
        aXYZ[ip*3+0] += aW[jp]*dp.x;
        aXYZ[ip*3+1] += aW[jp]*dp.y;
        aXYZ[ip*3+2] += aW[jp]*dp.z;
      }
       */
    }
  }
  aNorm.resize(aXYZ.size());
  delfem2::Normal_MeshTri3D(aNorm.data(),
                            aXYZ.data(), aXYZ.size()/3,
                            aTri.data(), aTri.size()/3);
}

// ---------------------------------------------

unsigned int delfem2::AddPointEdge(
    unsigned int ie_div,
    double ratio_edge,
    std::vector<CCad3D_Vertex>& aVertex,
    std::vector<CCad3D_Edge>& aEdge,
    std::vector<CCad3D_Face>& aFace,
    double elen)
{
  if( ie_div >= aEdge.size() ) return UINT_MAX;
  if (ratio_edge < 0.01 || ratio_edge > 0.99) return UINT_MAX;
  const unsigned int iv_new = static_cast<unsigned int>(aVertex.size());
  {
    CVec3d nv;
    for (const auto & fc : aFace){
      for (unsigned int ie = 0; ie<fc.aIE.size(); ++ie){
        if (fc.aIE[ie].first!=ie_div) continue;         
        CVec3d cg, nf;
        cad3d::FaceCenterNormal(cg, nf, fc.aIE, aEdge);
        nv += nf;        
      }
    }
    nv.normalize();
    // --------------
    CVec3d p = aEdge[ie_div].GetPosInEdge(ratio_edge);
    CCad3D_Vertex v(p);
    {
      int ien = aEdge[ie_div].inorm;
      v.isConst[ien] = true;
    }
    v.norm = nv;
    aVertex.push_back(v);
  }
  const int iv0 = aEdge[ie_div].iv0;
  const int iv1 = aEdge[ie_div].iv1;
  {
    aEdge[ie_div].iv0 = iv0;
    aEdge[ie_div].iv1 = iv_new;
    aEdge[ie_div].Initialize(aVertex,elen);
  }
  const unsigned int ie_new = static_cast<unsigned int>(aEdge.size());
  aEdge.emplace_back(iv_new,iv1,aEdge[ie_div].is_sim,aEdge[ie_div].inorm );
  aEdge[ie_new].Initialize(aVertex,elen);
  
  for(auto & fc : aFace){
    for(unsigned int ie=0;ie<fc.aIE.size();++ie){
      if(fc.aIE[ie].first!=ie_div) continue;
      if(fc.aIE[ie].second){
        fc.aIE.insert(fc.aIE.begin()+ie+1,std::make_pair(ie_new,true));
      }
      else{
        fc.aIE.insert(fc.aIE.begin()+ie,std::make_pair(ie_new,false));
      }
      fc.Initialize(aVertex, aEdge, elen);
      break;
    }
  }
  return iv_new;
}

void delfem2::ConectEdge
(int iv0, int iv1, int iface_div, int inorm_new,
 std::vector<CCad3D_Vertex>& aVertex,
 std::vector<CCad3D_Edge>& aEdge,
 std::vector<CCad3D_Face>& aFace,
 double elen)
{
  if( iface_div < 0 || iface_div >= (int)aFace.size() ) return;
  int iie0 = aFace[iface_div].findIIE_CP(iv0,aEdge);
  int iie1 = aFace[iface_div].findIIE_CP(iv1,aEdge);
  if( iie0 == -1 || iie1 == -1 ) return;
  { // move iv0 and iv1
    CVec3d p0 = aVertex[iv0].pos;
    CVec3d p1 = aVertex[iv1].pos;
    CVec3d mid = (p0+p1)*0.5;
    CVec3d n(0,0,0); n[inorm_new] =1;
    aVertex[iv0].pos = p0-((p0-mid).dot(n))*n;
    aVertex[iv1].pos = p1-((p1-mid).dot(n))*n;
  }
  if( inorm_new >= 0 && inorm_new < 3 ){
    aVertex[iv0].isConst[inorm_new] = true;
    aVertex[iv1].isConst[inorm_new] = true;
  }
  for(auto & iie : aFace[iface_div].aIE){
    int ie0 = iie.first;
    int jv0 = aEdge[ie0].iv0;
    int jv1 = aEdge[ie0].iv1;
    if( (jv0==iv0&&jv1==iv1) || (jv0==iv1&&jv1==iv0) ) return;
  }
  ////////////
  const int ie_new = (int)aEdge.size();
  aEdge.emplace_back(iv0,iv1,false,inorm_new );
  aEdge[ie_new].Initialize(aVertex,elen);
  
  const std::vector< std::pair<unsigned int,bool> > aIE = aFace[iface_div].aIE;
  const int nie = (int)aIE.size();
  { // modify exisiting
    std::vector< std::pair<unsigned int,bool> > aIE0;
    aIE0.emplace_back(ie_new,true );
    for(int iie=iie1;iie%nie!=iie0;++iie){
      aIE0.push_back( aIE[iie%nie] );
    }
    aFace[iface_div].aIE = aIE0;
    aFace[iface_div].Initialize(aVertex, aEdge, elen);
  }
  { // make new
    std::vector< std::pair<unsigned int,bool> > aIE0;
    aIE0.emplace_back(ie_new,false );
    for(int iie=iie0;iie%nie!=iie1;++iie){
      aIE0.push_back( aIE[iie%nie] );
    }
    CCad3D_Face face(aIE0);
    face.Initialize(aVertex, aEdge, elen);
    aFace.push_back(face);
  }
  for(auto & ie : aEdge){
    ie.MovePoints(aVertex);
  }
  for(auto & ifc : aFace){
    ifc.MovePoints(aVertex,aEdge);
  }
}

void delfem2::MakeItSmooth
(std::vector<CCad3D_Vertex>& aVertex,
 std::vector<CCad3D_Edge>& aEdge,
 std::vector<CCad3D_Face>& aFace)
{
  for(auto & iv : aVertex){
    iv.norm.setZero();
  }
  for(auto & ifc : aFace){
    const std::vector< std::pair<unsigned int,bool> >& aIE = ifc.aIE;
    int nIE = (int)aIE.size();
    CVec3d nf,cg; cad3d::FaceCenterNormal(cg,nf,aIE,aEdge);
    nf.normalize();
    for(int iie=0;iie<nIE;++iie){
      int ie0 = aIE[iie].first;
      bool dir0 = aIE[iie].second;
      int ipA = dir0 ? aEdge[ie0].iv0 : aEdge[ie0].iv1;
      aVertex[ipA].norm += nf;
    }
  }
  for(auto & iv : aVertex){
//    if( aVertex[iv].isConst[0] ){ aVertex[iv].norm.x = 0; }
//    if( aVertex[iv].isConst[1] ){ aVertex[iv].norm.y = 0; }
//    if( aVertex[iv].isConst[2] ){ aVertex[iv].norm.z = 0; }
    if( iv.norm.norm() < 0.1 ){
      iv.norm.setZero();
      continue;
    }
    iv.norm.normalize();
  }
  for(auto & ie : aEdge){
    ie.MovePoints(aVertex); // ie0+0
  }
  for(auto & ifc : aFace){
    ifc.MovePoints(aVertex,aEdge); // ie0+0
  }
}

void delfem2::findEdgeGroup(
    std::vector< std::pair<unsigned int,bool> >& aIE,
    unsigned int iedge0,
    [[maybe_unused]] std::vector<CCad3D_Vertex>& aVertex,
    std::vector<CCad3D_Edge>& aEdge)
{
  aIE.clear();
  if( iedge0 >= aEdge.size() ) return;
  std::deque< std::pair<unsigned int,bool> > deqIE;
  deqIE.emplace_back(iedge0,true );
  bool is_loop = false;
  for(;;){
    const unsigned int ie0 = deqIE.back().first;
    int iv0 = deqIE.back().second ? aEdge[ie0].iv1 : aEdge[ie0].iv0;
    int ine0 = aEdge[ie0].inorm;
    if( iv0 == aEdge[iedge0].iv0 ){ is_loop = true; break; }
    const unsigned int ndeqIE = static_cast<unsigned int>(deqIE.size()); // prev
    for(unsigned int ie=0;ie<aEdge.size();++ie){
      if( ie == ie0 ) continue;
      if(      aEdge[ie].iv0 == iv0 && aEdge[ie].inorm == ine0){
        deqIE.emplace_back(ie,true );
        break;
      }
      else if( aEdge[ie].iv1 == iv0 && aEdge[ie].inorm == ine0 ){
        deqIE.emplace_back(ie,false );
        break;
      }
    }
    if( deqIE.size() == ndeqIE ) break; // couldn't find new one
  }
  if( is_loop ){ aIE.assign(deqIE.begin(), deqIE.end()); return; }
  //
  for(;;){
    const unsigned int ie0 = deqIE.front().first;
    int iv0 = deqIE.front().second ? aEdge[ie0].iv0 : aEdge[ie0].iv1;
    int ine0 = aEdge[ie0].inorm;
    assert( iv0 != aEdge[iedge0].iv1 ); // this should not be loop
    const unsigned int ndeqIE = static_cast<unsigned int>(deqIE.size()); // prev
    for(unsigned int ie=0;ie<aEdge.size();++ie){
      if( ie == ie0 ) continue;
      if(      aEdge[ie].iv0 == iv0 && aEdge[ie].inorm == ine0){
        deqIE.push_front( std::make_pair(ie,false ) );
        break;
      }
      else if( aEdge[ie].iv1 == iv0 && aEdge[ie].inorm == ine0 ){
        deqIE.push_front( std::make_pair(ie,true) );
        break;
      }
    }
    if( deqIE.size() == ndeqIE ) break; // couldn't find new one
  }
  aIE.assign(deqIE.begin(), deqIE.end());
}

void delfem2::AddSphere_XSym
(std::vector<CCad3D_Vertex>& aVertex,
 std::vector<CCad3D_Edge>& aEdge,
 std::vector<CCad3D_Face>& aFace,
 double elen)
{
  int icp0 = (int)aVertex.size();
  aVertex.emplace_back(CVec3d(0,-1, 0) ); // icp0+0
  aVertex.emplace_back(CVec3d(0, 0,+1) ); // icp0+1
  aVertex.emplace_back(CVec3d(0,+1, 0) ); // icp0+2
  aVertex.emplace_back(CVec3d(0, 0,-1) ); // icp0+3
  aVertex.emplace_back(CVec3d(1, 0, 0) ); // icp0+4
  //------------------
  int ie0 = (int)aEdge.size();
  CCad3D_Edge e0(icp0+0,icp0+1,true,0);
  CCad3D_Edge e1(icp0+1,icp0+2,true,0);
  CCad3D_Edge e2(icp0+2,icp0+3,true,0);
  CCad3D_Edge e3(icp0+3,icp0+0,true,0);
  CCad3D_Edge e4(icp0+4,icp0+0,false,2);
  CCad3D_Edge e5(icp0+4,icp0+1,false,1);
  CCad3D_Edge e6(icp0+2,icp0+4,false,2);
  CCad3D_Edge e7(icp0+3,icp0+4,false,1);
  aEdge.push_back(e0);
  aEdge.push_back(e1);
  aEdge.push_back(e2);
  aEdge.push_back(e3);
  aEdge.push_back(e4);
  aEdge.push_back(e5);
  aEdge.push_back(e6);
  aEdge.push_back(e7);
  //
  int ifc0 = (int)aFace.size();
  { // face041
    std::vector< std::pair<unsigned int,bool> > aIE;
    aIE.emplace_back(ie0+0,false );
    aIE.emplace_back(ie0+4,false );
    aIE.emplace_back(ie0+5,true );
    aFace.emplace_back(aIE);
  }
  { // face041
    std::vector< std::pair<unsigned int,bool> > aIE;
    aIE.emplace_back(ie0+1,false );
    aIE.emplace_back(ie0+5,false );
    aIE.emplace_back(ie0+6,false );
    aFace.emplace_back(aIE);
  }
  { // face041
    std::vector< std::pair<unsigned int,bool> > aIE;
    aIE.emplace_back(ie0+2,false );
    aIE.emplace_back(ie0+6,true );
    aIE.emplace_back(ie0+7,false );
    aFace.emplace_back(aIE);
  }
  { // face041
    std::vector< std::pair<unsigned int,bool> > aIE;
    aIE.emplace_back(ie0+3,false );
    aIE.emplace_back(ie0+7,true );
    aIE.emplace_back(ie0+4,true );
    aFace.emplace_back(aIE);
  }
  {
    for(auto & iv : aVertex){
      iv.isConst[0] = false;
      iv.isConst[1] = false;
      iv.isConst[2] = false;
    }
    for(auto & ie : aEdge){
      int iv0 = ie.iv0;
      int iv1 = ie.iv1;
      int inorm = ie.inorm;
      if( inorm < 0 || inorm >= 3 ){ continue; }
      aVertex[iv0].isConst[inorm] = true;
      aVertex[iv1].isConst[inorm] = true;
    }
  }
  elen = 0.1;
  for(int ie=ie0;ie<ie0+8;ie++){
    aEdge[ie].Initialize(aVertex,elen); // ie0+0
  }
  for(int ifc=ifc0;ifc<ifc0+4;ifc++){
    aFace[ifc].Initialize(aVertex,aEdge, elen); // ie0+0
  }
  MakeItSmooth(aVertex,aEdge,aFace);
}



bool delfem2::FindFittingPoint(
    CVec2d& p2d_near,
    CVec2d& p2d_norm,
    const CVec2d& p2d_org,
    const std::vector<CVec2d>& aP2D,
    bool isConstX,
    bool isConstY,
    [[maybe_unused]] double half_view_height)
{
  bool isHit = false;
  if( isConstX &&  isConstY ){ return false; }
  else if( isConstX && !isConstY ){
    for(std::size_t iq=0;iq<aP2D.size()-1;++iq){
      CVec2d q0 = aP2D[iq+0];
      CVec2d q1 = aP2D[iq+1];
      if( (q0.x-p2d_org.x)*(q1.x-p2d_org.x) < 0 ){
        p2d_near.p[0] = p2d_org.x;
        p2d_near.p[1] = ((q0+q1)*0.5).y;
        p2d_norm.p[0] = (q1-q0).y;
        p2d_norm.p[1] = (q0-q1).x;
        isHit = true;
        break;
      }
    }
  }
  else if( !isConstX && isConstY ){
    assert( !aP2D.empty() );
    for(std::size_t iq=0;iq<aP2D.size()-1;++iq){
      CVec2d q0 = aP2D[iq+0];
      CVec2d q1 = aP2D[iq+1];
      if( (q0.y-p2d_org.y)*(q1.y-p2d_org.y) < 0 ){
        p2d_near.p[0] = ((q0+q1)*0.5).x;
        p2d_near.p[1] = p2d_org.y;
        p2d_norm.p[0] = (q1-q0).y;
        p2d_norm.p[1] = (q0-q1).x;
        isHit = true;
        break;
      }
    }
  }
  else{
    double min_dist = -1;
    for(const auto & iq : aP2D){
      double len = (iq-p2d_org).norm();
      if( min_dist < 0 || len < min_dist ){
        min_dist = len;
        p2d_near = iq;
        p2d_norm = iq-p2d_org;
        isHit = true;
      }
    }
  }
  if( !isHit ) return false;
  double dist = (p2d_near-p2d_org).norm();
  if( dist > 0.4 ) return false;
  p2d_norm.normalize();
  return true;
}

std::vector<int> delfem2::getPointsInEdges(
    const std::vector< std::pair<unsigned int,bool> >& aIE_picked,
    const std::vector<CCad3D_Edge>& aEdge)
{
  std::vector<int> aIP;
  for(size_t iie=0;iie<aIE_picked.size()+1;++iie){
    int iv0;
    if( iie != aIE_picked.size() ){
      int ie0 = aIE_picked[iie].first;
      iv0 = aIE_picked[iie].second ? aEdge[ie0].iv0 : aEdge[ie0].iv1;
    }
    else{
      int ie0 = aIE_picked[iie-1].first;
      iv0 = aIE_picked[iie-1].second ? aEdge[ie0].iv1 : aEdge[ie0].iv0;
    }
    aIP.push_back(iv0);
  }
  if( aIP.front() == aIP.back() ){ aIP.pop_back(); }
  return aIP;
}

bool delfem2::MovePointsAlongSketch(
    std::vector<CCad3D_Vertex>& aVertex,
    std::vector<CCad3D_Edge>& aEdge,
    std::vector<CCad3D_Face>& aFace,
    // --
    const std::vector<CVec2d>& aStroke,
    const std::vector< std::pair<unsigned int,bool> >& aIE_picked,
    const CVec3d& plane_org, int inorm,
    float mMV[16],
    float mPj[16],
    double view_height)
{
  // resampling
  std::vector<CVec2d> aStroke1 = Polyline_Resample_Polyline(aStroke,0.025);
  //
  CVec3d plane_nrm(0,0,0); plane_nrm[inorm] = 1;
  CVec3d plane_ex(0,0,0); plane_ex[(inorm+1)%3] = 1;
  CVec3d plane_ey(0,0,0); plane_ey[(inorm+2)%3] = 1;
  std::vector<CVec2d> aP2D;
  for(const auto& sp0 : aStroke1){
    CVec3d src = screenUnProjection(CVec3d(sp0.x,sp0.y,0), mMV, mPj);
    CVec3d dir = screenUnProjection(CVec3d(0,0,1), mMV, mPj);
    CVec3d p = intersection_Plane_Line(plane_org, plane_nrm, src,dir);
    aP2D.emplace_back((p-plane_org).dot(plane_ex),(p-plane_org).dot(plane_ey));
  }
  bool is_moved = false;
  std::vector<int> aIP = getPointsInEdges(aIE_picked,aEdge);
  for(int iv0 : aIP){
    CCad3D_Vertex& v = aVertex[iv0];
    CVec2d p2d_org((v.pos-plane_org).dot(plane_ex), (v.pos-plane_org).dot(plane_ey));
    const bool isConstX = v.isConst[(inorm+1)%3];
    const bool isConstY = v.isConst[(inorm+2)%3];
    CVec2d p2d_near, p2d_norm;
    bool res = FindFittingPoint(p2d_near,p2d_norm,
        p2d_org, aP2D, isConstX,isConstY,view_height*0.2);
    if( res ){
      CVec3d p3d_near = plane_org + p2d_near.x*plane_ex + p2d_near.y*plane_ey;
      CVec3d n3d_near = p2d_norm.x*plane_ex + p2d_norm.y*plane_ey;
      v.pos = p3d_near;
      v.norm = n3d_near;
      is_moved = true;
    }
  }
  ////
  if( is_moved ){
    for(auto & ie : aEdge){
      ie.MovePoints(aVertex);
    }
    for(auto & ifc : aFace){
      ifc.MovePoints(aVertex,aEdge);
    }
  }
  return is_moved;
}


void delfem2::DivideFace
(int ifc,
 const CVec3d& org, int inorm,
 std::vector<CCad3D_Vertex>& aVertex,
 std::vector<CCad3D_Edge>& aEdge,
 std::vector<CCad3D_Face>& aFace,
 double elen)
{
  if( inorm == -1 ) return;
//  const CVector3& plane_ex = CVector3::Axis((inorm+1)%3);
//  const CVector3& plane_ey = CVector3::Axis((inorm+2)%3);
  const CVec3d n01 = CVec3d::Axis(inorm);
  const std::vector< std::pair<unsigned int,bool> > aIE = aFace[ifc].aIE;
  const int nie = (int)aIE.size();
  std::set<int> setIV_new;
  for(int iie=0;iie<nie;++iie){
    const int ie0 = aIE[iie].first;
    double ratio;
    if (!aEdge[ie0].GetParameterIntersection(ratio, org, n01)) continue;
    int iv0 = -1;
    if(      fabs(ratio-1.0)<1.0e-5 ){ iv0 = aEdge[ie0].iv1; }
    else if( fabs(ratio-0.0)<1.0e-5 ){ iv0 = aEdge[ie0].iv0; }
    else if( ratio < 0.01 || ratio > 0.99 ) continue;
    else{ iv0 = AddPointEdge(ie0, ratio, aVertex, aEdge, aFace, elen); }
    setIV_new.insert(iv0);
  }
  if( setIV_new.size() != 2 ) return;
  int iv0 = *(setIV_new.begin());
  int iv1 = *(++setIV_new.begin());
  {
    CVec3d p0 = aVertex[iv0].pos;
    CVec3d p1 = aVertex[iv1].pos;
    CVec3d n0 = aVertex[iv0].norm;
    CVec3d n1 = aVertex[iv1].norm;
    CVec3d v = Cross(n0+n1,n01);
    if( v.dot(p1-p0) > 0 ){ ConectEdge(iv0,iv1, ifc, inorm, aVertex, aEdge, aFace, elen); }
    else{                   ConectEdge(iv1,iv0, ifc, inorm, aVertex, aEdge, aFace, elen); }
  }
}

void delfem2::BuildTriMesh
(std::vector<double>& aXYZ,
 std::vector<unsigned int>& aTri,
 std::vector<unsigned int>& aTriSuTri,
 std::vector<double>& aNorm,
 std::vector<CCad3D_Vertex>& aVertex,
 std::vector<CCad3D_Edge>& aEdge,
 std::vector<CCad3D_Face>& aFace,
 int isym)
{
  std::vector<int> aIsSymVtx(aVertex.size(),0);
  for(const auto & e : aEdge){
    if( e.is_sim ){
      assert( e.inorm == isym );
      int iv0 = e.iv0;
      int iv1 = e.iv1;
      aIsSymVtx[iv0] = 1;
      aIsSymVtx[iv1] = 1;
    }
  }
  aXYZ.resize(0);
  for(std::size_t iv=0;iv<aVertex.size();++iv){
    CCad3D_Vertex& v = aVertex[iv];
    int iq0 = (int)aXYZ.size()/3;
    aXYZ.push_back(+v.pos.x);
    aXYZ.push_back(+v.pos.y);
    aXYZ.push_back(+v.pos.z);
    v.iq_right = iq0;
    if( aIsSymVtx[iv] == 1 ){
      aVertex[iv].iq_left = iq0;
    }
    else{
      if( isym == 2 ){
        aXYZ.push_back(+v.pos.x);
        aXYZ.push_back(+v.pos.y);
        aXYZ.push_back(-v.pos.z);
      }
      else if( isym == 0 ){
        aXYZ.push_back(-v.pos.x);
        aXYZ.push_back(+v.pos.y);
        aXYZ.push_back(+v.pos.z);
      }
      aVertex[iv].iq_left = iq0+1;
    }
  }
  for(auto & e : aEdge){
    const int iv0 = e.iv0;
    const int iv1 = e.iv1;
    assert( iv0>=0 && iv0<(int)aVertex.size() );
    assert( iv1>=0 && iv1<(int)aVertex.size() );
    const std::size_t np = e.aP.size();
    e.aIQ_RightLeft.resize(np*2);
    e.aIQ_RightLeft[0*2+0] = aVertex[iv0].iq_right;
    e.aIQ_RightLeft[0*2+1] = aVertex[iv0].iq_left;
    e.aIQ_RightLeft[(np-1)*2+0] = aVertex[iv1].iq_right;
    e.aIQ_RightLeft[(np-1)*2+1] = aVertex[iv1].iq_left;
    assert(np>=2);
    for(unsigned int ip=1;ip<np-1;++ip){
      int iq0 = (int)aXYZ.size()/3;
      aXYZ.push_back(+e.aP[ip].x);
      aXYZ.push_back(+e.aP[ip].y);
      aXYZ.push_back(+e.aP[ip].z);
      e.aIQ_RightLeft[ip*2+0] = iq0;
      if( e.is_sim ){
        e.aIQ_RightLeft[ip*2+1] = iq0;
      }
      else{
        if( isym == 2 ){
          aXYZ.push_back(+e.aP[ip].x);
          aXYZ.push_back(+e.aP[ip].y);
          aXYZ.push_back(-e.aP[ip].z);
        }
        else if( isym == 0 ){
          aXYZ.push_back(-e.aP[ip].x);
          aXYZ.push_back(+e.aP[ip].y);
          aXYZ.push_back(+e.aP[ip].z);
        }
        e.aIQ_RightLeft[ip*2+1] = iq0+1;
      }
    }
  }
  for(auto & fc : aFace){
    std::size_t np = fc.aPInfo.size();
    for(std::size_t ip=0;ip<np;++ip){
      CCad3D_Face::CFacePointInfo& pinfo = fc.aPInfo[ip];
      if( pinfo.itype == 0 ){
        int iv0 = pinfo.iv;
        pinfo.iq_right = aVertex[iv0].iq_right;
        pinfo.iq_left  = aVertex[iv0].iq_left;
      }
      else if( pinfo.itype == 1 ){
        int ie0 = pinfo.ie;
        int iep0 = pinfo.iep;
        pinfo.iq_right = aEdge[ie0].aIQ_RightLeft[iep0*2+0];
        pinfo.iq_left  = aEdge[ie0].aIQ_RightLeft[iep0*2+1];
      }
      else if( pinfo.itype == 2 ){
        int iq0 = (int)aXYZ.size()/3;
        aXYZ.push_back(+fc.aXYZ[ip*3+0]);
        aXYZ.push_back(+fc.aXYZ[ip*3+1]);
        aXYZ.push_back(+fc.aXYZ[ip*3+2]);
        ////
        if( isym == 2 ){
          aXYZ.push_back(+fc.aXYZ[ip*3+0]);
          aXYZ.push_back(+fc.aXYZ[ip*3+1]);
          aXYZ.push_back(-fc.aXYZ[ip*3+2]);
        }
        else if( isym == 0 ){
          aXYZ.push_back(-fc.aXYZ[ip*3+0]);
          aXYZ.push_back(+fc.aXYZ[ip*3+1]);
          aXYZ.push_back(+fc.aXYZ[ip*3+2]);
        }
        pinfo.iq_right = iq0;
        pinfo.iq_left  = iq0+1;
      }
    }
  }
  aTri.resize(0);
  for(auto & fc : aFace){
    for(unsigned int it=0;it<fc.aTri.size()/3;++it){
      int ip0 = fc.aTri[it*3+0];
      int ip1 = fc.aTri[it*3+1];
      int ip2 = fc.aTri[it*3+2];
      aTri.push_back(fc.aPInfo[ip0].iq_right);
      aTri.push_back(fc.aPInfo[ip1].iq_right);
      aTri.push_back(fc.aPInfo[ip2].iq_right);
      ///
      aTri.push_back(fc.aPInfo[ip0].iq_left);
      aTri.push_back(fc.aPInfo[ip2].iq_left);
      aTri.push_back(fc.aPInfo[ip1].iq_left);
    }
  }
  for(const auto & e : aEdge){
    if( !e.is_sim ){ continue; }
    for(std::size_t ip=0;ip<e.aP.size();++ip){
      int iq0 = e.aIQ_RightLeft[ip*2+0];
      aXYZ[iq0*3+2] += (double)rand()/(RAND_MAX+1.0)*1.0e-5;
    }
  }
  ElSuEl_MeshElem(aTriSuTri,
                  aTri.data(),aTri.size()/3,
                  delfem2::MESHELEM_TRI,
                  (int)aXYZ.size()/3);
  aNorm.resize(aXYZ.size());
  delfem2::Normal_MeshTri3D(aNorm.data(),
                            aXYZ.data(), aXYZ.size()/3,
                            aTri.data(), aTri.size()/3);

}

void delfem2::UpdateTriMesh(
	std::vector<double>& aXYZ,
	std::vector<double>& aNorm,
	// ---------------------
	const std::vector<unsigned int>& aTri,
	[[maybe_unused]] const std::vector<CCad3D_Vertex>& aVertex,
	const std::vector<CCad3D_Edge>& aEdge,
	const std::vector<CCad3D_Face>& aFace,
	int isym)
{
  for(const auto & fc : aFace){
    for(std::size_t ip=0;ip<fc.aPInfo.size();++ip){
      int iq0 = fc.aPInfo[ip].iq_right;
      int iq1 = fc.aPInfo[ip].iq_left;
      assert( iq0 < (int)aXYZ.size()/3 );
      assert( iq1 < (int)aXYZ.size()/3 );
      aXYZ[iq0*3+0] = +fc.aXYZ[ip*3+0];
      aXYZ[iq0*3+1] = +fc.aXYZ[ip*3+1];
      aXYZ[iq0*3+2] = +fc.aXYZ[ip*3+2];
      if( isym == 0 ){
        aXYZ[iq1*3+0] = -fc.aXYZ[ip*3+0];
        aXYZ[iq1*3+1] = +fc.aXYZ[ip*3+1];
        aXYZ[iq1*3+2] = +fc.aXYZ[ip*3+2];
      }
      else{
        aXYZ[iq1*3+0] = +fc.aXYZ[ip*3+0];
        aXYZ[iq1*3+1] = +fc.aXYZ[ip*3+1];
        aXYZ[iq1*3+2] = -fc.aXYZ[ip*3+2];
      }
    }
  }
  for(const auto & e : aEdge){
    if( !e.is_sim ){ continue; }
    for(std::size_t ip=0;ip<e.aP.size();++ip){
      int iq0 = e.aIQ_RightLeft[ip*2+0];
      aXYZ[iq0*3+2] += (double)rand()/(RAND_MAX+1.0)*1.0e-5;
    }
  }
  aNorm.resize(aXYZ.size());
  delfem2::Normal_MeshTri3D(aNorm.data(),
                            aXYZ.data(), aXYZ.size()/3,
                            aTri.data(), aTri.size()/3);

}


void delfem2::CCad3D::Pick(
    const CVec3d& src_pick,
    const CVec3d& dir_pick,
    const CVec2d& sp0,
    float mMV[16],
    float mPj[16],
    double view_height)
{
  if ( ivtx_picked<aVertex.size() ){
    ielem_vtx_picked = 0;
    for (int iaxis = 0; iaxis<3; iaxis++){
      if( aVertex[ivtx_picked].isConst[iaxis] ) continue;
      CVec3d axis(0, 0, 0);
      axis[iaxis] = 1;
      if (isPick_AxisHandler(sp0, aVertex[ivtx_picked].pos,
                             axis, 0.5,
                             mMV, mPj, 0.03)){
        ielem_vtx_picked = iaxis+1;
        break;
      }
    }
  }
  if( ielem_vtx_picked == 0 || ivtx_picked >= aVertex.size() ){
    ivtx_picked = UINT_MAX;
    for(unsigned int icp=0;icp<aVertex.size();++icp){
      const CVec3d& pos = aVertex[icp].pos;
      CVec3d pn = nearest_Line_Point(pos, src_pick, dir_pick);
      if( (pn-pos).norm() < 0.05 ){
        ivtx_picked = icp;
        ielem_vtx_picked = 0;
        break;
      }
    }
  }
  if( ivtx_picked < aVertex.size() ){
    iedge_picked = -1;
    plane_inorm = -1;
    iface_picked = UINT_MAX;
    aIE_picked.clear();
    return;
  }
  
  // edge pick
  if( iedge_picked != -1 ){ // edge was picked
    ielem_edge_picked = 0;
    {
      CVec2d sp = screenXYProjection(aEdge[iedge_picked].q0, mMV, mPj);
      if( (sp0-sp).norm() < 0.05 ){
        ielem_edge_picked = 2;
        iface_picked = UINT_MAX;
        return;
      }
    }
    {
      CVec2d sp = screenXYProjection(aEdge[iedge_picked].q1, mMV, mPj);
      if( (sp0-sp).norm() < 0.05 ){
        ielem_edge_picked = 3;
        iface_picked = UINT_MAX;
        return;
      }
    }
    if( plane_inorm>=0 && plane_inorm<3 ){ // plane pick
      CVec3d plane_ex = CVec3d::Axis((plane_inorm+1)%3);
      CVec3d plane_ey = CVec3d::Axis((plane_inorm+2)%3);
      CVec3d aP[4] = {
        plane_org-plane_sizeX*plane_ex-plane_sizeY*plane_ey,
        plane_org+plane_sizeX*plane_ex-plane_sizeY*plane_ey,
        plane_org+plane_sizeX*plane_ex+plane_sizeY*plane_ey,
        plane_org-plane_sizeX*plane_ex+plane_sizeY*plane_ey };
      CVec2d sp[4]  = {
        screenXYProjection(aP[0], mMV, mPj),
        screenXYProjection(aP[1], mMV, mPj),
        screenXYProjection(aP[2], mMV, mPj),
        screenXYProjection(aP[3], mMV, mPj) };
      double d01 = GetDist_LineSeg_Point(sp0, sp[0],sp[1]);
      double d12 = GetDist_LineSeg_Point(sp0, sp[1],sp[2]);
      double d23 = GetDist_LineSeg_Point(sp0, sp[2],sp[3]);
      double d30 = GetDist_LineSeg_Point(sp0, sp[3],sp[0]);
      if( d01 < 0.05 || d12 < 0.05 || d23 < 0.05 || d30 < 0.05 ) {
        ielem_edge_picked = 1;
        iface_picked = UINT_MAX;
        return;
      }
    }
  }
  //
  plane_inorm = -1;
  iedge_picked = -1;
  aIE_picked.clear();
  {
    std::map<double, std::pair<int, double> > mapDepthEdge;
    for(unsigned int ie=0;ie<aEdge.size();++ie){
      double ratio_edge;
      bool res = aEdge[ie].isPick(ratio_edge, sp0, mMV, mPj);
      if( res ){
        CVec3d p = aEdge[ie].GetPosInEdge(ratio_edge);
        double depth = -p.dot(dir_pick);
        mapDepthEdge.insert( std::make_pair(depth, std::make_pair(ie,ratio_edge) ) );
      }
    }
    if( !mapDepthEdge.empty() ){
      iedge_picked = mapDepthEdge.begin()->second.first;
      ratio_edge_picked = mapDepthEdge.begin()->second.second;
      ielem_vtx_picked = 0;
      iface_picked = UINT_MAX;
    }
  }
  if( iedge_picked != -1 ){ // make plane
    findEdgeGroup(aIE_picked, iedge_picked, aVertex, aEdge);
    plane_inorm = aEdge[iedge_picked].inorm;
    const CVec3d n = CVec3d::Axis(plane_inorm);
    const CVec3d plane_ex = CVec3d::Axis((plane_inorm+1)%3);
    const CVec3d plane_ey = CVec3d::Axis((plane_inorm+2)%3);
    int iv0 = aEdge[iedge_picked].iv0;
    int iv1 = aEdge[iedge_picked].iv1;
    plane_org = (aVertex[iv0].pos+aVertex[iv1].pos)*0.5;
    double minX=0, maxX=0, minY=0, maxY=0;
    for(std::size_t iie=0;iie<aIE_picked.size();++iie){
      int ie = aIE_picked[iie].first;
      std::vector<CVec3d>& aP = aEdge[ie].aP;
      for(std::size_t ip=0;ip<aP.size();++ip){
        const CVec3d& p = aP[ip];
        double x0 = (p-plane_org).dot(plane_ex);
        double y0 = (p-plane_org).dot(plane_ey);
        if( iie==0 && ip==0 ){
          minX = x0;
          minY = y0;
          maxX = x0;
          maxY = y0;
          continue;
        }
        minX = (x0<minX) ? x0 : minX;
        maxX = (x0>maxX) ? x0 : maxX;
        minY = (y0<minY) ? y0 : minY;
        maxY = (y0>maxY) ? y0 : maxY;
      }
    }
    plane_org += ((minX+maxX)*plane_ex + (minY+maxY)*plane_ey)*0.5;
    plane_sizeX = 0.5*(maxX-minX) + view_height*0.3;
    plane_sizeY = 0.5*(maxY-minY) + view_height*0.3;
    return;
  }
  
  iface_picked = UINT_MAX;
  for(unsigned int ifc=0;ifc<aFace.size();++ifc){
    if( aFace[ifc].isPick(src_pick, dir_pick) ){
      iface_picked = ifc;
      return;
    }
  }
}

void delfem2::CCad3D::MouseUp
 (float mMV[16],
  float mPj[16], double view_height)
{
  
  if( imode_edit == EDIT_SKETCH ){
    if( aStroke.size() > 3 && iedge_picked != -1 ){
      bool res = MovePointsAlongSketch(aVertex,aEdge,aFace,
          aStroke,aIE_picked,
          plane_org,aEdge[iedge_picked].inorm,
          mMV,mPj,view_height);
      if( res ){
        UpdateTriMesh(aXYZ,aNorm, aTri,aVertex,aEdge,aFace,isym);
      }
      else{
        plane_inorm = -1;
        iedge_picked = -1;
        aIE_picked.clear();
      }
    }
    aStroke.clear();
  }
}

bool delfem2::CCad3D::ReflectChangeForCurveAndSurface
(std::vector<int>& aIsMoved_Edge,
 const std::vector<int>& aIsMoved_Vtx)
{
  bool is_edit = false;
  for(std::size_t ie=0;ie<aEdge.size();++ie){
    int iv0 = aEdge[ie].iv0;
    int iv1 = aEdge[ie].iv1;
    if( aIsMoved_Vtx[iv0]==0 && aIsMoved_Vtx[iv1]==0 && aIsMoved_Edge[ie]==0 ) continue;
    aIsMoved_Edge[ie] = 1;
    aEdge[ie].MovePoints(aVertex);
    is_edit = true;
  }
  for(auto & ifc : aFace){
    bool is_edit_face = false;
    for(auto & iie : ifc.aIE){
      int ie0 = iie.first;
      if( aIsMoved_Edge[ie0] == 0 ) continue;
      is_edit_face = true;
      break;
    }
    if( !is_edit_face ) continue;
    ifc.MovePoints(aVertex,aEdge);
    is_edit = true;
  }
  if( is_edit ){
    UpdateTriMesh(aXYZ,aNorm,
        aTri,aVertex,aEdge,aFace,isym);
  }
  return is_edit;
}

bool delfem2::CCad3D::MouseMotion
(const CVec3d& src_pick, const CVec3d& dir_pick,
 const CVec2d& sp0, const CVec2d& sp1,
 float mMV[16], float mPj[16])
{
  if( imode_edit == EDIT_MOVE ){
    std::vector<int> aIsMoved_Vtx(aVertex.size(),0);
    std::vector<int> aIsMoved_Edge(aEdge.size(),0);
    if( ivtx_picked < aVertex.size() ){ // move vtx
      if( ielem_vtx_picked <= 0 || ielem_vtx_picked > 3 ){ return false; }
      int iaxis = ielem_vtx_picked-1;
      CVec3d axis(0, 0, 0); axis[iaxis] = 1;
      CVec3d d0 = drag_AxisHandler(sp0, sp1, aVertex[ivtx_picked].pos, axis, 0.5, mMV, mPj);
      aVertex[ivtx_picked].pos += d0;
      aIsMoved_Vtx[ivtx_picked] = 1;
    }
    if( iedge_picked>=0 && iedge_picked<(int)aEdge.size() && (ielem_edge_picked==2 || ielem_edge_picked==3) ){ // moved edge ctrl point
      const int ie0 = iedge_picked;
      CVec3d axis = CVec3d::Axis(aEdge[ie0].inorm);
      CVec3d qe = intersection_Plane_Line(plane_org, axis, src_pick, dir_pick);
      const int iv = (ielem_edge_picked==2) ? aEdge[ie0].iv0 : aEdge[ie0].iv1;
      CVec3d  qs = (ielem_edge_picked==2) ? aEdge[ie0].q0  : aEdge[ie0].q1;
      CVec3d  pv = (ielem_edge_picked==2) ? aEdge[ie0].p0  : aEdge[ie0].p1;
      const double len0 = (aEdge[ie0].p0-aEdge[ie0].p1).norm();
      if( (qs-pv).norm() < 0.01*len0 ) return false;
      if( (qe-pv).norm() < 0.01*len0 ) return false;
      CMat3d R = Mat3_MinimumRotation(qs-pv, qe-pv);
      aVertex[iv].norm = R*(aVertex[iv].norm);
      double len1 = (qe-pv).norm();
      if( ielem_edge_picked==2 ){ aEdge[ie0].r0 = len1/len0; }
      else{                       aEdge[ie0].r1 = len1/len0; }
      aIsMoved_Vtx[iv] = 1;
      aIsMoved_Edge[iedge_picked] = 1;
    }
    if( iedge_picked>=0 && iedge_picked<(int)aEdge.size() && ielem_edge_picked == 1 ){ // move vtx on the plane
      if( !aEdge[iedge_picked].is_sim  ){
        int iaxis = aEdge[iedge_picked].inorm;
        CVec3d axis(0, 0, 0); axis[iaxis] = 1;
        CVec2d spa0 = screenXYProjection(plane_org+axis, mMV, mPj);
        CVec2d spa1 = screenXYProjection(plane_org-axis, mMV, mPj);
        double r = (spa0-spa1).dot(sp1-sp0)/(spa0-spa1).squaredNorm();
        CVec3d d = r*axis;
        plane_org += d;
        std::vector<int> aIP = getPointsInEdges(aIE_picked, aEdge);
        for(int ip0 : aIP){
          aVertex[ip0].pos += d;
          aIsMoved_Vtx[ip0] = 1;
        }
      }
    }
    return ReflectChangeForCurveAndSurface(aIsMoved_Edge,aIsMoved_Vtx);
  }
  else if( imode_edit == EDIT_SKETCH ){
    aStroke.push_back(sp0);
    return false;
  }
  else if( imode_edit == EDIT_ADD_CROSS_SECTION ){
    if( plane_inorm >= 0 && plane_inorm < 3 ){
      int iaxis = aEdge[iedge_picked].inorm;
      if( iaxis>=0 && iaxis<3 ){
        CVec3d axis(0, 0, 0); axis[iaxis] = 1;
        CVec2d spa0 = screenXYProjection(plane_org+axis, mMV, mPj);
        CVec2d spa1 = screenXYProjection(plane_org-axis, mMV, mPj);
        double r = (spa0-spa1).dot(sp1-sp0)/(spa0-spa1).squaredNorm();
        CVec3d d = r*axis;
        plane_org += d;
        return false;
      }
    }
  }
  return false;
}

void delfem2::CCad3D::MouseDown
(const CVec3d& src_pick, const CVec3d& dir_pick,
 const CVec2d& sp0, float mMV[16], float mPj[16],
 double view_height)
{
  if( imode_edit == EDIT_ADD_POINT_EDGE ){
    Pick(src_pick, dir_pick, sp0, mMV, mPj, view_height);
    AddPointEdge(iedge_picked, ratio_edge_picked, aVertex, aEdge, aFace, elen);
    BuildTriMesh(aXYZ,aTri,aTriSuTri,aNorm, aVertex,aEdge,aFace, isym);
    this->iedge_picked = -1;
    this->aIE_picked.clear();
    this->plane_inorm = -1;
    imode_edit = EDIT_NONE;
  }
  else if( imode_edit == EDIT_SKETCH ){
    if( iedge_picked == -1 ){
      Pick(src_pick, dir_pick, sp0, mMV, mPj, view_height);
    }
  }
  else if( imode_edit == EDIT_ADD_CROSS_SECTION ){
    if( plane_inorm != -1 ){
      for(unsigned int ifc=0;ifc<aFace.size();++ifc){
        if( aFace[ifc].isPick(src_pick, dir_pick) ){
          DivideFace(ifc,plane_org,plane_inorm,
              aVertex,aEdge,aFace, elen);
          BuildTriMesh(aXYZ,aTri,aTriSuTri,aNorm, aVertex,aEdge,aFace, isym);
          return;
        }
      }
    }
    Pick(src_pick, dir_pick, sp0, mMV, mPj, view_height);
  }
  else if( imode_edit == EDIT_MOVE ){
    Pick(src_pick, dir_pick, sp0, mMV, mPj, view_height);
  }
}

bool delfem2::CCad3D::isSym(int iv) const{
  for(const auto & ie : aEdge){
    if( !ie.is_sim ) continue;
    if( ie.iv0 == iv ) return true;
    if( ie.iv1 == iv ) return true;
  }
  return false;
}

void delfem2::CCad3D::WriteFile(std::ofstream& fout) const
{
  fout << aVertex.size() << std::endl;
  for(std::size_t iv=0;iv<aVertex.size();++iv){
    fout << " " << iv << std::endl;
    aVertex[iv].WriteFile(fout);
  }
  ////
  fout << aEdge.size() << std::endl;
  for(std::size_t ie=0;ie<aEdge.size();++ie){
    fout << " " << ie << std::endl;
    aEdge[ie].WriteFile(fout);
  }
  ///
  fout << aFace.size() << std::endl;
  for(std::size_t ifc=0;ifc<aFace.size();++ifc){
    fout << " " << ifc << std::endl;
    aFace[ifc].WriteFile(fout);
  }
}

void delfem2::CCad3D::ReadFile(std::ifstream& fin)
{
  int nv = 0;
  fin >> nv;
  aVertex.resize(nv);
  for(std::size_t iv=0;iv<aVertex.size();++iv){
    unsigned int iv0;
    fin >> iv0;
    assert( iv0 == iv );
    aVertex[iv].ReadFile(fin);
  }
  ////
  int ne;
  fin >> ne;
  aEdge.resize(ne);
  for(std::size_t ie=0;ie<aEdge.size();++ie){
    unsigned int ie0;
    fin >> ie0;
    assert( ie0 == ie );
    aEdge[ie].ReadFile(fin);
  }
  ////
  int nfc;
  fin >> nfc;
  aFace.resize(nfc);
  for(std::size_t ifc=0;ifc<aFace.size();++ifc){
    unsigned int ifc0;
    fin >> ifc0;
    assert( ifc0 == ifc );
    aFace[ifc].ReadFile(fin);
  }
  {
    for(auto & iv : aVertex){
      iv.isConst[0] = false;
      iv.isConst[1] = false;
      iv.isConst[2] = false;
    }
    for(auto & ie : aEdge){
      int iv0 = ie.iv0;
      int iv1 = ie.iv1;
      int inorm = ie.inorm;
      if( inorm < 0 || inorm >= 3 ){ continue; }
      aVertex[iv0].isConst[inorm] = true;
      aVertex[iv1].isConst[inorm] = true;
    }
  }
  for(auto & ie : aEdge){
    ie.Initialize(aVertex,elen); // ie0+0
  }
  for(auto & ifc : aFace){
    ifc.Initialize(aVertex,aEdge, elen); // ie0+0
  }
  BuildTriMesh(aXYZ,aTri,aTriSuTri,aNorm, aVertex,aEdge,aFace, isym);
}







