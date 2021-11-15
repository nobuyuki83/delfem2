/*
 * Copyright (c) 2019 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */


#include "delfem2/cad2_mesher.h"

#include <deque>
#include <climits>

#include "delfem2/dtri2_v2dtri.h"
#include "delfem2/geo_polyline2.h"
#include "delfem2/str.h"
#include "delfem2/cad2.h"

// ------------------------------------------

namespace delfem2::cad2_mesher {

DFM2_INLINE delfem2::CBoundingBox2<double> BB_LoopEdgeCad2D(
    const std::vector<CCad2D_EdgeGeo> &aEdge) {
  CBoundingBox2<double> bb;
  for (const auto &ie : aEdge) {
    CBoundingBox2<double> bb0 = ie.BB();
    bb += bb0;
  }
  return bb;
}

DFM2_INLINE void GetBound(
    double bound_2d[4],
    unsigned int ifc0,
    const CadTopo &topo,
    const std::vector<CCad2D_EdgeGeo> &aEdgeGeo) {
  assert(ifc0 < topo.faces.size());
  const int il0 = topo.faces[ifc0].aIL[0];
  const std::vector<std::pair<unsigned int, bool> > &aIE = topo.loops[il0].aIE;
  {
    int ie0 = aIE[0].first;
    const CVec2d &p0 = aEdgeGeo[ie0].p0;
    bound_2d[0] = p0.x;
    bound_2d[1] = p0.x;
    bound_2d[2] = p0.y;
    bound_2d[3] = p0.y;
  }
  for (const auto &iie : aIE) {
    const unsigned int ie0 = iie.first;
    {
      const CVec2d &p0 = aEdgeGeo[ie0].p0;
      if (p0.x < bound_2d[0]) { bound_2d[0] = p0.x; }
      if (p0.x > bound_2d[1]) { bound_2d[1] = p0.x; }
      if (p0.y < bound_2d[2]) { bound_2d[2] = p0.y; }
      if (p0.y > bound_2d[3]) { bound_2d[3] = p0.y; }
    }
    /*
    for (const auto &p0 : aEdgeGeo[ie0].aP) {
      if (p0.x < bound_2d[0]) { bound_2d[0] = p0.x; }
      if (p0.x > bound_2d[1]) { bound_2d[1] = p0.x; }
      if (p0.y < bound_2d[2]) { bound_2d[2] = p0.y; }
      if (p0.y > bound_2d[3]) { bound_2d[3] = p0.y; }
    }
     */
  }
}

DFM2_INLINE void GenMeshCadFace(
	std::vector<CVec2d> &aVec2,
    std::vector<CDynTri> &aETri,
	unsigned int iface0,
    const CadTopo &topo,
    [[maybe_unused]] const std::vector<CCad2D_VtxGeo> &aVtxGeo,
    [[maybe_unused]] const std::vector<CCad2D_EdgeGeo> &aEdgeGeo,
    const std::vector< std::vector<std::pair<unsigned int, double> > >& edge_point){
  assert(iface0 < topo.faces.size());
  std::vector<CDynPntSur> aPo2D;
  {
    aPo2D.resize(aVec2.size());
    for (size_t ixys = 0; ixys < aVec2.size(); ixys++) {
      aPo2D[ixys].e = UINT_MAX;
      aPo2D[ixys].d = 0;
    }
  }
  {
    double bound_2d[4];
    GetBound(
        bound_2d,
        iface0, topo, aEdgeGeo);
    MakeSuperTriangle(
        aVec2, aPo2D, aETri,
        bound_2d);
  }
  std::vector<int> aIP;
  { // make list of index of point involved this face
    for (size_t iil = 0; iil < topo.faces[iface0].aIL.size(); ++iil) {
      const int il0 = topo.faces[iface0].aIL[iil];
      if (topo.loops[il0].iv != UINT_MAX) {
        aIP.push_back(topo.loops[il0].iv);
        continue;
      }
      const std::vector<std::pair<unsigned int, bool> > &aIE = topo.loops[il0].aIE;
      for (const auto &iie : aIE) {
        const int ie = iie.first;
        aIP.push_back(topo.edges[ie].iv0);
        for( const auto& pair : edge_point[ie] ) {
          aIP.push_back( pair.first );
        }
      }
    }
  }
  {
    const double MIN_TRI_AREA = 1.0e-10;
    for (int ip : aIP) {
      AddPointsMesh(aVec2, aPo2D, aETri,
                    ip, MIN_TRI_AREA);
      DelaunayAroundPoint(ip, aPo2D, aETri, aVec2);
    }
  }
  {
    for (size_t iil = 0; iil < topo.faces[iface0].aIL.size(); ++iil) {
      const int il0 = topo.faces[iface0].aIL[iil];
      if (topo.loops[il0].iv != UINT_MAX) { continue; }
      const std::vector<std::pair<unsigned int, bool> > &aIE = topo.loops[il0].aIE;
      for (const auto &iie : aIE) {
        const int ie0 = iie.first;
        const size_t np = edge_point[ie0].size();
        const size_t nseg = np + 1;
        for (unsigned int iseg = 0; iseg < nseg; ++iseg) {
          const int ip0 = (iseg == 0) ? topo.edges[ie0].iv0 : edge_point[ie0][iseg-1].first;
          const int ip1 = (iseg == nseg - 1) ? topo.edges[ie0].iv1 : edge_point[ie0][iseg].first;
          EnforceEdge(aPo2D, aETri,
                      ip0, ip1, aVec2);
        }
      }
    }
  }
  std::vector<int> aFlgTri(aETri.size(), -1);
  {
    int ip0 = -1, ip1 = -1;
    {
      const int il0 = topo.faces[iface0].aIL[0];
      const std::vector<std::pair<unsigned int, bool> > &aIE = topo.loops[il0].aIE;
      unsigned int ie0 = aIE[0].first;
      const size_t np = edge_point[ie0].size();
      ip0 = topo.edges[ie0].iv0;
      ip1 = (np == 0) ? topo.edges[ie0].iv1 : edge_point[ie0][0].first;
    }
    assert(ip0 != -1 && ip1 != -1);
    unsigned int itri0_ker=UINT_MAX, iedtri;
    FindEdge_LookAllTriangles(itri0_ker, iedtri,
                              ip0, ip1, aETri);
    assert( itri0_ker < aETri.size() );
    FlagConnected(aFlgTri,
                  aETri, itri0_ker, (int) iface0);
  }
  { // Delete Outer Triangles & Points
    assert(aFlgTri.size() == aETri.size());
    DeleteTriFlag(aETri, aFlgTri,
                  -1);
    aPo2D.resize(aPo2D.size() - 3);
    aVec2.resize(aVec2.size() - 3);
  }
}

} // delfem2::cad2

// =========================================================

DFM2_INLINE void delfem2::CMesher_Cad2D::Meshing(
	CMeshDynTri2D& dmsh,
	const CCad2D& cad)
{
  std::vector< std::vector<double> > aaXYW(cad.nEdge());
  for(unsigned int ie=0;ie<cad.aEdge.size();++ie){
    unsigned int ndiv = 1;
    if( edge_length > 0 ) {
      if (this->mapIdEd_NDiv.find(ie) == this->mapIdEd_NDiv.end()) {
        const double len0 = cad.aEdge[ie].LengthNDiv(20);
        ndiv = static_cast<unsigned int>(len0 / edge_length + 1);
      } else {
        ndiv = this->mapIdEd_NDiv.find(ie)->second;
      }
    }
    else{
      ndiv = 20;
    }
    aaXYW[ie] = cad.aEdge[ie].GenMesh(ndiv);
  }
  //
  this->edge_point.resize(cad.aEdge.size());
  aFlgPnt.clear();
  dmsh.Clear();
  for(unsigned int iv=0;iv<cad.aVtx.size();++iv){
    dmsh.aVec2.push_back( cad.aVtx[iv].pos );
    aFlgPnt.push_back(iv);
  }
  for(unsigned int ie=0;ie<cad.aEdge.size();++ie){
    edge_point[ie].clear();
    for(unsigned int ip=0;ip<aaXYW[ie].size()/3;++ip){
      const unsigned int iv = dmsh.aVec2.size();  //  id of mesh vtx
      dmsh.aVec2.push_back({aaXYW[ie][ip*3+0],aaXYW[ie][ip*3+1]});
      edge_point[ie].push_back( std::make_pair(iv,aaXYW[ie][ip*3+2]) );
      aFlgPnt.push_back(static_cast<unsigned int>(cad.aVtx.size()+ie));
    }
  }
  //
  aFlgTri.clear();
  for(unsigned int ifc=0;ifc<cad.topo.faces.size();++ifc){ // add face to dmsh
    std::vector<CDynTri> aDTri;
    cad2_mesher::GenMeshCadFace(dmsh.aVec2, aDTri,
        ifc,
        cad.topo, cad.aVtx, cad.aEdge, edge_point);
    const unsigned int ntri0 = static_cast<unsigned int>(dmsh.aETri.size());
    for(auto & tri : aDTri){
      if( tri.s2[0]!=UINT_MAX ){ tri.s2[0] += ntri0; }
      if( tri.s2[1]!=UINT_MAX ){ tri.s2[1] += ntri0; }
      if( tri.s2[2]!=UINT_MAX ){ tri.s2[2] += ntri0; }
      dmsh.aETri.push_back(tri);
    }
    aFlgTri.resize(dmsh.aETri.size(),ifc);
  }
  { // make EPo
    dmsh.aEPo.resize(dmsh.aVec2.size());
    for(unsigned int it=0;it<dmsh.aETri.size();++it){
      for(unsigned int inotri=0;inotri<3;++inotri){
        dmsh.aEPo[ dmsh.aETri[it].v[inotri] ].e = it;
        dmsh.aEPo[ dmsh.aETri[it].v[inotri] ].d = inotri;
      }
    }
  }
  if( edge_length > 1.0e-10 ){
    CInputTriangulation_Uniform param(1.0);
    MeshingInside(
        dmsh.aEPo,dmsh.aETri,dmsh.aVec2,
        aFlgPnt,aFlgTri,
        dmsh.aVec2.size(),
        static_cast<unsigned int>(cad.aVtx.size()+cad.aEdge.size()),
        edge_length, param);
  }
  nvtx = cad.aVtx.size();
  nedge = cad.aEdge.size();
  nface = cad.topo.faces.size();
}


DFM2_INLINE std::vector<unsigned int>
    delfem2::CMesher_Cad2D::IndPoint_IndEdge(
    const unsigned int iedge,
    bool is_end_point,
    const CCad2D& cad2d)
{
  std::vector<unsigned int> res;
  if( iedge >= cad2d.aEdge.size() ){ return res; }
  std::vector<int> aflg(nvtx+nedge+nface,0);
  {
    aflg[nvtx+iedge] = 1;
  }
  std::pair<unsigned int, unsigned int> aIP_E = cad2d.topo.VertexIndexs_EdgeEndPoints(iedge);
  if( is_end_point ){ res.push_back(aIP_E.first); }
  for(unsigned int ip=0;ip<this->aFlgPnt.size();++ip){
    int iflg = aFlgPnt[ip]; assert(iflg<int(nvtx+nedge+nface));
    if( iflg >= (int)(nvtx+nedge) ){ break; }
    if( aflg[iflg] == 1 ){ res.push_back(ip); }
  }
  if( is_end_point ){ res.push_back(aIP_E.second); }
  return res;
}

DFM2_INLINE std::vector<unsigned int>
    delfem2::CMesher_Cad2D::IndPoint_IndEdgeArray(
    const std::vector<int>& aIndEd,
    const CCad2D& cad2d)
{
  //    std::cout << nvtx << " " << nedge << " " << nface << std::endl;
  std::vector<int> aflg(nvtx+nedge+nface,0);
  {
    for(unsigned int iedge : aIndEd){
      assert(iedge<nedge);
      aflg[iedge+nvtx] = 1;
      {
        const auto pair = cad2d.topo.VertexIndexs_EdgeEndPoints(iedge);
        aflg[pair.first] = 1;
        aflg[pair.second] = 1;
      }
    }
  }
  std::vector<unsigned int> res;
  for(unsigned int ip=0;ip<this->aFlgPnt.size();++ip){
    int iflg = aFlgPnt[ip]; assert(iflg<int(nvtx+nedge+nface));
    if( aflg[iflg] == 1 ){ res.push_back(ip); }
  }
  return res;
}


DFM2_INLINE std::vector<int> delfem2::CMesher_Cad2D::IndPoint_IndFaceArray(
    const std::vector<int>& aIndFc,
    const CCad2D& cad2d)
{
  std::vector<int> aflg(nvtx+nedge+nface,0);
  {
    for(unsigned int iface : aIndFc){
      assert(iface<nface);
      aflg[nvtx+nedge+iface] = 1;
      for(const auto & iie : cad2d.topo.EdgeIndexes_Face(iface)){
        aflg[nvtx+iie.first] = 1;
      }
      for(unsigned int iv0 : cad2d.topo.VertexIndexs_Face(iface)){
        aflg[iv0] = 1;
      }
    }
  }
  std::vector<int> res;
  for(unsigned int ip=0;ip<this->aFlgPnt.size();++ip){
    int iflg = aFlgPnt[ip]; assert(iflg<int(nvtx+nedge+nface));
    if( aflg[iflg] == 1 ){ res.push_back(ip); }
  }
  return res;
}
