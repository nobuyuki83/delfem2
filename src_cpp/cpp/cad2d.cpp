/*
 * Copyright (c) 2019 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */


#include <stdio.h>
#include <deque>
#include <set>

#include "delfem2/mshtopo.h"
#include "delfem2/msh.h"
#include "delfem2/bv.h"

#include "delfem2/vec2.h"

#include "delfem2/dyntri.h"
#include "delfem2/dyntri_v2.h"

#include "delfem2/cad2d.h"



////////////////////////////////////////////////////////////////////////////////////

void CCad2D::Pick(double x0, double y0,
                  double view_height)
{
  CVector2 pin(x0,y0);
  if( this->iedge_picked != -1 ){
    const CCad2D_EdgeGeo& edge = aEdge[iedge_picked];
    if( edge.type_edge == 1 ){
      assert( edge.param.size() == 4 );
      const CVector2 lx = (edge.p1 - edge.p0).Normalize();
      const CVector2 ly = CVector2(lx.y,-lx.x);
      const CVector2 q0 = edge.p0 + edge.param[0]*lx + edge.param[1]*ly;
      const CVector2 q1 = edge.p1 + edge.param[2]*lx + edge.param[3]*ly;
      if( Distance(pin, q0) < view_height*0.05 ){ this->ipicked_elem = 1; return; }
      if( Distance(pin, q1) < view_height*0.05 ){ this->ipicked_elem = 2; return; }
    }
  }
  this->ipicked_elem = 0;
  this->ivtx_picked = -1;
  this->iedge_picked = -1;
  this->iface_picked = -1;
  for(unsigned int ivtx=0;ivtx<aVtx.size();++ivtx){
    double x1 = aVtx[ivtx].pos.x;
    double y1 = aVtx[ivtx].pos.y;
    double dist = sqrt( (x0-x1)*(x0-x1)+(y0-y1)*(y0-y1) );
    if( dist < view_height*0.05 ){
      this->ivtx_picked = ivtx;
      return;
    }
  }
  ////
  for(unsigned int iedge=0;iedge<aEdge.size();++iedge){
    double dist = aEdge[iedge].Distance(x0, y0);
    if( dist < view_height*0.05 ){
      this->iedge_picked = iedge;
      return;
    }
  }
  ////
  for(unsigned int iface=0;iface<aFace.size();++iface){
    bool is_inside = aFace[iface].IsInside(x0, y0);
    if( is_inside ){
      this->iface_picked = iface;
      return;
    }
  }
}



void CCad2D::DragPicked(double p1x, double p1y, double p0x, double p0y)
{
  if( ivtx_picked >= 0 && ivtx_picked < (int)aVtx.size() ){
    aVtx[ivtx_picked].pos.x = p1x;
    aVtx[ivtx_picked].pos.y = p1y;
    Tessellation();
    return;
  }
  if( iedge_picked >= 0 && iedge_picked < (int)aEdge.size() ){
    if( ipicked_elem == 0 ){
      int iv0 = topo.aEdge[iedge_picked].iv0;
      int iv1 = topo.aEdge[iedge_picked].iv1;
      aVtx[iv0].pos.x += p1x-p0x;
      aVtx[iv0].pos.y += p1y-p0y;
      aVtx[iv1].pos.x += p1x-p0x;
      aVtx[iv1].pos.y += p1y-p0y;
    }
    else{
      CCad2D_EdgeGeo& edge = aEdge[iedge_picked];
      if( edge.type_edge == 1 ){
        assert( edge.param.size() == 4 );
        const CVector2 lx = (edge.p1 - edge.p0).Normalize();
        const CVector2 ly = CVector2(lx.y,-lx.x);
        if( ipicked_elem == 1 ){
          edge.param[0] = (CVector2(p1x,p1y)- edge.p0 )*lx;
          edge.param[1] = (CVector2(p1x,p1y)- edge.p0 )*ly;
        }
        else if( ipicked_elem == 2 ){
          edge.param[2] = (CVector2(p1x,p1y)- edge.p1 )*lx;
          edge.param[3] = (CVector2(p1x,p1y)- edge.p1 )*ly;
        }
      }
    }
    Tessellation();
    return;
  }
  if( iface_picked >= 0 && iface_picked < (int)aFace.size() ){
    std::vector<int> aIdV = topo.aLoop[iface_picked].GetArray_IdVertex(topo.aEdge);
    for(unsigned int iiv=0;iiv<aIdV.size();++iiv){
      const int iv1 = aIdV[iiv];
      aVtx[iv1].pos.x += p1x-p0x;
      aVtx[iv1].pos.y += p1y-p0y;
    }
    Tessellation();
    return;
  }
}

bool CCad2D::Check() const
{
  if( !this->topo.Check() ){ assert(0); return false; }
  if( aVtx.size() != topo.nVertex ){ assert(0); return false; }
  if( aEdge.size() != topo.aEdge.size() ){ assert(0); return false; }
  if( aFace.size() != topo.aFace.size() ){ assert(0); return false; }
  return true;
}

void CCad2D::AddPolygon(const std::vector<double>& aXY)
{
  const int np = aXY.size()/2;
  topo.AddPolygon(np);
  for(int ip=0;ip<np;++ip){
    aVtx.push_back(CVector2(aXY[ip*2+0], aXY[ip*2+1]));
  }
  ////
//  const unsigned int iedge0 = aEdge.size();
//  const unsigned int iface0 = aFace.size();
  for(int ie=0;ie<np;++ie){
    aEdge.push_back(CCad2D_EdgeGeo());
  }
  aFace.push_back(CCad2D_FaceGeo());
  ////
  assert( this->Check() );
  Tessellation();
}

void CCad2D::AddVtxFace(double x0, double y0, unsigned int ifc_add)
{
  if( ifc_add >= topo.aFace.size() ){ return; }
  topo.AddVtx_Face(ifc_add);
  assert( topo.Check() );
  aVtx.push_back( CCad2D_VtxGeo(CVector2(x0,y0)) );
  assert( this->Check() );
  Tessellation();
}

void CCad2D::AddVtxEdge(double x, double y, unsigned int ie_add)
{
  if( ie_add >= topo.aEdge.size() ){ return; }
  topo.AddVtx_Edge(ie_add);
  assert( topo.Check() );
  aVtx.push_back( CCad2D_VtxGeo(CVector2(x,y)) );
  aEdge.push_back( CCad2D_EdgeGeo() );
  Tessellation();
}

std::vector<double> CCad2D::MinMaxXYZ() const
{
  CBV3D_AABB aabb;
  for(unsigned int iv=0;iv<aVtx.size();++iv){
    const CVector2& v = aVtx[iv].pos;
    aabb.AddPoint(v.x, v.y, 0.0, 0.0);
  }
  return aabb.MinMaxXYZ();
}

void CCad2D::GetPointsEdge
(std::vector<int>& aIdP,
 const double* pXY, int np,
 const std::vector<int>& aIE,
 double tolerance ) const
{
  aIdP.clear();
  for(int ip=0;ip<np;++ip){
    const double x = pXY[ip*2+0];
    const double y = pXY[ip*2+1];
    for(unsigned int iie=0;iie<aIE.size();++iie){
      const int ie0 = aIE[iie];
      const CCad2D_EdgeGeo& eg = this->aEdge[ie0];
      const double dist = eg.Distance(x,y);
      if( dist > tolerance ){ continue; }
      aIdP.push_back(ip);
    }
  }
}

std::vector<double> CCad2D::XY_VtxCtrl_Face
(int iface) const
{
  std::vector<double> aXY;
  for(unsigned int iie=0;iie<topo.aLoop[iface].aIE.size();++iie){
    int ie0 = topo.aLoop[iface].aIE[iie].first;
    bool dir = topo.aLoop[iface].aIE[iie].second;
    if( dir ){
      const int iv = topo.aEdge[ie0].iv0;
      aXY.push_back( aVtx[iv].pos.x );
      aXY.push_back( aVtx[iv].pos.y );
      if( aEdge[ie0].type_edge == 1 ){
        const CCad2D_EdgeGeo& edge = aEdge[ie0];
        const CVector2 lx = (edge.p1 - edge.p0).Normalize();
        const CVector2 ly = CVector2(lx.y,-lx.x);
        const CVector2 q0 = edge.p0 + edge.param[0]*lx + edge.param[1]*ly;
        const CVector2 q1 = edge.p1 + edge.param[2]*lx + edge.param[3]*ly;
        aXY.push_back( q0.x );
        aXY.push_back( q0.y );
        aXY.push_back( q1.x );
        aXY.push_back( q1.y );
      }
    }
    else{
      const int iv = topo.aEdge[ie0].iv1;
      aXY.push_back( aVtx[iv].pos.x );
      aXY.push_back( aVtx[iv].pos.y );
    }
  }
  /*
  std::vector<int> aIdV = topo.aFace[iface].GetArray_IdVertex(topo.aEdge);
  for(unsigned int iiv=0;iiv<aIdV.size();++iiv){
    const int iv = aIdV[iiv];
    aXY.push_back( aVtx[iv].pos.x );
    aXY.push_back( aVtx[iv].pos.y );
  }
   */
  return aXY;
}

std::vector<int> CCad2D::Ind_Vtx_Face(int iface) const
{
  return topo.aLoop[iface].GetArray_IdVertex(topo.aEdge);
}

std::vector<int> CCad2D::Ind_Vtx_Edge(int iedge) const
{
  std::vector<int> aRes;
  if( iedge < 0 || iedge > (int)topo.aEdge.size() ){ return aRes; }
  aRes.push_back(topo.aEdge[iedge].iv0);
  aRes.push_back(topo.aEdge[iedge].iv1);
  return aRes;
}

std::vector<std::pair<int,bool> > CCad2D::Ind_Edge_Face(int iface) const
{
//  std::vector<int> aIdE;
  std::vector<std::pair<int,bool> > aIdE;
  for(unsigned int iie=0;iie<topo.aLoop[iface].aIE.size();++iie){
    const int ie0 = topo.aLoop[iface].aIE[iie].first;
    const bool dir0 = topo.aLoop[iface].aIE[iie].second;
    aIdE.push_back( std::make_pair(ie0,dir0) );
  }
  return aIdE;
}

///////////////////////////////////////////////////////////

// for visualization
void CCad2D_EdgeGeo::GenMesh
(unsigned int iedge, const CCadTopo& topo,
 const std::vector<CCad2D_VtxGeo>& aVtxGeo,
 double elen)
{
  assert( iedge<topo.aEdge.size() );
  const int iv0 = topo.aEdge[iedge].iv0;
  const int iv1 = topo.aEdge[iedge].iv1;
  assert(iv0 >= 0 && iv0 < (int)aVtxGeo.size());
  assert(iv1 >= 0 && iv1 < (int)aVtxGeo.size());
  this->p0 = aVtxGeo[iv0].pos;
  this->p1 = aVtxGeo[iv1].pos;
  aP.clear();
  if( type_edge == 0 ){
    if( elen < 0 ){ return; }
    const int nadd = (int)( this->Length() / elen);
    for(int iadd=0;iadd<nadd;++iadd){
      double r2 = (double)(iadd+1)/(nadd+1);
      CVector2 v2 = (1-r2)*p0 + r2*p1;
      aP.push_back(v2);
    }
  }
  if( this->type_edge == 1 ){
    const CVector2 lx = (this->p1 - this->p0).Normalize();
    const CVector2 ly = CVector2(lx.y,-lx.x);
    const CVector2 q0 = p0 + param[0]*lx + param[1]*ly;
    const CVector2 q1 = p1 + param[2]*lx + param[3]*ly;
    int ndiv = 10;
    if( elen > 0 ){
      std::vector<CVector2> aP0;
      Polyline_BezierCubic(aP0,
                           20, p0,q0,q1,p1);
      const double len0 = Length_Polygon(aP0);
      ndiv = len0 / elen;
    }
    for(int ip=1;ip<ndiv;++ip){
      double t = (double)ip/ndiv;
      CVector2 pos = pointCurve_BezierCubic(t, p0, q0, q1, p1);
      aP.push_back(pos);
    }
  }
}

/*
// for meshing
void CCad2D_EdgeGeo::GetInternalPoints_ElemLen
(std::vector<CVector2>& aV, double elen) const
{
  aV.clear();
  const int nadd = (int)( this->Length() / elen);
  if( nadd == 0 ) return;
  /////
  if( type_edge == 0 ){
    for(int iadd=0;iadd<nadd;++iadd){
      double r2 = (double)(iadd+1)/(nadd+1);
      CVector2 v2 = (1-r2)*p0 + r2*p1;
      aV.push_back(v2);
    }
  }
  else if( type_edge == 1 ){
    const CVector2 lx = (this->p1 - this->p0).Normalize();
    const CVector2 ly = CVector2(lx.y,-lx.x);
    const CVector2 q0 = p0 + param[0]*lx + param[1]*ly;
    const CVector2 q1 = p1 + param[2]*lx + param[3]*ly;
    for(int iadd=0;iadd<nadd;++iadd){
      double t = (double)(iadd+1)/(nadd+1);
      CVector2 pos = pointCurve_BezierCubic(t, p0, q0, q1, p1);
      aV.push_back(pos);
    }
  }
}
 */

double CCad2D_EdgeGeo::Distance(double x, double y) const
{
  const CVector2 q(x,y);
  if( type_edge == 0 ){
    const CVector2 pn = GetNearest_LineSeg_Point(q,p0,p1);
    return ::Distance(pn,q);
  }
  else if( type_edge == 1 ){
    assert( param.size() == 4 );
    double min_dist = -1;
    for(unsigned int ie=0;ie<aP.size()+1;++ie){
      CVector2 q0 = (ie==0) ? p0 : aP[ie];
      CVector2 q1 = (ie==aP.size()-1) ? p1 : aP[ie+1];
      double dist = ::Distance(q, GetNearest_LineSeg_Point(q,q0,q1));
      if( min_dist < 0 || dist < min_dist ){ min_dist = dist; }
    }
    return min_dist;
  }
  assert(0);
  return 0;
}

double CCad2D_EdgeGeo::Length() const
{
  double len0 = 0.0;
  for(unsigned int ie=0;ie<aP.size()+1;++ie){
    const CVector2 q0 = (ie==0) ? p0 : aP[ie-1];
    const CVector2 q1 = (ie==aP.size()) ? p1 : aP[ie];
    double dist0 = ::Distance(q0, q1);
    len0 += dist0;
  }
  return len0;
}




///////////////////////////////////////////////////////////

void GetBound
(double bound_2d[4],
 unsigned int ifc0,
 const CCadTopo& topo,
 const std::vector<CCad2D_EdgeGeo>& aEdgeGeo)
{
  assert( ifc0 < topo.aFace.size() );
  const int il0 = topo.aFace[ifc0].aIL[0];
  const std::vector< std::pair<int,bool> >& aIE = topo.aLoop[il0].aIE;
  {
    int ie0 = aIE[0].first;
    const CVector2& p0 = aEdgeGeo[ie0].p0;
    bound_2d[0] = p0.x;
    bound_2d[1] = p0.x;
    bound_2d[2] = p0.y;
    bound_2d[3] = p0.y;
  }
  for(unsigned int iie=0;iie<aIE.size();++iie){
    const unsigned int ie0 = aIE[iie].first;
    {
      const CVector2& p0 = aEdgeGeo[ie0].p0;
      if( p0.x < bound_2d[0] ){ bound_2d[0] = p0.x; }
      if( p0.x > bound_2d[1] ){ bound_2d[1] = p0.x; }
      if( p0.y < bound_2d[2] ){ bound_2d[2] = p0.y; }
      if( p0.y > bound_2d[3] ){ bound_2d[3] = p0.y; }
    }
    for(unsigned int ip=0;ip<aEdgeGeo[ie0].aP.size();ip++){
      const CVector2& p0 = aEdgeGeo[ie0].aP[ip];
      if( p0.x < bound_2d[0] ){ bound_2d[0] = p0.x; }
      if( p0.x > bound_2d[1] ){ bound_2d[1] = p0.x; }
      if( p0.y < bound_2d[2] ){ bound_2d[2] = p0.y; }
      if( p0.y > bound_2d[3] ){ bound_2d[3] = p0.y; }
    }
  }
}

void GenMeshCadFace
(std::vector<CVector2>& aVec2,
 std::vector<ETri>& aETri,
 const CCad2D_FaceGeo& facegeo, unsigned int iface0,
 const CCadTopo& topo,
 const std::vector<CCad2D_VtxGeo>& aVtxGeo,
 const std::vector<CCad2D_EdgeGeo>& aEdgeGeo)
{
  assert( iface0<topo.aFace.size() );
  std::vector<CEPo2> aPo2D;
  {
    aPo2D.resize(aVec2.size());
    for(unsigned int ixys=0;ixys<aVec2.size();ixys++){
      aPo2D[ixys].e = -1;
      aPo2D[ixys].d = -1;
    }
  }
  {
    double bound_2d[4];
    GetBound(bound_2d,
             iface0, topo, aEdgeGeo);
    MakeSuperTriangle(aVec2, aPo2D, aETri,
                      bound_2d);
  }
  std::vector<int> aIP;
  { // make list of index of point involved this face
    for(int iil=0;iil<topo.aFace[iface0].aIL.size();++iil){
      const int il0 = topo.aFace[iface0].aIL[iil];
      if( topo.aLoop[il0].iv != -1 ){
        aIP.push_back(topo.aLoop [il0].iv );
        continue;
      }
      const std::vector< std::pair<int,bool> >& aIE = topo.aLoop[il0].aIE;
      for(unsigned int iie=0;iie<aIE.size();++iie){
        const int ie = aIE[iie].first;
        aIP.push_back(topo.aEdge[ie].iv0);
        for(unsigned int iip=0;iip<aEdgeGeo[ie].aP.size();++iip){
          aIP.push_back(aEdgeGeo[ie].ip0+iip);
        }
      }
    }
  }
  {
    const double MIN_TRI_AREA = 1.0e-10;
    for(unsigned int iip=0;iip<aIP.size();++iip){
      const int ip = aIP[iip];
      AddPointsMesh(aVec2,aPo2D,aETri,
                    ip,MIN_TRI_AREA);
      DelaunayAroundPoint(ip,aPo2D,aETri,aVec2);
    }
  }
  {
    for(unsigned int iil=0;iil<topo.aFace[iface0].aIL.size();++iil){
      const int il0 = topo.aFace[iface0].aIL[iil];
      if( topo.aLoop[il0].iv != -1 ){ continue; }
      const std::vector< std::pair<int,bool> >& aIE = topo.aLoop[il0].aIE;
      for(unsigned int iie=0;iie<aIE.size();++iie){
        const int ie0 = aIE[iie].first;
        const unsigned int np = aEdgeGeo[ie0].aP.size();
        const unsigned int nseg = np+1;
        for(unsigned int iseg=0;iseg<nseg;++iseg){
          const int ip0 = ( iseg == 0      ) ? topo.aEdge[ie0].iv0 : aEdgeGeo[ie0].ip0+iseg-1;
          const int ip1 = ( iseg == nseg-1 ) ? topo.aEdge[ie0].iv1 : aEdgeGeo[ie0].ip0+iseg;
          EnforceEdge(aPo2D,aETri,
                      ip0,ip1,aVec2);
        }
      }
    }
  }
  std::vector<int> aFlgTri(aETri.size(),-1);
  {
    int ip0 = -1, ip1 = -1;
    {
      const int il0 = topo.aFace[iface0].aIL[0];
      const std::vector< std::pair<int,bool> >& aIE = topo.aLoop[il0].aIE;
      unsigned int ie0 = aIE[0].first;
      const int np = aEdgeGeo[ie0].aP.size();
      ip0 = topo.aEdge[ie0].iv0;
      ip1 = ( np == 0 ) ? topo.aEdge[ie0].iv1 : aEdgeGeo[ie0].ip0;
    }
    assert( ip0 != -1 && ip1 != -1 );
    int itri0_ker, iedtri;
    FindEdge_LookAllTriangles(itri0_ker, iedtri,
                              ip0,ip1, aETri);
    assert(itri0_ker>=0&&itri0_ker<(int)aETri.size());
    FlagConnected(aFlgTri,
                  aETri, itri0_ker,iface0);
  }
  { // Delete Outer Triangles & Points
    assert(aFlgTri.size()==aETri.size());
    DeleteTriFlag(aETri,aFlgTri,
                  -1);
    aPo2D.resize(aPo2D.size()-3);
    aVec2.resize(aVec2.size()-3);
  }
}


void CCad2D::Tessellation()
{
  for(unsigned int ie=0;ie<topo.aEdge.size();++ie){
    aEdge[ie].GenMesh(ie,topo,aVtx,-1);
  }
  std::vector<CVector2> aVec2;
  for(unsigned int iv=0;iv<aVtx.size();++iv){
    aVec2.push_back( aVtx[iv].pos );
  }
  for(unsigned int ie=0;ie<aEdge.size();++ie){
    aEdge[ie].ip0 = aVec2.size();
    for(unsigned int ip=0;ip<aEdge[ie].aP.size();++ip){
      aVec2.push_back(aEdge[ie].aP[ip]);
    }
  }
  for(unsigned int ifc=0;ifc<topo.aFace.size();++ifc){
    std::vector<ETri> aETri;
    GenMeshCadFace(aVec2, aETri,
                   aFace[ifc],ifc,
                   topo,
                   aVtx,aEdge);
    MeshTri2D_Export(aFace[ifc].aXY, aFace[ifc].aTri,
                     aVec2, aETri);
    
  }
}

/////////////////////////////////////////////////////////////



void CMesher_Cad2D::Meshing
(CMeshDynTri2D& dmsh,
 const CCad2D& cad)
{
  std::vector<CCad2D_EdgeGeo> aEdgeGeo = cad.aEdge;
  for(unsigned int ie=0;ie<aEdgeGeo.size();++ie){
    aEdgeGeo[ie].GenMesh(ie, cad.topo, cad.aVtx,
                         edge_length);
  }
  ////
  aFlgPnt.clear();
  dmsh.Clear();
  for(unsigned int iv=0;iv<cad.aVtx.size();++iv){
    dmsh.aVec2.push_back( cad.aVtx[iv].pos );
    aFlgPnt.push_back(iv);
  }
  for(unsigned int ie=0;ie<aEdgeGeo.size();++ie){
    aEdgeGeo[ie].ip0 = dmsh.aVec2.size();
    for(unsigned int ip=0;ip<aEdgeGeo[ie].aP.size();++ip){
      dmsh.aVec2.push_back(aEdgeGeo[ie].aP[ip]);
      aFlgPnt.push_back(cad.aVtx.size()+ie);
    }
  }
  ///////
  aFlgTri.clear();
  { // add inital face
    { // face 0
      GenMeshCadFace(dmsh.aVec2, dmsh.aETri,
                     cad.aFace[0], 0,
                     cad.topo,cad.aVtx,aEdgeGeo);
      aFlgTri.resize(dmsh.aETri.size(),0);
    }
    for(unsigned int ifc=1;ifc<cad.aFace.size();++ifc){ // face index bigger than 0
      std::vector<ETri> aETri;
      GenMeshCadFace(dmsh.aVec2, aETri,
                     cad.aFace[ifc], ifc,
                     cad.topo,cad.aVtx,aEdgeGeo);
      const unsigned int ntri0 = dmsh.aETri.size();
      for(unsigned int it=0;it<aETri.size();++it){
        if( aETri[it].s2[0]>=0){ aETri[it].s2[0] += ntri0; }
        if( aETri[it].s2[1]>=0){ aETri[it].s2[1] += ntri0; }
        if( aETri[it].s2[2]>=0){ aETri[it].s2[2] += ntri0; }
        dmsh.aETri.push_back(aETri[it]);
      }
      aFlgTri.resize(dmsh.aETri.size(),ifc);
    }
  }
  { // make EPo
    dmsh.aEPo.resize(dmsh.aVec2.size());
    for(int it=0;it<dmsh.aETri.size();++it){
      for(int inotri=0;inotri<3;++inotri){
        dmsh.aEPo[ dmsh.aETri[it].v[inotri] ].e = it;
        dmsh.aEPo[ dmsh.aETri[it].v[inotri] ].d = inotri;
      }
    }
  }
  if( edge_length > 1.0e-10 ){
    CInputTriangulation_Uniform param(1.0);
    MeshingInside(dmsh.aEPo,dmsh.aETri,dmsh.aVec2,
                  aFlgPnt,aFlgTri,
                  dmsh.aVec2.size(), cad.aVtx.size()+cad.aEdge.size(),
                  edge_length, param);
  }
  nvtx = cad.aVtx.size();
  nedge = cad.aEdge.size();
  nface = cad.aFace.size();
}


std::vector<int> CMesher_Cad2D::IndPoint_IndEdge
(const int iedge,
 bool is_end_point,
 const CCad2D& cad2d)
{
  //    std::cout << nvtx << " " << nedge << " " << nface << std::endl;
  std::vector<int> aflg(nvtx+nedge+nface,0);
  {
    aflg[nvtx+iedge] = 1;
  }
  std::vector<int> aIP_E = cad2d.Ind_Vtx_Edge(iedge);
  std::vector<int> res;
  if( is_end_point ){ res.push_back(aIP_E[0]); }
  for(unsigned int ip=0;ip<this->aFlgPnt.size();++ip){
    int iflg = aFlgPnt[ip]; assert(iflg<int(nvtx+nedge+nface));
    if( iflg >= (int)(nvtx+nedge) ){ break; }
    if( aflg[iflg] == 1 ){ res.push_back(ip); }
  }
  if( is_end_point ){ res.push_back(aIP_E[1]); }
  return res;
}

std::vector<int> CMesher_Cad2D::IndPoint_IndEdgeArray
(const std::vector<int>& aIndEd,
 const CCad2D& cad2d)
{
  //    std::cout << nvtx << " " << nedge << " " << nface << std::endl;
  std::vector<int> aflg(nvtx+nedge+nface,0);
  {
    for(unsigned int iie=0;iie<aIndEd.size();++iie){
      const unsigned int iedge = aIndEd[iie]; assert(iedge<nedge);
      aflg[iedge+nvtx] = 1;
      {
        const std::vector<int>& aIV = cad2d.Ind_Vtx_Edge(iedge);
        for(unsigned int iiv=0;iiv<aIV.size();++iiv){
          const int iv0 = aIV[iiv];
          aflg[iv0] = 1;
        }
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


std::vector<int> CMesher_Cad2D::IndPoint_IndFaceArray
(const std::vector<int>& aIndFc,
 const CCad2D& cad2d)
{
  std::vector<int> aflg(nvtx+nedge+nface,0);
  {
    for(unsigned int iif=0;iif<aIndFc.size();++iif){
      const unsigned int iface = aIndFc[iif]; assert(iface<nface);
      aflg[nvtx+nedge+iface] = 1;
      {
        const std::vector<std::pair<int,bool> >& aIE = cad2d.Ind_Edge_Face(iface);
        for(unsigned int iie=0;iie<aIE.size();++iie){
          const int ie0 = aIE[iie].first;
          aflg[nvtx+ie0] = 1;
        }
      }
      {
        const std::vector<int>& aIV = cad2d.Ind_Vtx_Face(iface);
        for(unsigned int iiv=0;iiv<aIV.size();++iiv){
          const int iv0 = aIV[iiv];
          aflg[iv0] = 1;
        }
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
