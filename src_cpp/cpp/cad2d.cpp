/*
 * Copyright (c) 2019 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */


#include <stdio.h>
#include <deque>
#include <set>

#if defined(__APPLE__) && defined(__MACH__)
#include <OpenGL/gl.h>
#elif defined(__MINGW32__) // probably I'm using Qt and don't want to use GLUT
#include <GL/glu.h>
#elif defined(_WIN32)
#include <windows.h>
#include <GL/glu.h>
#else
#include <GL/glu.h>
#endif

#include "delfem2/mshtopo.h"
#include "delfem2/msh.h"
#include "delfem2/bv.h"

#include "delfem2/dyntri.h"
#include "delfem2/dyntri_v2.h"

#include "delfem2/funcs_gl.h"
#include "delfem2/color_gl.h"
#include "delfem2/v23q_gl.h"

#include "delfem2/cad2d.h"


void CCad2D::Draw() const
{
  ::glDisable(GL_LIGHTING);
  ::glPointSize(6);
  ::glBegin(GL_POINTS);
  for(unsigned int iv=0;iv<aVtx.size();++iv){
    if( (int)iv == ivtx_picked ){ ::glColor3d(1,1,0); }
    else{ ::glColor3d(1,0,0); }
    ::glVertex3d( aVtx[iv].pos.x, aVtx[iv].pos.y, 0.0);
  }
  ::glEnd();
  /////
  ::glLineWidth(3);
  ::glBegin(GL_LINES);
  for(unsigned int ie=0;ie<aEdge.size();++ie){
    if( (int)ie == iedge_picked ){ ::glColor3d(1,1,0); }
    else{ ::glColor3d(0,0,0); }
    int iv0 = topo.aEdge[ie].iv0;
    int iv1 = topo.aEdge[ie].iv1;
    ::glVertex3d( aVtx[iv0].pos.x, aVtx[iv0].pos.y, -0.1);
    ::glVertex3d( aVtx[iv1].pos.x, aVtx[iv1].pos.y, -0.1);
  }
  ::glEnd();
  //////
  if( is_draw_face ){
    ::glLineWidth(1);
    glTranslated(0,0,-0.2);
    for(unsigned int iface=0;iface<aFace.size();++iface){
      const CCad2D_FaceGeo& face = aFace[iface];
      if( (int)iface == iface_picked ){ ::glColor3d(1,1,0); }
      else{ ::glColor3d(0.8,0.8,0.8); }
      DrawMeshTri2D_Face(face.aTri, face.aXY);
      ::glColor3d(0.0,0.0,0.0);
      DrawMeshTri2D_Edge(face.aTri, face.aXY);
    }
    glTranslated(0,0,+0.2);
  }
}

void CCad2D::Pick(double x0, double y0,
                  double view_height)
{
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

void CCad2D::Tessellation()
{
  for(unsigned int ie=0;ie<topo.aEdge.size();++ie){
    aEdge[ie].GenMesh(ie,topo,aVtx);
  }
  for(unsigned int ifc=0;ifc<topo.aFace.size();++ifc){
    aFace[ifc].GenMesh(ifc, topo, aEdge);
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
    int iv0 = topo.aEdge[iedge_picked].iv0;
    int iv1 = topo.aEdge[iedge_picked].iv1;
    aVtx[iv0].pos.x += p1x-p0x;
    aVtx[iv0].pos.y += p1y-p0y;
    aVtx[iv1].pos.x += p1x-p0x;
    aVtx[iv1].pos.y += p1y-p0y;
    Tessellation();
    return;
  }
  if( iface_picked >= 0 && iface_picked < (int)aFace.size() ){
    std::vector<int> aIdV = topo.aFace[iface_picked].GetArray_IdVertex(topo.aEdge);
    for(unsigned int iiv=0;iiv<aIdV.size();++iiv){
      const int iv1 = aIdV[iiv];
      aVtx[iv1].pos.x += p1x-p0x;
      aVtx[iv1].pos.y += p1y-p0y;
    }
    Tessellation();
    return;
  }
}

void CCad2D::Check() const
{
  this->topo.Check();
}

void CCad2D::AddPolygon(const std::vector<double>& aXY)
{
  const int np = aXY.size()/2;
  topo.AddPolygon(np);
  for(int ip=0;ip<np;++ip){
    aVtx.push_back(CVector2(aXY[ip*2+0], aXY[ip*2+1]));
  }
  ////
  const int iedge0 = aEdge.size();
  const int iface0 = aFace.size();
  for(int ie=0;ie<np;++ie){
    aEdge.push_back(CCad2D_EdgeGeo());
  }
  aFace.push_back(CCad2D_FaceGeo());
  ////
  for(int ie=0;ie<np;++ie){
    aEdge[iedge0+ie].GenMesh(iedge0+ie,topo,aVtx);
  }
  aFace[iface0].GenMesh(iface0, topo, aEdge);
}


void CCad2D::AddVtxEdge(double x, double y, int ie_add)
{
  if( ie_add < 0 || ie_add >= (int)topo.aEdge.size() ){ return; }
  topo.AddPoint_Edge(ie_add);
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

void CCad2D::Meshing
(std::vector<double>& aXY,
 std::vector<int>& aTri,
 std::vector<int>& aFlgPnt,
 std::vector<int>& aFlgTri,
 double elen) const
{
  aXY.clear();
  aTri.clear();
  if( aVtx.empty() ){ return; }
  std::vector<CVector2> aVec2;
  for(unsigned int iv=0;iv<aVtx.size();++iv){
    aVec2.push_back( aVtx[iv].pos );
    aFlgPnt.push_back(iv);
  }
  std::vector<int> edgeIP_ind,edgeIP;
  edgeIP_ind.push_back(0);
  for(unsigned int ie=0;ie<aEdge.size();++ie){
    std::vector<CVector2> aP0;
    aEdge[ie].GetInternalPoints_ElemLen(aP0, elen);
    for(unsigned int ip=0;ip<aP0.size();++ip){
      const int ip0 = aVec2.size();
      edgeIP.push_back(ip0);
      aVec2.push_back(aP0[ip]);
      aFlgPnt.push_back(ie+aVtx.size());
    }
    edgeIP_ind.push_back(edgeIP_ind[ie]+aP0.size());
  }
  assert( aVec2.size() == aFlgPnt.size() );
  /*
  for(int ie=0;ie<edgeIP_ind.size()-1;++ie){
    std::cout << ie << " --> ";
    for(int iip=edgeIP_ind[ie];iip<edgeIP_ind[ie+1];++iip){
      std::cout << edgeIP[iip] << " ";
    }
    std::cout << std::endl;
  }
   */
  std::vector<CEPo2> aPo2D;
  std::vector<ETri> aETri;
  Meshing_Initialize(aPo2D,aETri,aVec2);
  {
    aFlgPnt.push_back(-1);
    aFlgPnt.push_back(-1);
    aFlgPnt.push_back(-1);
    assert(aFlgPnt.size() == aVec2.size() );
    assert(aFlgPnt.size() == aPo2D.size() );
  }
  ////
  for(unsigned int ie=0;ie<aEdge.size();++ie){
    const int iv0 = topo.aEdge[ie].iv0;
    const int iv1 = topo.aEdge[ie].iv1;
    const int nseg = edgeIP_ind[ie+1]-edgeIP_ind[ie]+1;
    for(int iseg=0;iseg<nseg;++iseg){
      int ip0 = iv0; if( iseg != 0      ){ ip0 = edgeIP[edgeIP_ind[ie]+iseg-1]; }
      int ip1 = iv1; if( iseg != nseg-1 ){ ip1 = edgeIP[edgeIP_ind[ie]+iseg+0]; }
      EnforceEdge(aPo2D,aETri,
                  ip0,ip1,aVec2);
    }
  }
  //////////////////////////////////////
  // Make Flag for Triangles
  aFlgTri.assign(aETri.size(),-1);
  for(unsigned int ifc=0;ifc<aFace.size();++ifc){
    int ip0 = -1, ip1 = -1;
    {
      int ie = topo.aFace[ifc].aIE[0].first;
      bool dir = topo.aFace[ifc].aIE[0].second;
      ip0 = (dir) ? topo.aEdge[ie].iv0 : topo.aEdge[ie].iv1;
      ip1 = (dir) ? topo.aEdge[ie].iv1 : topo.aEdge[ie].iv0;
      if( edgeIP_ind[ie+1]-edgeIP_ind[ie] > 0 ){
        if(dir){ ip1 = edgeIP[edgeIP_ind[ie]]; }
        else{    ip1 = edgeIP[edgeIP_ind[ie+1]-1]; }
      }
    }
    assert( ip0 != -1 && ip1 != -1 );
    int itri0_ker, iedtri;
    FindEdge_LookAllTriangles(itri0_ker, iedtri,
                              ip0,ip1, aETri);
    assert(itri0_ker>=0&&itri0_ker<(int)aETri.size());
    FlagConnected(aFlgTri,
                  aETri, itri0_ker,ifc);
  }
  //////////////////////////////////////
  { // Delete Outer Triangles & Points
    assert(aFlgTri.size()==aETri.size());
    assert(aFlgPnt.size()==aVec2.size());
    assert(aFlgPnt.size()==aPo2D.size());
    DeleteTriFlag(aETri,aFlgTri,
                  -1);
    DeletePointsFlag(aVec2,aPo2D,aETri, aFlgPnt,
                     -1);
    assert(aFlgTri.size()==aETri.size());
    assert(aFlgPnt.size()==aVec2.size());
    assert(aFlgPnt.size()==aPo2D.size());
  }
  if( elen > 1.0e-10 ){
    CInputTriangulation_Uniform param(1.0);
    MeshingInside(aPo2D,aETri,aVec2, aFlgPnt,aFlgTri,
                  aVec2.size(), aVtx.size()+aEdge.size(),  elen, param);
  }
  assert(aFlgTri.size()==aETri.size());
  assert(aFlgPnt.size()==aVec2.size());
  assert(aFlgPnt.size()==aPo2D.size());
  MeshTri2D_Export(aXY,aTri, aVec2,aETri);
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

std::vector<double> CCad2D::XY_Vtx_Face
(int iface) const
{
  std::vector<double> aXY;
  std::vector<int> aIdV = topo.aFace[iface].GetArray_IdVertex(topo.aEdge);
  for(unsigned int iiv=0;iiv<aIdV.size();++iiv){
    const int iv = aIdV[iiv];
    aXY.push_back( aVtx[iv].pos.x );
    aXY.push_back( aVtx[iv].pos.y );
  }
  return aXY;
}

std::vector<int> CCad2D::Ind_Vtx_Face(int iface) const
{
  return topo.aFace[iface].GetArray_IdVertex(topo.aEdge);
}

std::vector<std::pair<int,bool> > CCad2D::Ind_Edge_Face(int iface) const
{
//  std::vector<int> aIdE;
  std::vector<std::pair<int,bool> > aIdE;
  for(unsigned int iie=0;iie<topo.aFace[iface].aIE.size();++iie){
    const int ie0 = topo.aFace[iface].aIE[iie].first;
    const bool dir0 = topo.aFace[iface].aIE[iie].second;
    aIdE.push_back( std::make_pair(ie0,dir0) );
  }
  return aIdE;
}

///////////////////////////////////////////////////////////

void CCad2D_EdgeGeo::GenMesh
(unsigned int iedge, const CCadTopo& topo,
 std::vector<CCad2D_VtxGeo>& aVtxGeo)
{
  assert( iedge<topo.aEdge.size() );
  const int iv0 = topo.aEdge[iedge].iv0;
  const int iv1 = topo.aEdge[iedge].iv1;
  assert(iv0 >= 0 && iv0 < (int)aVtxGeo.size());
  assert(iv1 >= 0 && iv1 < (int)aVtxGeo.size());
  this->p0 = aVtxGeo[iv0].pos;
  this->p1 = aVtxGeo[iv1].pos;
}

void CCad2D_EdgeGeo::GetInternalPoints_ElemLen
(std::vector<CVector2>& aV, double elen) const
{
  aV.clear();
  const int nadd = (int)( ::Distance(p0,p1) / elen);
  if( nadd == 0 ) return;
  for(int iadd=0;iadd<nadd;++iadd){
    double r2 = (double)(iadd+1)/(nadd+1);
    CVector2 v2 = (1-r2)*p0 + r2*p1;
    aV.push_back(v2);
  }
}

double CCad2D_EdgeGeo::Distance(double x, double y) const
{
  CVector2 pn = GetNearest_LineSeg_Point(CVector2(x,y),
                                         this->p0,this->p1);
  return ::Distance(pn,CVector2(x,y));
}

///////////////////////////////////////////////////////////

void CCad2D_FaceGeo::GenMesh
(unsigned int iface0, const CCadTopo& topo,
 std::vector<CCad2D_EdgeGeo>& aEdgeGeo)
{
  assert( iface0<topo.aFace.size() );
  const std::vector< std::pair<int,bool> >& aIE = topo.aFace[iface0].aIE;
  std::vector< std::vector<double> > aaXY;
  aaXY.resize(1);
  for(unsigned int iie=0;iie<aIE.size();++iie){
    const unsigned int ie0 = (unsigned int)aIE[iie].first;
    assert( ie0<topo.aEdge.size() );
    const bool dir0 = aIE[iie].second;
    const CCad2D_EdgeGeo& eg0 = aEdgeGeo[ie0];
    CVector2 p0 = (dir0) ? eg0.p0 : eg0.p1;
    aaXY[0].push_back(p0.x);
    aaXY[0].push_back(p0.y);
  }
  std::vector<int> loopIP_ind,loopIP;
  std::vector<CVector2> aVec2;
  {
    JArray_FromVecVec_XY(loopIP_ind,loopIP, aVec2,
                         aaXY);
    if( !CheckInputBoundaryForTriangulation(loopIP_ind,aVec2) ){
      std::cout << "loop invalid" << std::endl;
      return;
    }
    FixLoopOrientation(loopIP,
                       loopIP_ind,aVec2);
  }
  {
    std::vector<CEPo2> aPo2D;
    std::vector<ETri> aETri;
    Meshing_SingleConnectedShape2D(aPo2D, aVec2, aETri,
                                   loopIP_ind, loopIP);
    MeshTri2D_Export(aXY,aTri, aVec2,aETri);
  }
}

