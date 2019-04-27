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
  ::glColor3d(0,0,0);
  ::glLineWidth(3);
  ::glBegin(GL_LINES);
  for(unsigned int ie=0;ie<aEdge.size();++ie){
    int iv0 = topo.aEdge[ie].iv0;
    int iv1 = topo.aEdge[ie].iv1;
    ::glVertex3d( aVtx[iv0].pos.x, aVtx[iv0].pos.y, -0.1);
    ::glVertex3d( aVtx[iv1].pos.x, aVtx[iv1].pos.y, -0.1);
  }
  ::glEnd();
  //////
  if( is_draw_face ){
    ::glColor3d(0.8,0.8,0.8);
    ::glLineWidth(1);
    glTranslated(0,0,-0.2);
    for(unsigned int iface=0;iface<aFace.size();++iface){
      const CCad2D_FaceGeo& face = aFace[iface];
      DrawMeshTri2D_Face(face.aTri, face.aXY);
      DrawMeshTri2D_Edge(face.aTri, face.aXY);
    }
    glTranslated(0,0,+0.2);
  }
}

void CCad2D::Mouse(int btn, int action, int mods,
                   const std::vector<double>& src,
                   const std::vector<double>& dir,
                   double view_height)
{
//  std::cout << "mouse called btn:" << btn << "   action:" << action << " mods:" << mods << std::endl;
//  std::cout << src[0] << " " << src[1] << " " << src[2] << std::endl;
//  std::cout << dir[0] << " " << dir[1] << " " << dir[2] << std::endl;
  if( action == 1 ){ //press
    const double x0 = src[0];
    const double y0 = src[1];
    this->ivtx_picked = -1;
    for(unsigned int ivtx=0;ivtx<aVtx.size();++ivtx){
      double x1 = aVtx[ivtx].pos.x;
      double y1 = aVtx[ivtx].pos.y;
      double dist = sqrt( (x0-x1)*(x0-x1)+(y0-y1)*(y0-y1) );
      if( dist < view_height*0.05 ){
        this->ivtx_picked = ivtx;
        return;
      }
    }
  }
  else{
    ivtx_picked = -1;
  }
}

void CCad2D::Motion(const std::vector<double>& src0,
                    const std::vector<double>& src1,
                    const std::vector<double>& dir)
{
//  std::cout << src0[0] << " " << src0[1] << " " << src0[2] << std::endl;
//  std::cout << src1[0] << " " << src1[1] << " " << src1[2] << std::endl;
//  std::cout << dir[0] << " " << dir[1] << " " << dir[2] << std::endl;
  if( ivtx_picked >= 0 && ivtx_picked < (int)aVtx.size() ){
    aVtx[ivtx_picked].pos.x = src1[0];
    aVtx[ivtx_picked].pos.y = src1[1];
  }
  for(unsigned int ie=0;ie<topo.aEdge.size();++ie){
    const CCadTopo::CEdge& e = topo.aEdge[ie];
    if( e.iv0 == ivtx_picked || e.iv1 == ivtx_picked ){
      aEdge[ie].GenMesh(ie,topo,aVtx);
    }
  }
  /*
  for(unsigned int ifc=0;ifc<topo.aFace.size();++ifc){
    const CCadTopo::CFace& fc = topo.aFace[ifc];
  }
   */
  for(unsigned int ifc=0;ifc<topo.aFace.size();++ifc){
    aFace[ifc].GenMesh(ifc, topo, aEdge);
  }
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
 double elen) const
{
  const int iface0 = 0;
  assert( iface0<topo.aFace.size() );
  const std::vector< std::pair<int,bool> >& aIE = topo.aFace[iface0].aIE;
  std::vector< std::vector<double> > aaXY;
  aaXY.resize(1);
  for(unsigned int iie=0;iie<aIE.size();++iie){
    const unsigned int ie0 = (unsigned int)aIE[iie].first;
    assert( ie0<topo.aEdge.size() );
    const bool dir0 = aIE[iie].second;
    //    int iv0 = (dir0) ? topo.aEdge[ie0].iv0 : topo.aEdge[ie0].iv1;
    {
      const CCad2D_EdgeGeo& eg0 = aEdge[ie0];
      CVector2 p0 = (dir0) ? eg0.p0 : eg0.p1;
      aaXY[0].push_back(p0.x);
      aaXY[0].push_back(p0.y);
    }
  }
  std::vector<int> loopIP_ind,loopIP;
  std::vector<CVector2> aVec2;
  {
    JArray_FromVecVec_XY(loopIP_ind,loopIP, aVec2,
                         aaXY);
    if( !CheckInputBoundaryForTriangulation(loopIP_ind,aVec2) ){
      return;
    }
    FixLoopOrientation(loopIP,
                       loopIP_ind,aVec2);
    if( elen > 10e-10 ){
      ResamplingLoop(loopIP_ind,loopIP,aVec2,
                     elen );
    }
  }
  {
    std::vector<CEPo2> aPo2D;
    std::vector<ETri> aETri;
    Meshing_SingleConnectedShape2D(aPo2D, aVec2, aETri,
                                   loopIP_ind,loopIP);
    if( elen > 1.0e-10 ){
      CInputTriangulation_Uniform param(1.0);
      MeshingInside(aPo2D,aETri,aVec2, loopIP,
                    elen, param);
    }
    MeshTri2D_Export(aXY,aTri, aVec2,aETri);
  }
}

void CCad2D::GetPointsEdge
(std::vector<int>& aIdP,
 const double* pXY, int np,
 const std::vector<int>& aIE,
 double tolerance ) const
{
  aIdP.clear();
  for(int ip=0;ip<np;++ip){
//    if( pBC[ip*nDimVal+idimVal] == iflag ){ continue; } // flag already set for this point
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

std::vector<double> CCad2D::GetVertexXY_Face
(int iface) const
{
  std::vector<double> aXY;
  std::vector<int> aIdV = topo.aFace[iface].GetArray_IdVertex(topo.aEdge);
  for(unsigned int iv=0;iv<aIdV.size();++iv){
    aXY.push_back( aVtx[iv].pos.x );
    aXY.push_back( aVtx[iv].pos.y );
  }
  return aXY;
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
//    int iv0 = (dir0) ? topo.aEdge[ie0].iv0 : topo.aEdge[ie0].iv1;
    {
      const CCad2D_EdgeGeo& eg0 = aEdgeGeo[ie0];
      CVector2 p0 = (dir0) ? eg0.p0 : eg0.p1;
      aaXY[0].push_back(p0.x);
      aaXY[0].push_back(p0.y);
    }
  }
  std::vector<int> loopIP_ind,loopIP;
  std::vector<CVector2> aVec2;
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
    std::vector<CEPo2> aPo2D;
    std::vector<ETri> aETri;
    Meshing_SingleConnectedShape2D(aPo2D, aVec2, aETri,
                                   loopIP_ind, loopIP);
    MeshTri2D_Export(aXY,aTri, aVec2,aETri);
  }
}

