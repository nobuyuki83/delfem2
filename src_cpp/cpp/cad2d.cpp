#include <stdio.h>
#include <deque>
#include <set>

#if defined(__APPLE__) && defined(__MACH__)
#include <GLUT/glut.h>
#else
#include <GL/glut.h>
#endif

#include "delfem2/mshtopo.h"
#include "delfem2/msh.h"

#include "delfem2/dyntri_v3.h"

#include "delfem2/funcs_gl.h"
#include "delfem2/color_gl.h"
#include "delfem2/v23_gl.h"

#include "delfem2/cad2d.h"


void CCad2D::Draw() const
{
  ::glDisable(GL_LIGHTING);
  ::glColor3d(1,0,0);
  ::glPointSize(6);
  ::glBegin(GL_POINTS);
  for(int iv=0;iv<aVtx.size();++iv){
    ::glVertex3d( aVtx[iv].pos.x, aVtx[iv].pos.y, 0.0);
  }
  ::glEnd();
  /////
  ::glColor3d(0,0,0);
  ::glLineWidth(3);
  ::glBegin(GL_LINES);
  for(int ie=0;ie<aEdge.size();++ie){
    int iv0 = topo.aEdge[ie].iv0;
    int iv1 = topo.aEdge[ie].iv1;
    ::glVertex3d( aVtx[iv0].pos.x, aVtx[iv0].pos.y, -0.1);
    ::glVertex3d( aVtx[iv1].pos.x, aVtx[iv1].pos.y, -0.1);
  }
  ::glEnd();
  //////
  ::glColor3d(0.8,0.8,0.8);
  ::glLineWidth(1);
  glTranslated(0,0,-0.2);
  for(int iface=0;iface<aFace.size();++iface){
    const CCad2D_FaceGeo& face = aFace[iface];
    DrawMeshTri2D_Face(face.aTri, face.aXY);
    DrawMeshTri2D_Edge(face.aTri, face.aXY);
  }
  glTranslated(0,0,+0.2);
}

void CCad2D::Initialize_Square(){
  Clear();
  topo.AddPolygon(4);
  aVtx.push_back(CVector2(-1.0, -1.0));
  aVtx.push_back(CVector2(+1.0, -1.0));
  aVtx.push_back(CVector2(+1.0, +1.0));
  aVtx.push_back(CVector2(-1.0, +1.0));
  int iedge0 = aEdge.size();
  aEdge.push_back(CCad2D_EdgeGeo());
  aEdge.push_back(CCad2D_EdgeGeo());
  aEdge.push_back(CCad2D_EdgeGeo());
  aEdge.push_back(CCad2D_EdgeGeo());
  int iface0 = aFace.size();
  aFace.push_back(CCad2D_FaceGeo());
  ////
  aEdge[iedge0+0].GenMesh(iedge0+0,topo,aVtx);
  aEdge[iedge0+1].GenMesh(iedge0+1,topo,aVtx);
  aEdge[iedge0+2].GenMesh(iedge0+2,topo,aVtx);
  aEdge[iedge0+3].GenMesh(iedge0+3,topo,aVtx);
  aFace[iface0].GenMesh(topo,iface0,aEdge);
}

///////////////////////////////////////////////////////////

void CCad2D_EdgeGeo::GenMesh
(int iedge, const CCadTopo& topo,
 std::vector<CCad2D_VtxGeo>& aVtxGeo)
{
  assert( iedge>=0 && iedge<topo.aEdge.size() );
  const int iv0 = topo.aEdge[iedge].iv0;
  const int iv1 = topo.aEdge[iedge].iv1;
  this->p0 = aVtxGeo[iv0].pos;
  this->p1 = aVtxGeo[iv1].pos;
}

///////////////////////////////////////////////////////////

void CCad2D_FaceGeo::GenMesh
(const CCadTopo& topo,int iface0,
 std::vector<CCad2D_EdgeGeo>& aEdgeGeo)
{
  assert( iface0>=0 && iface0<topo.aFace.size() );
  const std::vector< std::pair<int,bool> >& aIE = topo.aFace[iface0].aIE;
  std::vector<double> aXY_corner;
  for(int iie=0;iie<aIE.size();++iie){
    int ie0 = aIE[iie].first;
    assert( ie0>=0 && ie0<topo.aEdge.size() );
    const bool dir0 = aIE[iie].second;
    int iv0 = (dir0) ? topo.aEdge[ie0].iv0 : topo.aEdge[ie0].iv1;
    {
      const CCad2D_EdgeGeo& eg0 = aEdgeGeo[ie0];
      CVector2 p0 = (dir0) ? eg0.p0 : eg0.p1;
      aXY_corner.push_back(p0.x);
      aXY_corner.push_back(p0.y);
    }
  }
  /*
   for(int ixy=0;ixy<aXY_corner.size()/2;++ixy){
   std::cout << aXY_corner[ixy*2+0] << " " << aXY_corner[ixy*2+1] << std::endl;
   }
   */
  {
    std::vector<int> aPtrVtxInd,aVtxInd;
    std::vector< std::vector<double> > aVecAry0;
    aVecAry0.push_back(aXY_corner);
    GenerateTesselation2(aTri, aXY,  aPtrVtxInd,aVtxInd,
                         -1, false, aVecAry0);
  }
  //      std::cout << aTri.size() << std::endl;
}

