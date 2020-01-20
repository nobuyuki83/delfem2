/*
 * Copyright (c) 2019 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */


#if defined(__APPLE__) && defined(__MACH__)
  #include <OpenGL/gl.h>
#elif defined(_WIN32) // windows
  #include <windows.h>
  #include <GL/gl.h>
#else
  #include <GL/gl.h>
#endif

#include "delfem2/opengl/glold_v23.h"
#include "delfem2/opengl/glold_v23dtricad.h"
#include "delfem2/opengl/glold_funcs.h"

namespace dfm2 = delfem2;

// -------------------------------------------------

void dfm2::opengl::DrawMeshDynTri_FaceNorm
(const std::vector<CDynTri>& aSTri,
 const std::vector<CVector2>& aVec2)
{
  ::glDisable(GL_LIGHTING);
  ::glBegin(GL_TRIANGLES);
  for (const auto & itri : aSTri){
    const int i0 = itri.v[0];
    const int i1 = itri.v[1];
    const int i2 = itri.v[2];
    if( i0 == -1 ){
      assert( i1 == -1 );
      assert( i2 == -1 );
      continue;
    }
    const CVector2& p0 = aVec2[i0];
    const CVector2& p1 = aVec2[i1];
    const CVector2& p2 = aVec2[i2];
    ::glVertex2d(p0.x,p0.y);
    ::glVertex2d(p1.x,p1.y);
    ::glVertex2d(p2.x,p2.y);
  }
  ::glEnd();
}

void dfm2::opengl::DrawMeshDynTri_Edge
(const std::vector<CDynTri>& aSTri,
 const std::vector<CVector2>& aVec2)
{
  ::glDisable(GL_LIGHTING);
  ::glLineWidth(1);
  ::glColor3d(0,0,0);
  ::glBegin(GL_LINES);
  for (const auto & itri : aSTri){
    const int i0 = itri.v[0];
    const int i1 = itri.v[1];
    const int i2 = itri.v[2];
    if( i0 == -1 ){
      assert( i1 == -1 );
      assert( i2 == -1 );
    }
    const CVector2& p0 = aVec2[i0];
    const CVector2& p1 = aVec2[i1];
    const CVector2& p2 = aVec2[i2];
    glVertex2d(p0.x,p0.y);  glVertex2d(p1.x,p1.y);
    glVertex2d(p1.x,p1.y);  glVertex2d(p2.x,p2.y);
    glVertex2d(p2.x,p2.y);  glVertex2d(p0.x,p0.y);
  }
  ::glEnd();
}


void dfm2::opengl::DrawMeshDynTri_FaceNorm
(const std::vector<CDynTri>& aSTri,
 const std::vector<CVector3>& aVec3)
{
  //  ::glPushAttrib(GL_ENABLE_BIT);
  
  ::glEnable(GL_LIGHTING);
  //  ::myGlColorDiffuse(CColor::Orange());
  ::glBegin(GL_TRIANGLES);
  for (const auto & itri : aSTri){
    const int i0 = itri.v[0];
    const int i1 = itri.v[1];
    const int i2 = itri.v[2];
    if( i0 == -1 ){
      assert( i1 == -1 );
      assert( i2 == -1 );
      continue;
    }
    {
      CVector3 n; UnitNormal(n, aVec3[i0], aVec3[i1], aVec3[i2]);
      ::glNormal3d(n.x(),n.y(),n.z());
    }
    {
      CVector3 p0 = aVec3[i0];
      ::glVertex3d(p0.x(),p0.y(),p0.z());
    }
    {
      CVector3 p1 = aVec3[i1];
      ::glVertex3d(p1.x(),p1.y(),p1.z());
    }
    {
      CVector3 p2 = aVec3[i2];
      ::glVertex3d(p2.x(),p2.y(),p2.z());
    }
  }
  ::glEnd();
}

void dfm2::opengl::DrawMeshDynTri_Edge
(const std::vector<CDynTri>& aSTri,
 const std::vector<CVector3>& aVec3)
{
  ::glDisable(GL_LIGHTING);
  ::glLineWidth(1);
  ::glColor3d(0,0,0);
  ::glBegin(GL_LINES);
  for (const auto & itri : aSTri){
    const int i0 = itri.v[0];
    const int i1 = itri.v[1];
    const int i2 = itri.v[2];
    if( i0 == -1 ){
      assert( i1 == -1 );
      assert( i2 == -1 );
    }
    const CVector3& p0 = aVec3[i0];
    const CVector3& p1 = aVec3[i1];
    const CVector3& p2 = aVec3[i2];
    glVertex3d(p0.x(),p0.y(),p0.z());
    glVertex3d(p1.x(),p1.y(),p1.z());
    glVertex3d(p1.x(),p1.y(),p1.z());
    glVertex3d(p2.x(),p2.y(),p2.z());
    glVertex3d(p2.x(),p2.y(),p2.z());
    glVertex3d(p0.x(),p0.y(),p0.z());
  }
  ::glEnd();
}


void dfm2::opengl::DrawMeshDynTri3D_Edge
(const std::vector<double>& aXYZ,
 const std::vector<CDynTri>& aSTri)
{
  ::glDisable(GL_LIGHTING);
  ::glLineWidth(1);
  ::glColor3d(0,0,0);
  ::glBegin(GL_LINES);
  for (const auto & itri : aSTri){
    const int i0 = itri.v[0];
    const int i1 = itri.v[1];
    const int i2 = itri.v[2];
    if( i0 == -1 ){
      assert( i1 == -1 );
      assert( i2 == -1 );
    }
    glVertex3d(aXYZ[i0*3+0],aXYZ[i0*3+1],aXYZ[i0*3+2]);
    glVertex3d(aXYZ[i1*3+0],aXYZ[i1*3+1],aXYZ[i1*3+2]);
    glVertex3d(aXYZ[i1*3+0],aXYZ[i1*3+1],aXYZ[i1*3+2]);
    glVertex3d(aXYZ[i2*3+0],aXYZ[i2*3+1],aXYZ[i2*3+2]);
    glVertex3d(aXYZ[i2*3+0],aXYZ[i2*3+1],aXYZ[i2*3+2]);
    glVertex3d(aXYZ[i0*3+0],aXYZ[i0*3+1],aXYZ[i0*3+2]);
  }
  ::glEnd();
}

// -------------------------------------------------------------------------

void dfm2::opengl::Draw_CCad2DEdge
 (const dfm2::CCad2D_EdgeGeo& edge,
  bool is_selected,
  int ipicked_elem)
{
  if( is_selected ){ ::glColor3d(1,1,0); }
  else{ ::glColor3d(0,0,0); }
  ::glBegin(GL_LINE_STRIP);
  dfm2::opengl::myGlVertex( edge.p0 );
  for(const auto & ip : edge.aP){
    dfm2::opengl::myGlVertex( ip );
  }
  dfm2::opengl::myGlVertex( edge.p1 );
  ::glEnd();
  //
  if( is_selected ){
    if( edge.type_edge == 1 ){
      assert( edge.param.size() == 4 );
      const CVector2 lx = (edge.p1 - edge.p0).Normalize();
      const CVector2 ly = CVector2(lx.y,-lx.x);
      const CVector2 q0 = edge.p0 + edge.param[0]*lx + edge.param[1]*ly;
      const CVector2 q1 = edge.p1 + edge.param[2]*lx + edge.param[3]*ly;
      ::glColor3d(0,1,0);
      ::glBegin(GL_LINES);
      dfm2::opengl::myGlVertex(edge.p0);
      dfm2::opengl::myGlVertex(q0);
      dfm2::opengl::myGlVertex(edge.p1);
      dfm2::opengl::myGlVertex(q1);
      ::glEnd();
      if( ipicked_elem == 1 ){ ::glColor3d(0.8, 0.0, 0.0 ); }
      else{ ::glColor3d(0.0, 0.8, 0.0 ); }
      ::glBegin(GL_POINTS);
      dfm2::opengl::myGlVertex(q0);
      ::glEnd();
      if( ipicked_elem == 2 ){ ::glColor3d(0.8, 0.0, 0.0 ); }
      else{ ::glColor3d(0.0, 0.8, 0.0 ); }
      ::glBegin(GL_POINTS);
      dfm2::opengl::myGlVertex(q1);
      ::glEnd();
    }
  }
}

void dfm2::opengl::Draw_CCad2D(const dfm2::CCad2D& cad2d)
{
  const std::vector<dfm2::CCad2D_VtxGeo>& aVtx = cad2d.aVtx;
  const std::vector<dfm2::CCad2D_EdgeGeo>& aEdge = cad2d.aEdge;
  const std::vector<dfm2::CCad2D_FaceGeo>& aFace = cad2d.aFace;
  int ivtx_picked = cad2d.ivtx_picked;
  int iedge_picked = cad2d.iedge_picked;
  int iface_picked = cad2d.iface_picked;
  int ipicked_elem = cad2d.ipicked_elem;
  bool is_draw_face = cad2d.is_draw_face;
  ///
  ::glDisable(GL_LIGHTING);
  ::glPointSize(6);
  ::glBegin(GL_POINTS);
  for(size_t iv=0;iv<aVtx.size();++iv){
    if( (int)iv == ivtx_picked ){ ::glColor3d(1,1,0); }
    else{ ::glColor3d(1,0,0); }
    ::glVertex3d( aVtx[iv].pos.x, aVtx[iv].pos.y, 0.0);
  }
  ::glEnd();
  //
  ::glLineWidth(3);
  for(size_t ie=0;ie<aEdge.size();++ie){
    Draw_CCad2DEdge(aEdge[ie],
                    (int)ie == iedge_picked,
                    ipicked_elem);
  }
  //
  if( is_draw_face ){
    ::glLineWidth(1);
    glTranslated(0,0,-0.2);
    for(size_t iface=0;iface<aFace.size();++iface){
      const dfm2::CCad2D_FaceGeo& face = aFace[iface];
      if( (int)iface == iface_picked ){ ::glColor3d(1,1,0); }
      else{ ::glColor3d(0.8,0.8,0.8); }
      dfm2::opengl::Draw_MeshTri(cad2d.aVec2_Tessellation, face.aTri);
      ::glColor3d(0.0,0.0,0.0);
      dfm2::opengl::Draw_MeshTri_Edge(cad2d.aVec2_Tessellation, face.aTri);
    }
    glTranslated(0,0,+0.2);
  }
}
