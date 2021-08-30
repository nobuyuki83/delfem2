/*
 * Copyright (c) 2019 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */

#include "delfem2/opengl/old/cad2dtriv2.h"

#include <climits>

#if defined(_WIN32) // windows
  #define NOMINMAX   // to remove min,max macro
  #include <windows.h>
#endif

#if defined(__APPLE__) && defined(__MACH__)
  #define GL_SILENCE_DEPRECATION
  #include <OpenGL/gl.h>
#else
  #include <GL/gl.h>
#endif

#include "delfem2/opengl/old/v2.h"
#include "delfem2/opengl/old/funcs.h"

// -------------------------------------------------

DFM2_INLINE void delfem2::opengl::DrawMeshDynTri_FaceNorm(
    const std::vector<CDynTri>& aSTri,
    const std::vector<CVec2d>& aVec2)
{
  ::glDisable(GL_LIGHTING);
  ::glBegin(GL_TRIANGLES);
  for (const auto & itri : aSTri){
    const unsigned int i0 = itri.v[0];
    const unsigned int i1 = itri.v[1];
    const unsigned int i2 = itri.v[2];
    if( i0 == UINT_MAX ){
      assert( i1 == UINT_MAX );
      assert( i2 == UINT_MAX );
      continue;
    }
    const CVec2d& p0 = aVec2[i0];
    const CVec2d& p1 = aVec2[i1];
    const CVec2d& p2 = aVec2[i2];
    ::glVertex2d(p0.x,p0.y);
    ::glVertex2d(p1.x,p1.y);
    ::glVertex2d(p2.x,p2.y);
  }
  ::glEnd();
}

DFM2_INLINE void delfem2::opengl::DrawMeshDynTri_Edge(
    const std::vector<CDynTri>& aSTri,
    const std::vector<CVec2d>& aVec2)
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
    const CVec2d& p0 = aVec2[i0];
    const CVec2d& p1 = aVec2[i1];
    const CVec2d& p2 = aVec2[i2];
    glVertex2d(p0.x,p0.y);  glVertex2d(p1.x,p1.y);
    glVertex2d(p1.x,p1.y);  glVertex2d(p2.x,p2.y);
    glVertex2d(p2.x,p2.y);  glVertex2d(p0.x,p0.y);
  }
  ::glEnd();
}

// -------------------------------------------------------------------------

DFM2_INLINE void delfem2::opengl::Draw_CCad2DEdge(
    const CCad2D_EdgeGeo& edge,
    bool is_selected,
    int ipicked_elem)
{
  if( is_selected ){ ::glColor3d(1,1,0); }
  else{ ::glColor3d(0,0,0); }
  ::glBegin(GL_LINE_STRIP);
  myGlVertex( edge.p0 );
  for(const auto & ip : edge.aP){
    myGlVertex( ip );
  }
  myGlVertex( edge.p1 );
  ::glEnd();
  //
  if( is_selected ){
    if( edge.type_edge == CCad2D_EdgeGeo::BEZIER_CUBIC ){
      assert( edge.param.size() == 4 );
      const CVec2d lx = (edge.p1 - edge.p0);
      const CVec2d ly = CVec2d(lx.y,-lx.x);
      const CVec2d q0 = edge.p0 + edge.param[0]*lx + edge.param[1]*ly;
      const CVec2d q1 = edge.p0 + edge.param[2]*lx + edge.param[3]*ly;
      ::glColor3d(0,1,0);
      ::glBegin(GL_LINES);
      myGlVertex(edge.p0);
      myGlVertex(q0);
      myGlVertex(edge.p1);
      myGlVertex(q1);
      ::glEnd();
      if( ipicked_elem == 1 ){ ::glColor3d(0.8, 0.0, 0.0 ); }
      else{ ::glColor3d(0.0, 0.8, 0.0 ); }
      ::glBegin(GL_POINTS);
      myGlVertex(q0);
      ::glEnd();
      if( ipicked_elem == 2 ){ ::glColor3d(0.8, 0.0, 0.0 ); }
      else{ ::glColor3d(0.0, 0.8, 0.0 ); }
      ::glBegin(GL_POINTS);
      myGlVertex(q1);
      ::glEnd();
    }
    else if( edge.type_edge == CCad2D_EdgeGeo::BEZIER_QUADRATIC ){
      assert( edge.param.size() == 2 );
      const CVec2d lx = (edge.p1 - edge.p0);
      const CVec2d ly = CVec2d(lx.y,-lx.x);
      const CVec2d q0 = edge.p0 + edge.param[0]*lx + edge.param[1]*ly;
      ::glColor3d(0,1,0);
      ::glBegin(GL_LINES);
      myGlVertex(edge.p0);
      myGlVertex(q0);
      myGlVertex(edge.p1);
      myGlVertex(q0);
      ::glEnd();
      if( ipicked_elem == 1 ){ ::glColor3d(0.8, 0.0, 0.0 ); }
      else{ ::glColor3d(0.0, 0.8, 0.0 ); }
      ::glBegin(GL_POINTS);
      myGlVertex(q0);
      ::glEnd();
    }
  }
}

DFM2_INLINE void delfem2::opengl::Draw_CCad2D(const CCad2D& cad2d)
{
  const std::vector<CCad2D_VtxGeo>& aVtx = cad2d.aVtx;
  const std::vector<CCad2D_EdgeGeo>& aEdge = cad2d.aEdge;
  const std::vector<CCad2D_FaceGeo>& aFace = cad2d.aFace;
  int ivtx_picked = cad2d.ivtx_picked;
  int iedge_picked = cad2d.iedge_picked;
  int iface_picked = cad2d.iface_picked;
  int ipicked_elem = cad2d.ipicked_elem;
  bool is_draw_face = cad2d.is_draw_face;
  ///
  ::glDisable(GL_LIGHTING);
  ::glDisable(GL_TEXTURE_2D);
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
      const CCad2D_FaceGeo& face = aFace[iface];
      if( (int)iface == iface_picked ){ ::glColor3d(1,1,0); }
      else{ ::glColor3d(0.8,0.8,0.8); }
      Draw_MeshTri(cad2d.aVec2_Tessellation, face.aTri);
      ::glColor3d(0.0,0.0,0.0);
      Draw_MeshTri_Edge(cad2d.aVec2_Tessellation, face.aTri);
    }
    glTranslated(0,0,+0.2);
  }
}

