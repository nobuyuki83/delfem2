/*
 * Copyright (c) 2019 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */

#include <cmath>
#include <algorithm>
#include "delfem2/cad3d.h"

// ----------
#include <GLFW/glfw3.h>
#include "delfem2/opengl/funcs_glold.h"
#include "delfem2/opengl/v3q_glold.h"
#include "delfem2/opengl/glfw/viewer_glfw.h"

#ifndef M_PI
#  define M_PI 3.141592653589793
#endif

namespace dfm2 = delfem2;

// ------------------------------

/*




void CCad3D::DrawFace_LeftRight() const
{
  ::glEnable(GL_LIGHTING);
  dfm2::opengl::myGlColorDiffuse(color_face);
  dfm2::opengl::DrawMeshTri3D_FaceNorm(aXYZ, aTri, aNorm);
}

void CCad3D_Face::DrawBackFace() const
{
  ::glEnable(GL_LIGHTING);
  dfm2::opengl::myGlColorDiffuse(dfm2::CColor::Gray(0.9));
  dfm2::opengl::DrawMeshTri3D_FaceNorm_XYsym(aXYZ, aTri);
}


*/

void DrawFace(const dfm2::CCad3D_Face& cf)
{
  dfm2::opengl::DrawMeshTri3D_FaceNorm(cf.aXYZ, cf.aTri, cf.aNorm);
}

void DrawEdge(const dfm2::CCad3D_Face& cf)
{
  ::glLineWidth(1);
  dfm2::opengl::DrawMeshTri3D_Edge(cf.aXYZ, cf.aTri);
}

void DrawFace_RightSelected(
    const dfm2::CCad3D& cad,
    bool is_edge)
{
  {
    float specular[4] ={ 0.797357f, 0.723991f, 0.208006f, 1.0f};
    float shine =83.2f ;
    glMaterialfv(GL_FRONT, GL_SPECULAR, specular);
    glMaterialf(GL_FRONT, GL_SHININESS, shine);
  }
  glEnable(GL_NORMALIZE);
  
  ::glEnable(GL_LIGHTING);
  for(std::size_t iface=0;iface<cad.aFace.size();++iface){
    if( iface == cad.iface_picked ){
      dfm2::opengl::myGlColorDiffuse(cad.color_face_selected);
    }
    else{
      dfm2::opengl::myGlColorDiffuse(cad.color_face);
    }
    DrawFace(cad.aFace[iface]);
    if( is_edge ){ DrawEdge(cad.aFace[iface]); }
  }
}


void DrawHandler
 (const dfm2::CCad3D_Edge& e,
  int ielem_picked, double view_height)
{
  ::glColor3d(0,1,0);
  dfm2::opengl::DrawCylinder(e.p0, e.q0, view_height*0.01);
  dfm2::opengl::DrawCylinder(e.p1, e.q1, view_height*0.01);
  {
    if( ielem_picked == 2 ){ ::glColor3d(1,0.0,0); }
    else{                    ::glColor3d(0,0.9,0); }
    ::glMatrixMode(GL_MODELVIEW);
    ::glPushMatrix();
    dfm2::opengl::myGlTranslate(e.q0);
    ::glScaled(view_height*0.02, view_height*0.02, view_height*0.02);
    dfm2::opengl::DrawSphere(32, 32);
    ::glPopMatrix();
  }
  {
    if( ielem_picked == 3 ){ ::glColor3d(1,0.0,0); }
    else{                    ::glColor3d(0,0.9,0); }
    ::glMatrixMode(GL_MODELVIEW);
    ::glPushMatrix();
    dfm2::opengl::myGlTranslate(e.q1);
    ::glScaled(view_height*0.02, view_height*0.02, view_height*0.02);
    dfm2::opengl::DrawSphere(32, 32);
    ::glPopMatrix();
  }
  ::glEnd();
}

void DrawLine
 (const dfm2::CCad3D_Edge& ed,
  bool is_picked, double view_height)
{
  ::glDisable(GL_LIGHTING);
  if( is_picked ){
    ::glColor3d(1,1,0);
  }
  else if( ed.is_sim ){
    ::glColor3d(0,0,0);
  }
  else{
    ::glColor3d(0.8,0.8,0.8);
  }
  ::glLineWidth(2);
  ::glBegin(GL_LINE_STRIP);
  for(const auto & ip : ed.aP){ dfm2::opengl::myGlVertex(ip); }
  ::glEnd();
}

void Draw
 (const dfm2::CCad3D_Vertex& vtx,
  bool is_selected, int ielem, double view_height)
{
  ::glDisable(GL_LIGHTING);
  if( is_selected ){
    double s = view_height*0.3;
    ::glDisable(GL_LIGHTING);
    if( !vtx.isConst[0] ){
      ::glColor3d(1, 0, 0);
      dfm2::opengl::DrawArrow(vtx.pos, dfm2::CVec3d(+s, 0, 0));
      dfm2::opengl::DrawArrow(vtx.pos, dfm2::CVec3d(-s, 0, 0));
    }
    if( !vtx.isConst[1] ){
      ::glColor3d(0, 1, 0);
      dfm2::opengl::DrawArrow(vtx.pos, dfm2::CVec3d(0, +s, 0));
      dfm2::opengl::DrawArrow(vtx.pos, dfm2::CVec3d(0, -s, 0));
    }
    if( !vtx.isConst[2] ){
      ::glColor3d(0, 0, 1);
      dfm2::opengl::DrawArrow(vtx.pos, dfm2::CVec3d(0, 0, +s));
      dfm2::opengl::DrawArrow(vtx.pos, dfm2::CVec3d(0, 0, -s));
    }
  }
  else{
    ::glDisable(GL_LIGHTING);
    ::glColor3d(0,0,1);
    ::glMatrixMode(GL_MODELVIEW);
    ::glPushMatrix();
    dfm2::opengl::myGlTranslate(vtx.pos);
    ::glDisable(GL_CULL_FACE);
    ::glScaled(view_height*0.03, view_height*0.03, view_height*0.03);
    dfm2::opengl::DrawSphere(32, 32);
    ::glPopMatrix();
  }
    ////
  dfm2::opengl::myGlColorDiffuse(dfm2::CColor::Red());
    //  DrawArrow(pos, +norm*view_height*0.1);
    //  ::DrawArrow(pos, -norm*0.15);
}

void DrawVtxEdgeHandler
 (const dfm2::CCad3D& cad,
  double view_height)
{
  {
    const std::vector<dfm2::CCad3D_Vertex>& aVertex = cad.aVertex;
    for(std::size_t icp=0;icp<aVertex.size();++icp){
      Draw(aVertex[icp], icp==cad.ivtx_picked, cad.ielem_vtx_picked, view_height);
    }
  }
  
  glEnable(GL_BLEND);
  glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
  
  glEnable(GL_LINE_SMOOTH);
  glHint(GL_LINE_SMOOTH_HINT, GL_NICEST);
    //  glEnable(GL_MULTISAMPLE);
  
  {
    const std::vector<dfm2::CCad3D_Edge>& aEdge = cad.aEdge;
    dfm2::opengl::myGlColorDiffuse(dfm2::CColor::Blue());
    for(unsigned int ie=0;ie<aEdge.size();++ie){
      bool is_loop0 = std::find(cad.aIE_picked.begin(), cad.aIE_picked.end(), std::make_pair(ie,true ) ) != cad.aIE_picked.end();
      bool is_loop1 = std::find(cad.aIE_picked.begin(), cad.aIE_picked.end(), std::make_pair(ie,false) ) != cad.aIE_picked.end();
      bool is_loop = is_loop0 || is_loop1;
      DrawLine(cad.aEdge[ie], is_loop, view_height);
    }
    if( cad.iedge_picked != -1 ){
      DrawHandler(aEdge[cad.iedge_picked], cad.ielem_edge_picked, view_height);
    }
  }
  
  ::glDisable(GL_LIGHTING);
  ::glDisable(GL_CULL_FACE);
  if( cad.plane_inorm >= 0 && cad.plane_inorm < 3 ){
    dfm2::CVec3d plane_org = cad.plane_org;
    double plane_sizeX = cad.plane_sizeX;
    double plane_sizeY = cad.plane_sizeY;
    dfm2::opengl::myGlColorDiffuse(dfm2::CColor(0.0, 0.0, 1.0, 0.5));
    dfm2::CVec3d plane_ex = dfm2::CVec3d::Axis((cad.plane_inorm+1)%3);
    dfm2::CVec3d plane_ey = dfm2::CVec3d::Axis((cad.plane_inorm+2)%3);
    dfm2::CVec3d p0 = plane_org-plane_sizeX*plane_ex-plane_sizeY*plane_ey;
    dfm2::CVec3d p1 = plane_org+plane_sizeX*plane_ex-plane_sizeY*plane_ey;
    dfm2::CVec3d p2 = plane_org+plane_sizeX*plane_ex+plane_sizeY*plane_ey;
    dfm2::CVec3d p3 = plane_org-plane_sizeX*plane_ex+plane_sizeY*plane_ey;
    ::glBegin(GL_QUADS);
    dfm2::opengl::myGlVertex(p0);
    dfm2::opengl::myGlVertex(p1);
    dfm2::opengl::myGlVertex(p2);
    dfm2::opengl::myGlVertex(p3);
    ::glEnd();
      ////
    ::glLineWidth(5);
    if( cad.ielem_edge_picked == 1 ){ dfm2::opengl::myGlColorDiffuse(dfm2::CColor(1.f, 0.f, 0.f, 0.9f)); }
    else{                             dfm2::opengl::myGlColorDiffuse(dfm2::CColor(0.f, 0.f, 1.f, 0.9f)); }
    dfm2::opengl::DrawCylinder(p0, p1, view_height*0.01);
    dfm2::opengl::DrawCylinder(p1, p2, view_height*0.01);
    dfm2::opengl::DrawCylinder(p2, p3, view_height*0.01);
    dfm2::opengl::DrawCylinder(p3, p0, view_height*0.01);
  }
  
  {
    ::glDisable(GL_DEPTH_TEST);
    ::glMatrixMode(GL_PROJECTION);
    ::glPushMatrix();
    ::glLoadIdentity();
    ::glMatrixMode(GL_MODELVIEW);
    ::glPushMatrix();
    ::glLoadIdentity();
    //
    /*
    ::glDisable(GL_LIGHTING);
    ::glColor3d(1,0,0);
    ::glLineWidth(3);
    ::glBegin(GL_LINE_STRIP);
    for(const auto & ist : cad.aStroke){
      dfm2::opengl::myGlVertex(ist);
    }
    ::glEnd();
     */
    //
    ::glMatrixMode(GL_PROJECTION);
    ::glPopMatrix();
    ::glMatrixMode(GL_MODELVIEW);
    ::glPopMatrix();
    ::glEnable(GL_DEPTH_TEST);
  }
  
  glDisable(GL_LINE_SMOOTH);
  
}

// ------------------------------

int main(int argc,char* argv[])
{
  class CCAD3DViewer : public delfem2::opengl::CViewer_GLFW {
  public:
    CCAD3DViewer(){
      cad.Initialize_Sphere();
      camera.view_height = 1.5;
      camera.camera_rot_mode = delfem2::CCam3_OnAxisZplusLookOrigin<double>::CAMERA_ROT_MODE::YTOP;
      cad.imode_edit = dfm2::CCad3D::EDIT_MOVE;
    }
    void Draw(){
      this->DrawBegin_oldGL();
      DrawFace_RightSelected(cad,false);
      DrawVtxEdgeHandler(cad,camera.view_height);
      this->SwapBuffers();
    }
    void mouse_press(const float src[3], const float dir[3]) override{
      const dfm2::CVec3d src_pick(src), dir_pick(dir);
      int nw, nh; glfwGetFramebufferSize(this->window,&nw,&nh);
      const float asp = (float)nw/nh;
      float mMV[16], mPrj[16]; camera.Mat4_MVP_OpenGL(mMV, mPrj, asp);
      cad.MouseDown(src_pick, dir_pick,
                    dfm2::CVec2d(nav.mouse_x,nav.mouse_y),
                    mMV,mPrj,
                    camera.view_height);
    }
    void mouse_drag(const float src0[3], const float src1[3], const float dir[3]) override{
      dfm2::CVec2d sp0(nav.mouse_x-nav.dx, nav.mouse_y-nav.dy);
      dfm2::CVec2d sp1(nav.mouse_x, nav.mouse_y);
      const dfm2::CVec3d src_pick(src1), dir_pick(dir);
      int nw, nh; glfwGetFramebufferSize(this->window,&nw,&nh);
      const float asp = (float)nw/nh;
      float mMV[16], mPrj[16]; camera.Mat4_MVP_OpenGL(mMV, mPrj, asp);
      cad.MouseMotion(src_pick,dir_pick, sp0,sp1, mMV, mPrj);
    }
  public:
    dfm2::CCad3D cad;
  };
  // -------------
  CCAD3DViewer viewer;
  viewer.Init_oldGL();
  delfem2::opengl::setSomeLighting();
  while(!glfwWindowShouldClose(viewer.window)){
    viewer.Draw();
    glfwPollEvents();
  }
  glfwDestroyWindow(viewer.window);
  glfwTerminate();
  exit(EXIT_SUCCESS);
}


